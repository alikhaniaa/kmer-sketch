//! Command-line interface. Each subcommand parses its arguments with clap
//! and then calls into the library modules.

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use kmer_sketch::cluster::{
    label_from_name, make_report, read_matrix, single_linkage, write_report, Linkage,
};
use kmer_sketch::compare::{compare_all, write_matrix, Metric};
use kmer_sketch::count::{count_files, write_counts};
use kmer_sketch::seqio::sample_name;
use kmer_sketch::simulate::{self, Genome, SimConfig, LABEL_SEPARATOR};
use kmer_sketch::sketch::Sketch;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "kmer-sketch",
    version,
    about = "Canonical k-mer counting and FracMinHash sketches for bacterial reads"
)]
struct Cli {
    /// Number of worker threads (0 = use all cores)
    #[arg(short = 't', long, default_value_t = 0, global = true)]
    threads: usize,

    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Count canonical k-mers in FASTA/FASTQ files (plain or .gz)
    Count(CountArgs),
    /// Build one FracMinHash sketch per input file
    Sketch(SketchArgs),
    /// Compare sketches pairwise and write a TSV matrix
    Compare(CompareArgs),
    /// Simulate single-cell read sets (one FASTQ per cell) from genomes
    Simulate(SimulateArgs),
    /// Single-linkage clustering of cells from a comparison matrix
    Cluster(ClusterArgs),
}

#[derive(Args)]
struct CountArgs {
    /// k-mer length (1..=31)
    #[arg(short = 'k', long, default_value_t = 21)]
    ksize: usize,
    /// Only report k-mers seen at least this many times
    #[arg(long, default_value_t = 1)]
    min_count: u64,
    /// Output TSV file (default: standard output)
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Input FASTA/FASTQ files; counts are summed over all files
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
}

#[derive(Args)]
struct SketchArgs {
    /// k-mer length (1..=31)
    #[arg(short = 'k', long, default_value_t = 21)]
    ksize: usize,
    /// Keep about one k-mer hash in `scaled`
    #[arg(short, long, default_value_t = 1000)]
    scaled: u64,
    /// Directory for the .sketch files (one per input, named after the input)
    #[arg(short, long)]
    out_dir: PathBuf,
    /// Input FASTA/FASTQ files
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
}

#[derive(Args)]
struct CompareArgs {
    /// Value to report for each pair
    #[arg(short, long, value_enum, default_value_t = Metric::Jaccard)]
    metric: Metric,
    /// Output TSV file (default: standard output)
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Sketch files written by `sketch`
    #[arg(required = true)]
    sketches: Vec<PathBuf>,
}

#[derive(Args)]
struct SimulateArgs {
    /// Total number of cells; genomes are used in turn
    #[arg(short, long, default_value_t = 40)]
    cells: usize,
    /// Read length in bases
    #[arg(long, default_value_t = 150)]
    read_length: usize,
    /// Smallest per-cell coverage
    #[arg(long, default_value_t = 2.0)]
    min_coverage: f64,
    /// Largest per-cell coverage
    #[arg(long, default_value_t = 30.0)]
    max_coverage: f64,
    /// Per-base substitution error probability
    #[arg(short, long, default_value_t = 0.001)]
    error_rate: f64,
    /// Random seed (same seed and inputs give the same reads)
    #[arg(long, default_value_t = 42)]
    seed: u64,
    /// Write gzip-compressed .fastq.gz files
    #[arg(long)]
    gzip: bool,
    /// Output directory for the per-cell FASTQ files and cells.tsv
    #[arg(short, long)]
    out_dir: PathBuf,
    /// Genome FASTA files; the file name (without extension) is the species label
    #[arg(required = true)]
    genomes: Vec<PathBuf>,
}

#[derive(Args)]
struct ClusterArgs {
    /// Link two cells when their similarity is >= threshold
    /// (or distance <= threshold with --distance)
    #[arg(long)]
    threshold: f64,
    /// The matrix holds distances instead of similarities
    #[arg(long)]
    distance: bool,
    /// True label = part of the cell name before this separator
    #[arg(long, default_value = LABEL_SEPARATOR)]
    label_separator: String,
    /// Also write per-cell assignments (cell, cluster, label) to this TSV file
    #[arg(long)]
    assignments: Option<PathBuf>,
    /// Matrix TSV written by `compare`
    matrix: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .context("cannot set up the thread pool")?;
    }
    match cli.command {
        Command::Count(args) => run_count(args),
        Command::Sketch(args) => run_sketch(args),
        Command::Compare(args) => run_compare(args),
        Command::Simulate(args) => run_simulate(args),
        Command::Cluster(args) => run_cluster(args),
    }
}

/// Open the output file, or standard output when no path is given.
fn open_output(path: Option<&Path>) -> Result<Box<dyn Write>> {
    match path {
        Some(path) => {
            let file =
                File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
            Ok(Box::new(BufWriter::new(file)))
        }
        None => Ok(Box::new(BufWriter::new(io::stdout().lock()))),
    }
}

fn run_count(args: CountArgs) -> Result<()> {
    let counter = count_files(&args.inputs, args.ksize)?;
    let counts = counter.into_sorted_counts(args.min_count);
    let mut out = open_output(args.output.as_deref())?;
    write_counts(&mut out, &counts, args.ksize)?;
    out.flush()?;
    Ok(())
}

fn run_sketch(args: SketchArgs) -> Result<()> {
    // Two inputs with the same sample name would overwrite each other's sketch.
    let mut seen = std::collections::HashSet::new();
    for input in &args.inputs {
        if !seen.insert(sample_name(input)) {
            bail!("two inputs share the sample name '{}'", sample_name(input));
        }
    }
    std::fs::create_dir_all(&args.out_dir)
        .with_context(|| format!("cannot create {}", args.out_dir.display()))?;
    // One sketch per input file, files processed in parallel.
    args.inputs.par_iter().try_for_each(|input| -> Result<()> {
        let name = sample_name(input);
        let sketch = Sketch::from_file(input, &name, args.ksize, args.scaled)?;
        warn_if_empty(&sketch);
        sketch.save(&args.out_dir.join(format!("{name}.sketch")))
    })?;
    eprintln!(
        "wrote {} sketches to {}",
        args.inputs.len(),
        args.out_dir.display()
    );
    Ok(())
}

/// An empty sketch (empty input, or too few k-mers for the chosen `scaled`)
/// is allowed, but every comparison with it gives Jaccard 0 and distance 1.
fn warn_if_empty(sketch: &Sketch) {
    if sketch.hashes.is_empty() {
        eprintln!(
            "warning: sketch '{}' has 0 hashes; its Jaccard is 0 and Mash distance 1 against everything",
            sketch.name
        );
    }
}

fn run_compare(args: CompareArgs) -> Result<()> {
    let sketches: Vec<Sketch> = args
        .sketches
        .par_iter()
        .map(|path| Sketch::load(path))
        .collect::<Result<_>>()?;
    sketches.iter().for_each(warn_if_empty);
    let matrix = compare_all(&sketches, args.metric)?;
    let names: Vec<String> = sketches.iter().map(|s| s.name.clone()).collect();
    let mut out = open_output(args.output.as_deref())?;
    write_matrix(&mut out, &names, &matrix)?;
    out.flush()?;
    Ok(())
}

fn run_simulate(args: SimulateArgs) -> Result<()> {
    let genomes: Vec<Genome> = args
        .genomes
        .iter()
        .map(|path| Genome::load(path))
        .collect::<Result<_>>()?;
    let config = SimConfig {
        cells: args.cells,
        read_length: args.read_length,
        min_coverage: args.min_coverage,
        max_coverage: args.max_coverage,
        error_rate: args.error_rate,
        seed: args.seed,
        gzip: args.gzip,
    };
    let plans = simulate::run(&genomes, &config, &args.out_dir)?;
    let total_reads: usize = plans.iter().map(|p| p.n_reads).sum();
    eprintln!(
        "simulated {} cells ({} reads) into {}",
        plans.len(),
        total_reads,
        args.out_dir.display()
    );
    Ok(())
}

fn run_cluster(args: ClusterArgs) -> Result<()> {
    let file = File::open(&args.matrix)
        .with_context(|| format!("cannot open {}", args.matrix.display()))?;
    let matrix = read_matrix(BufReader::new(file))?;
    let linkage = if args.distance {
        Linkage::Distance
    } else {
        Linkage::Similarity
    };
    let clusters = single_linkage(&matrix.values, args.threshold, linkage);
    let labels: Vec<String> = matrix
        .names
        .iter()
        .map(|name| label_from_name(name, &args.label_separator))
        .collect();

    if let Some(path) = &args.assignments {
        let mut out = open_output(Some(path))?;
        writeln!(out, "cell\tcluster\tlabel")?;
        for ((name, cluster), label) in matrix.names.iter().zip(&clusters).zip(&labels) {
            writeln!(out, "{name}\t{cluster}\t{label}")?;
        }
        out.flush()?;
    }

    let report = make_report(&labels, &clusters);
    let mut out = open_output(None)?;
    write_report(&mut out, &report)?;
    out.flush()?;
    Ok(())
}
