//! Simulate single-cell read sets from reference genomes.
//!
//! Each simulated "cell" is drawn from one genome. It gets its own coverage,
//! drawn uniformly between `min_coverage` and `max_coverage`, so cells are
//! unevenly covered. Reads are random substrings of the genome taken from a
//! random strand, with random substitution errors. One FASTQ file is written
//! per cell, named `<species>__cell<NNN>.fastq`, so the true species can be
//! recovered from the file name later.
//!
//! Every cell uses its own random generator seeded from (seed, cell index),
//! so the output does not depend on how rayon schedules the cells.

use crate::seqio::{read_all_sequences, sample_name};
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Separator between species and cell number in file names.
pub const LABEL_SEPARATOR: &str = "__";

/// Refuse to simulate more reads than this for one cell (all reads of a cell
/// are held in memory before writing).
pub const MAX_READS_PER_CELL: usize = 1_000_000_000;

/// Settings for a simulation run.
#[derive(Debug, Clone)]
pub struct SimConfig {
    pub cells: usize,
    pub read_length: usize,
    pub min_coverage: f64,
    pub max_coverage: f64,
    pub error_rate: f64,
    pub seed: u64,
    pub gzip: bool,
}

/// A reference genome: a species name and one or more sequences
/// (chromosome, plasmids).
#[derive(Debug, Clone)]
pub struct Genome {
    pub species: String,
    pub records: Vec<Vec<u8>>,
}

impl Genome {
    /// Load a genome; the species name is the file name without extensions.
    pub fn load(path: &Path) -> Result<Self> {
        let species = sample_name(path);
        if species.contains(LABEL_SEPARATOR) {
            bail!("genome name '{species}' must not contain '{LABEL_SEPARATOR}'");
        }
        let records = read_all_sequences(path)?;
        Ok(Genome { species, records })
    }

    pub fn total_length(&self) -> usize {
        self.records.iter().map(Vec::len).sum()
    }
}

/// What will be simulated for one cell.
#[derive(Debug, Clone, PartialEq)]
pub struct CellPlan {
    pub index: usize,
    pub genome_index: usize,
    pub coverage: f64,
    pub n_reads: usize,
    pub file_name: String,
}

/// Decide genome, coverage and read count for every cell. Cells are given
/// to genomes in turn (cell 0 -> genome 0, cell 1 -> genome 1, ...).
pub fn plan_cells(genomes: &[Genome], config: &SimConfig) -> Result<Vec<CellPlan>> {
    validate(genomes, config)?;
    let mut rng = StdRng::seed_from_u64(config.seed);
    let mut plans = Vec::with_capacity(config.cells);
    for index in 0..config.cells {
        let genome_index = index % genomes.len();
        let genome = &genomes[genome_index];
        let coverage = rng.random_range(config.min_coverage..=config.max_coverage);
        let n_reads = reads_for_coverage(coverage, genome.total_length(), config.read_length);
        if n_reads > MAX_READS_PER_CELL {
            bail!(
                "cell {index} would need {n_reads} reads (limit {MAX_READS_PER_CELL}); lower the coverage"
            );
        }
        let extension = if config.gzip { "fastq.gz" } else { "fastq" };
        let file_name = format!(
            "{}{LABEL_SEPARATOR}cell{index:03}.{extension}",
            genome.species
        );
        plans.push(CellPlan {
            index,
            genome_index,
            coverage,
            n_reads,
            file_name,
        });
    }
    Ok(plans)
}

fn validate(genomes: &[Genome], config: &SimConfig) -> Result<()> {
    if genomes.is_empty() {
        bail!("need at least one genome");
    }
    if config.read_length == 0 {
        bail!("read length must be at least 1");
    }
    // Written so that NaN also fails (every comparison with NaN is false).
    let coverage_ok = config.min_coverage > 0.0
        && config.min_coverage <= config.max_coverage
        && config.max_coverage.is_finite();
    if !coverage_ok {
        bail!("need 0 < min coverage <= max coverage (finite)");
    }
    if !(0.0..=1.0).contains(&config.error_rate) {
        bail!("error rate must be between 0 and 1");
    }
    for genome in genomes {
        if !genome.records.iter().any(|r| r.len() >= config.read_length) {
            bail!(
                "genome '{}' has no sequence as long as the read length",
                genome.species
            );
        }
    }
    Ok(())
}

/// Number of reads needed for a coverage: coverage * genome length / read length.
pub fn reads_for_coverage(coverage: f64, genome_length: usize, read_length: usize) -> usize {
    let reads = (coverage * genome_length as f64 / read_length as f64).round() as usize;
    reads.max(1)
}

/// Random generator for one cell. Mixing the cell index into the seed gives
/// each cell an independent, reproducible stream.
fn cell_rng(seed: u64, cell_index: usize) -> StdRng {
    let mixed = seed ^ (cell_index as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    StdRng::seed_from_u64(mixed)
}

/// Simulate all reads of one cell.
pub fn simulate_cell(genome: &Genome, plan: &CellPlan, config: &SimConfig) -> Vec<Vec<u8>> {
    let mut rng = cell_rng(config.seed, plan.index);
    let source = ReadSource::new(genome, config.read_length);
    (0..plan.n_reads)
        .map(|_| simulate_read(&mut rng, &source, config.error_rate))
        .collect()
}

/// The records of a genome that are long enough to hold a read. Built once
/// per cell, so picking a record for each read needs no allocation.
pub struct ReadSource<'a> {
    records: Vec<&'a [u8]>,
    total_length: usize,
    read_length: usize,
}

impl<'a> ReadSource<'a> {
    /// `validate` (called by `plan_cells`) guarantees at least one usable record.
    pub fn new(genome: &'a Genome, read_length: usize) -> Self {
        let records: Vec<&[u8]> = genome
            .records
            .iter()
            .filter(|r| r.len() >= read_length)
            .map(Vec::as_slice)
            .collect();
        let total_length = records.iter().map(|r| r.len()).sum();
        ReadSource {
            records,
            total_length,
            read_length,
        }
    }

    /// Pick a record with probability proportional to its length.
    fn pick_record(&self, rng: &mut StdRng) -> &'a [u8] {
        let mut target = rng.random_range(0..self.total_length);
        for &record in &self.records {
            if target < record.len() {
                return record;
            }
            target -= record.len();
        }
        // Unreachable because target < total_length, but return something valid.
        self.records[self.records.len() - 1]
    }
}

/// One read: a random substring from a random strand, with substitutions.
pub fn simulate_read(rng: &mut StdRng, source: &ReadSource, error_rate: f64) -> Vec<u8> {
    let record = source.pick_record(rng);
    let read_length = source.read_length;
    let start = rng.random_range(0..=record.len() - read_length);
    let mut read = record[start..start + read_length].to_vec();
    if rng.random_bool(0.5) {
        read = reverse_complement_bases(&read);
    }
    add_substitutions(rng, &mut read, error_rate);
    read
}

/// Reverse complement of a DNA string; bases other than ACGT become N.
pub fn reverse_complement_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Replace each A/C/G/T base, with probability `error_rate`, by one of the
/// three other bases (chosen uniformly).
fn add_substitutions(rng: &mut StdRng, read: &mut [u8], error_rate: f64) {
    const BASES: &[u8; 4] = b"ACGT";
    for base in read.iter_mut() {
        let Some(position) = BASES.iter().position(|b| b == &base.to_ascii_uppercase()) else {
            continue; // leave N and other symbols alone
        };
        if rng.random_bool(error_rate) {
            let shift = rng.random_range(1..4);
            *base = BASES[(position + shift) % 4];
        }
    }
}

/// Write reads as FASTQ, with a constant quality string.
pub fn write_fastq(out: &mut impl Write, cell_name: &str, reads: &[Vec<u8>]) -> Result<()> {
    for (number, read) in reads.iter().enumerate() {
        writeln!(out, "@{cell_name}:{number}")?;
        out.write_all(read)?;
        out.write_all(b"\n+\n")?;
        out.write_all(&vec![b'I'; read.len()])?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// Write one cell's FASTQ file, gzip-compressed if requested.
fn write_cell_file(path: &Path, cell_name: &str, reads: &[Vec<u8>], gzip: bool) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    if gzip {
        let mut encoder = GzEncoder::new(BufWriter::new(file), Compression::default());
        write_fastq(&mut encoder, cell_name, reads)?;
        encoder.finish()?.flush()?;
    } else {
        let mut out = BufWriter::new(file);
        write_fastq(&mut out, cell_name, reads)?;
        out.flush()?;
    }
    Ok(())
}

/// Simulate every cell (in parallel) and write one FASTQ per cell into
/// `out_dir`, plus a `cells.tsv` manifest. Returns the plans.
pub fn run(genomes: &[Genome], config: &SimConfig, out_dir: &Path) -> Result<Vec<CellPlan>> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("cannot create {}", out_dir.display()))?;
    let plans = plan_cells(genomes, config)?;

    plans.par_iter().try_for_each(|plan| -> Result<()> {
        let genome = &genomes[plan.genome_index];
        let reads = simulate_cell(genome, plan, config);
        let path: PathBuf = out_dir.join(&plan.file_name);
        write_cell_file(&path, &sample_name(&path), &reads, config.gzip)
    })?;

    write_manifest(&out_dir.join("cells.tsv"), genomes, &plans)?;
    Ok(plans)
}

/// Write `cells.tsv`: file, species, coverage, number of reads.
fn write_manifest(path: &Path, genomes: &[Genome], plans: &[CellPlan]) -> Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "file\tspecies\tcoverage\treads")?;
    for plan in plans {
        writeln!(
            out,
            "{}\t{}\t{:.3}\t{}",
            plan.file_name, genomes[plan.genome_index].species, plan.coverage, plan.n_reads
        )?;
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_genome(species: &str, len: usize, seed: u64) -> Genome {
        let mut rng = StdRng::seed_from_u64(seed);
        let seq = (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
        Genome {
            species: species.to_string(),
            records: vec![seq],
        }
    }

    fn config(error_rate: f64) -> SimConfig {
        SimConfig {
            cells: 6,
            read_length: 50,
            min_coverage: 2.0,
            max_coverage: 30.0,
            error_rate,
            seed: 7,
            gzip: false,
        }
    }

    #[test]
    fn plans_alternate_genomes_and_respect_coverage_range() {
        let genomes = vec![random_genome("a", 2_000, 1), random_genome("b", 3_000, 2)];
        let plans = plan_cells(&genomes, &config(0.0)).unwrap();
        assert_eq!(plans.len(), 6);
        for plan in &plans {
            assert_eq!(plan.genome_index, plan.index % 2);
            assert!((2.0..=30.0).contains(&plan.coverage));
            let genome_length = genomes[plan.genome_index].total_length();
            assert_eq!(
                plan.n_reads,
                reads_for_coverage(plan.coverage, genome_length, 50)
            );
        }
        assert_eq!(plans[0].file_name, "a__cell000.fastq");
        assert_eq!(plans[3].file_name, "b__cell003.fastq");
    }

    #[test]
    fn same_seed_gives_same_reads() {
        let genomes = vec![random_genome("a", 2_000, 1)];
        let cfg = config(0.01);
        let plans = plan_cells(&genomes, &cfg).unwrap();
        let first = simulate_cell(&genomes[0], &plans[0], &cfg);
        let second = simulate_cell(&genomes[0], &plans[0], &cfg);
        assert_eq!(first, second);
        // A different cell gets different reads.
        let other = simulate_cell(&genomes[0], &plans[1], &cfg);
        assert_ne!(first[..10], other[..10]);
    }

    #[test]
    fn error_free_reads_come_from_the_genome() {
        let genome = random_genome("a", 5_000, 3);
        let forward = String::from_utf8(genome.records[0].clone()).unwrap();
        let reverse = String::from_utf8(reverse_complement_bases(&genome.records[0])).unwrap();
        let mut rng = StdRng::seed_from_u64(11);
        let source = ReadSource::new(&genome, 100);
        for _ in 0..200 {
            let read = String::from_utf8(simulate_read(&mut rng, &source, 0.0)).unwrap();
            assert!(forward.contains(&read) || reverse.contains(&read));
        }
    }

    #[test]
    fn error_rate_is_roughly_respected() {
        let genome = random_genome("a", 10_000, 4);
        let mut rng = StdRng::seed_from_u64(12);
        let original = genome.records[0][..1000].to_vec();
        let mut read = original.clone();
        let mut changed = 0;
        for _ in 0..100 {
            read.copy_from_slice(&original);
            add_substitutions(&mut rng, &mut read, 0.05);
            changed += read.iter().zip(&original).filter(|(a, b)| a != b).count();
        }
        // 100 000 bases at 5% -> ~5000 errors.
        assert!((4_500..=5_500).contains(&changed), "changed {changed}");
    }

    #[test]
    fn reverse_complement_of_bases() {
        assert_eq!(reverse_complement_bases(b"AACGTN"), b"NACGTT".to_vec());
    }

    #[test]
    fn rejects_bad_settings() {
        let genomes = vec![random_genome("a", 100, 1)];
        let mut cfg = config(0.0);
        cfg.read_length = 500; // longer than the genome
        assert!(plan_cells(&genomes, &cfg).is_err());
        let mut cfg = config(0.0);
        cfg.min_coverage = 5.0;
        cfg.max_coverage = 1.0;
        assert!(plan_cells(&genomes, &cfg).is_err());
    }

    #[test]
    fn rejects_non_finite_coverage() {
        let genomes = vec![random_genome("a", 1_000, 1)];
        for (min, max) in [(1.0, f64::INFINITY), (f64::NAN, 2.0), (1.0, f64::NAN)] {
            let mut cfg = config(0.0);
            cfg.min_coverage = min;
            cfg.max_coverage = max;
            let err = plan_cells(&genomes, &cfg).unwrap_err();
            assert!(err.to_string().contains("finite"), "{err}");
        }
    }

    #[test]
    fn rejects_too_many_reads() {
        // 1e30x coverage is finite but would need far more than the read limit.
        let genomes = vec![random_genome("a", 1_000, 1)];
        let mut cfg = config(0.0);
        cfg.min_coverage = 1e30;
        cfg.max_coverage = 1e30;
        let err = plan_cells(&genomes, &cfg).unwrap_err();
        assert!(err.to_string().contains("limit"), "{err}");
    }

    #[test]
    fn picks_records_in_proportion_to_length() {
        // Records of length 300 and 100 (and one too short to use):
        // about 3/4 of the picks should come from the long one.
        let genome = Genome {
            species: "a".to_string(),
            records: vec![vec![b'A'; 300], vec![b'C'; 100], vec![b'G'; 10]],
        };
        let source = ReadSource::new(&genome, 50);
        let mut rng = StdRng::seed_from_u64(5);
        let long_picks = (0..10_000)
            .filter(|_| source.pick_record(&mut rng)[0] == b'A')
            .count();
        assert!((7_200..=7_800).contains(&long_picks), "{long_picks}");
    }
}
