//! End-to-end tests that run the compiled `kmer-sketch` binary.

use assert_cmd::Command;
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::fs;
use std::io::Write;
use std::path::Path;

fn kmer_sketch() -> Command {
    Command::cargo_bin("kmer-sketch").expect("binary is built by cargo test")
}

fn random_dna(seed: u64, len: usize) -> String {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..len)
        .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
        .collect()
}

fn write_fasta(path: &Path, name: &str, seq: &str) {
    fs::write(path, format!(">{name}\n{seq}\n")).unwrap();
}

#[test]
fn count_small_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("tiny.fa");
    // Two records; the N splits the second record.
    fs::write(&fasta, ">a\nACGTA\n>b\nAANTT\n").unwrap();

    // k=2: record a gives AC x2, CG, TA; record b gives AA (from AA) and AA (from TT).
    kmer_sketch()
        .args(["count", "-k", "2"])
        .arg(&fasta)
        .assert()
        .success()
        .stdout("kmer\tcount\nAA\t2\nAC\t2\nCG\t1\nTA\t1\n");

    kmer_sketch()
        .args(["count", "-k", "2", "--min-count", "2"])
        .arg(&fasta)
        .assert()
        .success()
        .stdout("kmer\tcount\nAA\t2\nAC\t2\n");
}

#[test]
fn gzip_fastq_gives_same_counts_as_plain() {
    let dir = tempfile::tempdir().unwrap();
    let fastq_text = "@r1\nACGTACGTTGCA\n+\nIIIIIIIIIIII\n@r2\nTTGCANNACGTA\n+\nIIIIIIIIIIII\n";
    let plain = dir.path().join("reads.fastq");
    fs::write(&plain, fastq_text).unwrap();
    let gz = dir.path().join("reads.fastq.gz");
    let mut encoder = GzEncoder::new(fs::File::create(&gz).unwrap(), Compression::default());
    encoder.write_all(fastq_text.as_bytes()).unwrap();
    encoder.finish().unwrap();

    let plain_out = kmer_sketch()
        .args(["count", "-k", "5"])
        .arg(&plain)
        .output()
        .unwrap();
    let gz_out = kmer_sketch()
        .args(["count", "-k", "5"])
        .arg(&gz)
        .output()
        .unwrap();
    assert!(plain_out.status.success() && gz_out.status.success());
    assert_eq!(plain_out.stdout, gz_out.stdout);
    assert!(plain_out.stdout.len() > "kmer\tcount\n".len());
}

#[test]
fn rejects_k_above_31() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("x.fa");
    write_fasta(&fasta, "x", "ACGT");
    kmer_sketch()
        .args(["count", "-k", "32"])
        .arg(&fasta)
        .assert()
        .failure();
}

#[test]
fn full_pipeline_recovers_three_random_species() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let genomes: Vec<_> = ["alpha", "beta", "gamma"]
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let path = root.join(format!("{name}.fa"));
            write_fasta(&path, name, &random_dna(100 + i as u64, 60_000));
            path
        })
        .collect();

    let cells = root.join("cells");
    kmer_sketch()
        .args(["simulate", "--cells", "9", "--seed", "1", "--gzip", "-o"])
        .arg(&cells)
        .args(&genomes)
        .assert()
        .success();
    assert!(cells.join("alpha__cell000.fastq.gz").exists());
    assert!(cells.join("cells.tsv").exists());

    let mut cell_files: Vec<_> = fs::read_dir(&cells)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.to_string_lossy().ends_with(".fastq.gz"))
        .collect();
    cell_files.sort();
    assert_eq!(cell_files.len(), 9);

    let sketches = root.join("sketches");
    kmer_sketch()
        .args(["sketch", "--scaled", "100", "-o"])
        .arg(&sketches)
        .args(&cell_files)
        .assert()
        .success();

    let mut sketch_files: Vec<_> = fs::read_dir(&sketches)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect();
    sketch_files.sort();
    let matrix = root.join("jaccard.tsv");
    kmer_sketch()
        .args(["compare", "-o"])
        .arg(&matrix)
        .args(&sketch_files)
        .assert()
        .success();

    let output = kmer_sketch()
        .args(["cluster", "--threshold", "0.1"])
        .arg(&matrix)
        .output()
        .unwrap();
    assert!(output.status.success());
    let report = String::from_utf8(output.stdout).unwrap();
    assert!(
        report.starts_with(
            "# cells=9 clusters=3 labels=3 ARI=1.0000 purity=1.0000 species_recovered=3/3"
        ),
        "unexpected report:\n{report}"
    );
}

#[test]
fn compare_rejects_mixed_k() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("g.fa");
    write_fasta(&fasta, "g", &random_dna(5, 5_000));
    for k in ["15", "21"] {
        kmer_sketch()
            .args(["sketch", "-k", k, "--scaled", "10", "-o"])
            .arg(dir.path().join(format!("k{k}")))
            .arg(&fasta)
            .assert()
            .success();
    }
    kmer_sketch()
        .arg("compare")
        .arg(dir.path().join("k15/g.sketch"))
        .arg(dir.path().join("k21/g.sketch"))
        .assert()
        .failure();
}

#[test]
fn empty_inputs_do_not_abort_the_batch() {
    let dir = tempfile::tempdir().unwrap();
    let full = dir.path().join("full.fa");
    write_fasta(&full, "full", &random_dna(9, 20_000));
    let empty = dir.path().join("empty.fastq");
    fs::write(&empty, b"").unwrap();
    let empty_gz = dir.path().join("emptygz.fastq.gz");
    GzEncoder::new(fs::File::create(&empty_gz).unwrap(), Compression::default())
        .finish()
        .unwrap();

    // count: empty files add nothing.
    kmer_sketch()
        .args(["count", "-k", "3"])
        .arg(&empty)
        .arg(&empty_gz)
        .assert()
        .success()
        .stdout("kmer\tcount\n");

    // sketch: all three files are sketched; the empty ones get a warning.
    let sketches = dir.path().join("sketches");
    kmer_sketch()
        .args(["sketch", "--scaled", "10", "-o"])
        .arg(&sketches)
        .args([&full, &empty, &empty_gz])
        .assert()
        .success()
        .stderr(predicates::str::contains("sketch 'empty' has 0 hashes"));

    // compare: Jaccard 0 against the empty sketch, also on its own diagonal.
    let output = kmer_sketch()
        .arg("compare")
        .arg(sketches.join("full.sketch"))
        .arg(sketches.join("empty.sketch"))
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        "name\tfull\tempty\nfull\t1.000000\t0.000000\nempty\t0.000000\t0.000000\n"
    );
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("0 hashes"));
}

#[test]
fn simulate_rejects_bad_coverage_without_panicking() {
    let dir = tempfile::tempdir().unwrap();
    let genome = dir.path().join("g.fa");
    write_fasta(&genome, "g", &random_dna(3, 2_000));
    for (min, max) in [("1", "inf"), ("1e30", "1e30")] {
        // Exit code 1 is an ordinary error; a panic would exit with 101.
        kmer_sketch()
            .args([
                "simulate",
                "--cells",
                "1",
                "--min-coverage",
                min,
                "--max-coverage",
                max,
            ])
            .arg("-o")
            .arg(dir.path().join("cells"))
            .arg(&genome)
            .assert()
            .code(1);
    }
}

#[test]
fn directory_as_input_is_a_clear_error() {
    let dir = tempfile::tempdir().unwrap();
    kmer_sketch()
        .args(["count", "-k", "5"])
        .arg(dir.path())
        .assert()
        .code(1)
        .stderr(predicates::str::contains("is not a file"));
}
