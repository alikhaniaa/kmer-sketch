//! FracMinHash ("scaled MinHash") sketches.
//!
//! A sketch keeps the hash of every distinct canonical k-mer whose hash is
//! at most `u64::MAX / scaled`. With `scaled = 1000` about one k-mer in a
//! thousand is kept. Because the same rule is used for every sample, a k-mer
//! that is kept in one sketch is kept in every sketch that contains it, so
//! set operations on sketches estimate set operations on the full k-mer sets.
//!
//! # File format (text, one sketch per file)
//!
//! ```text
//! # kmer-sketch FracMinHash v1
//! name<TAB>ecoli__cell001
//! ksize<TAB>21
//! scaled<TAB>1000
//! hash<TAB>splitmix64
//! count<TAB>3
//! hashes
//! 1234567
//! 2345678
//! 3456789
//! ```
//!
//! Lines starting with `#` are comments. After the `hashes` line come
//! `count` hash values in decimal, sorted ascending, one per line.

use crate::hash::hash_kmer;
use crate::kmer::{canonical_kmers, check_k};
use crate::seqio::for_each_sequence;
use anyhow::{bail, Context, Result};
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Name of the hash function, written into sketch files.
const HASH_NAME: &str = "splitmix64";

/// A FracMinHash sketch of one sample.
#[derive(Debug, Clone, PartialEq)]
pub struct Sketch {
    pub name: String,
    pub ksize: usize,
    pub scaled: u64,
    /// Kept hash values, sorted ascending, no duplicates.
    pub hashes: Vec<u64>,
}

/// Largest hash value that is kept for a given `scaled`.
pub fn max_hash(scaled: u64) -> u64 {
    u64::MAX / scaled
}

fn check_parameters(ksize: usize, scaled: u64) -> Result<()> {
    check_k(ksize)?;
    if scaled == 0 {
        bail!("scaled must be at least 1");
    }
    Ok(())
}

/// Hash each k-mer and keep the hashes that are <= `threshold`.
fn keep_small_hashes(
    kept: &mut FxHashSet<u64>,
    threshold: u64,
    kmers: impl IntoIterator<Item = u64>,
) {
    for kmer in kmers {
        let hash = hash_kmer(kmer);
        if hash <= threshold {
            kept.insert(hash);
        }
    }
}

impl Sketch {
    /// Build a sketch from an iterator of canonical k-mers.
    pub fn from_kmers(
        name: &str,
        ksize: usize,
        scaled: u64,
        kmers: impl IntoIterator<Item = u64>,
    ) -> Result<Self> {
        check_parameters(ksize, scaled)?;
        let mut kept = FxHashSet::default();
        keep_small_hashes(&mut kept, max_hash(scaled), kmers);
        Ok(Self::from_hash_set(name, ksize, scaled, kept))
    }

    /// Sketch every sequence of a FASTA/FASTQ file together.
    pub fn from_file(path: &Path, name: &str, ksize: usize, scaled: u64) -> Result<Self> {
        check_parameters(ksize, scaled)?;
        let threshold = max_hash(scaled);
        let mut kept = FxHashSet::default();
        for_each_sequence(path, |seq| {
            keep_small_hashes(&mut kept, threshold, canonical_kmers(seq, ksize));
        })?;
        Ok(Self::from_hash_set(name, ksize, scaled, kept))
    }

    fn from_hash_set(name: &str, ksize: usize, scaled: u64, kept: FxHashSet<u64>) -> Self {
        let mut hashes: Vec<u64> = kept.into_iter().collect();
        hashes.sort_unstable();
        Sketch {
            name: name.to_string(),
            ksize,
            scaled,
            hashes,
        }
    }

    /// Write the sketch in the text format described at the top of this file.
    pub fn write_to(&self, out: &mut impl Write) -> Result<()> {
        writeln!(out, "# kmer-sketch FracMinHash v1")?;
        writeln!(out, "name\t{}", self.name)?;
        writeln!(out, "ksize\t{}", self.ksize)?;
        writeln!(out, "scaled\t{}", self.scaled)?;
        writeln!(out, "hash\t{HASH_NAME}")?;
        writeln!(out, "count\t{}", self.hashes.len())?;
        writeln!(out, "hashes")?;
        for hash in &self.hashes {
            writeln!(out, "{hash}")?;
        }
        Ok(())
    }

    /// Read a sketch written by [`Sketch::write_to`].
    pub fn read_from(input: impl BufRead) -> Result<Self> {
        let mut name = None;
        let mut ksize = None;
        let mut scaled = None;
        let mut count = None;
        let mut hashes = Vec::new();
        let mut in_hash_section = false;

        for (line_index, line) in input.lines().enumerate() {
            let line = line?;
            let line_number = line_index + 1;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            if in_hash_section {
                let hash: u64 = line
                    .parse()
                    .with_context(|| format!("line {line_number}: bad hash value '{line}'"))?;
                hashes.push(hash);
                continue;
            }
            if line == "hashes" {
                in_hash_section = true;
                continue;
            }
            let (key, value) = line
                .split_once('\t')
                .with_context(|| format!("line {line_number}: expected 'key<TAB>value'"))?;
            let bad_value = || format!("line {line_number}: bad value '{value}' for '{key}'");
            match key {
                "name" => name = Some(value.to_string()),
                "ksize" => ksize = Some(value.parse::<usize>().with_context(bad_value)?),
                "scaled" => scaled = Some(value.parse::<u64>().with_context(bad_value)?),
                "count" => count = Some(value.parse::<usize>().with_context(bad_value)?),
                "hash" => {
                    if value != HASH_NAME {
                        bail!("line {line_number}: unsupported hash function '{value}'");
                    }
                }
                _ => bail!("line {line_number}: unknown key '{key}'"),
            }
        }

        let sketch = Sketch {
            name: name.context("missing 'name'")?,
            ksize: ksize.context("missing 'ksize'")?,
            scaled: scaled.context("missing 'scaled'")?,
            hashes,
        };
        check_parameters(sketch.ksize, sketch.scaled)?;
        let count = count.context("missing 'count'")?;
        if count != sketch.hashes.len() {
            bail!(
                "header says {count} hashes but file has {}",
                sketch.hashes.len()
            );
        }
        if !sketch.hashes.windows(2).all(|pair| pair[0] < pair[1]) {
            bail!("hashes must be sorted ascending without duplicates");
        }
        // Sorted, so only the last hash can be the largest.
        if let Some(&largest) = sketch.hashes.last() {
            if largest > max_hash(sketch.scaled) {
                bail!(
                    "hash {largest} is above the limit for scaled={}",
                    sketch.scaled
                );
            }
        }
        Ok(sketch)
    }

    /// Save the sketch to a file.
    pub fn save(&self, path: &Path) -> Result<()> {
        let file =
            File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
        let mut out = BufWriter::new(file);
        self.write_to(&mut out)?;
        out.flush()?;
        Ok(())
    }

    /// Load a sketch from a file.
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
        Self::read_from(BufReader::new(file))
            .with_context(|| format!("cannot parse sketch {}", path.display()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare::jaccard;
    use proptest::prelude::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::collections::HashSet;

    fn random_dna(rng: &mut StdRng, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }

    fn exact_jaccard(a: &[u8], b: &[u8], k: usize) -> f64 {
        let set_a: HashSet<u64> = canonical_kmers(a, k).collect();
        let set_b: HashSet<u64> = canonical_kmers(b, k).collect();
        let shared = set_a.intersection(&set_b).count() as f64;
        let union = set_a.union(&set_b).count() as f64;
        shared / union
    }

    #[test]
    fn scaled_one_gives_exact_jaccard() {
        let mut rng = StdRng::seed_from_u64(1);
        let a = random_dna(&mut rng, 5_000);
        let mut b = a[..3_000].to_vec();
        b.extend(random_dna(&mut rng, 2_000));
        let k = 15;
        let sa = Sketch::from_kmers("a", k, 1, canonical_kmers(&a, k)).unwrap();
        let sb = Sketch::from_kmers("b", k, 1, canonical_kmers(&b, k)).unwrap();
        let estimate = jaccard(&sa.hashes, &sb.hashes);
        assert!((estimate - exact_jaccard(&a, &b, k)).abs() < 1e-12);
    }

    #[test]
    fn scaled_estimate_is_close_to_exact_jaccard() {
        // Two related 200 kb sequences: b shares the first 60% of a.
        let mut rng = StdRng::seed_from_u64(2);
        let a = random_dna(&mut rng, 200_000);
        let mut b = a[..120_000].to_vec();
        b.extend(random_dna(&mut rng, 80_000));
        let k = 21;
        let exact = exact_jaccard(&a, &b, k);
        let sa = Sketch::from_kmers("a", k, 50, canonical_kmers(&a, k)).unwrap();
        let sb = Sketch::from_kmers("b", k, 50, canonical_kmers(&b, k)).unwrap();
        let estimate = jaccard(&sa.hashes, &sb.hashes);
        // About 8000 hashes are sampled from the union, so the standard error
        // is around 0.005. A tolerance of 0.03 is very safe.
        assert!(
            (estimate - exact).abs() < 0.03,
            "estimate {estimate} vs exact {exact}"
        );
    }

    #[test]
    fn sketch_keeps_about_one_in_scaled() {
        let mut rng = StdRng::seed_from_u64(3);
        let a = random_dna(&mut rng, 500_000);
        let sketch = Sketch::from_kmers("a", 21, 100, canonical_kmers(&a, 21)).unwrap();
        // ~500 000 distinct k-mers / 100 = ~5000 expected; allow +-10%.
        let n = sketch.hashes.len();
        assert!((4_500..=5_500).contains(&n), "kept {n}");
    }

    #[test]
    fn rejects_count_mismatch() {
        let text = "name\tx\nksize\t21\nscaled\t10\nhash\tsplitmix64\ncount\t2\nhashes\n5\n";
        assert!(Sketch::read_from(text.as_bytes()).is_err());
    }

    #[test]
    fn rejects_other_hash_function() {
        let text = "name\tx\nksize\t21\nscaled\t10\nhash\tmurmur3\ncount\t0\nhashes\n";
        assert!(Sketch::read_from(text.as_bytes()).is_err());
    }

    /// Header for a hand-written sketch file.
    fn header(ksize: &str, scaled: &str, count: usize) -> String {
        format!(
            "name\tx\nksize\t{ksize}\nscaled\t{scaled}\nhash\tsplitmix64\ncount\t{count}\nhashes\n"
        )
    }

    #[test]
    fn rejects_bad_ksize_and_scaled() {
        for (ksize, scaled) in [("0", "10"), ("32", "10"), ("21", "0")] {
            let text = header(ksize, scaled, 0);
            assert!(
                Sketch::read_from(text.as_bytes()).is_err(),
                "{ksize} {scaled}"
            );
        }
    }

    #[test]
    fn rejects_hash_above_threshold() {
        // With scaled=10 no kept hash can be larger than u64::MAX / 10.
        let text = header("21", "10", 1) + &format!("{}\n", u64::MAX);
        assert!(Sketch::read_from(text.as_bytes()).is_err());
    }

    #[test]
    fn parse_errors_name_the_line() {
        let text = header("twenty", "10", 0);
        let message = format!("{:#}", Sketch::read_from(text.as_bytes()).unwrap_err());
        assert!(message.contains("line 2"), "{message}");
    }

    #[test]
    fn empty_sketch_round_trips() {
        let text = header("21", "1000", 0);
        let sketch = Sketch::read_from(text.as_bytes()).unwrap();
        assert!(sketch.hashes.is_empty());
    }

    /// Strategy: a scaled value and a set of hashes that respect its threshold.
    fn scaled_and_hashes() -> impl Strategy<Value = (u64, std::collections::BTreeSet<u64>)> {
        (1u64..100_000).prop_flat_map(|scaled| {
            (
                Just(scaled),
                proptest::collection::btree_set(0..=max_hash(scaled), 0..200),
            )
        })
    }

    proptest! {
        #[test]
        fn write_then_read_round_trips(
            name in "[A-Za-z0-9_.]{1,20}",
            ksize in 1usize..=31,
            (scaled, raw) in scaled_and_hashes(),
        ) {
            let sketch = Sketch { name, ksize, scaled, hashes: raw.into_iter().collect() };
            let mut buffer = Vec::new();
            sketch.write_to(&mut buffer).unwrap();
            let parsed = Sketch::read_from(buffer.as_slice()).unwrap();
            prop_assert_eq!(parsed, sketch);
        }
    }
}
