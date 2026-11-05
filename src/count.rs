//! Canonical k-mer counting.
//!
//! Counts are kept in `NUM_SHARDS` separate hash maps ("shards"). Every k-mer
//! belongs to exactly one shard, chosen from its hash. This makes the work
//! easy to parallelise without locks:
//!
//! 1. In parallel over chunks of sequences, extract canonical k-mers and drop
//!    each one into the bucket of its shard.
//! 2. In parallel over shards, each shard's map takes the k-mers from its own
//!    buckets. No two threads ever touch the same map.

use crate::hash::hash_kmer;
use crate::kmer::{canonical_kmers, check_k, decode};
use crate::seqio::for_each_batch;
use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::io::Write;
use std::path::Path;

/// Number of shards. Must be a power of two (we take the top bits of the hash).
const NUM_SHARDS: usize = 64;
/// log2(NUM_SHARDS): how many top hash bits select the shard.
const SHARD_BITS: u32 = 6;
/// Roughly how many bases to read from a file before counting a batch.
const BATCH_BASES: usize = 16 * 1024 * 1024;

/// Which shard a k-mer belongs to.
fn shard_of(kmer: u64) -> usize {
    (hash_kmer(kmer) >> (64 - SHARD_BITS)) as usize
}

/// Canonical k-mer counts, split into shards.
pub struct KmerCounter {
    k: usize,
    shards: Vec<FxHashMap<u64, u64>>,
}

impl KmerCounter {
    /// Create an empty counter for k-mers of length `k` (1..=31).
    pub fn new(k: usize) -> Result<Self> {
        check_k(k)?;
        let shards = (0..NUM_SHARDS).map(|_| FxHashMap::default()).collect();
        Ok(KmerCounter { k, shards })
    }

    /// Add all canonical k-mers of one batch of sequences.
    pub fn add_sequences(&mut self, sequences: &[Vec<u8>]) {
        let k = self.k;

        // Phase 1: split the batch into a few chunks per thread, and for each
        // chunk collect its k-mers into one bucket per shard.
        let chunk_size = (sequences.len() / (4 * rayon::current_num_threads())).max(1);
        let buckets_per_chunk: Vec<Vec<Vec<u64>>> = sequences
            .par_chunks(chunk_size)
            .map(|chunk| bucket_kmers(chunk, k))
            .collect();

        // Phase 2: every shard map is updated by exactly one thread.
        self.shards
            .par_iter_mut()
            .enumerate()
            .for_each(|(shard_index, map)| {
                for buckets in &buckets_per_chunk {
                    for &kmer in &buckets[shard_index] {
                        *map.entry(kmer).or_insert(0) += 1;
                    }
                }
            });
    }

    /// Number of distinct k-mers counted so far.
    pub fn distinct_kmers(&self) -> usize {
        self.shards.iter().map(|map| map.len()).sum()
    }

    /// All (k-mer, count) pairs with count >= `min_count`, sorted by k-mer.
    /// Because of the 2-bit layout, this is alphabetical order of the k-mers.
    pub fn into_sorted_counts(self, min_count: u64) -> Vec<(u64, u64)> {
        let mut counts: Vec<(u64, u64)> = self
            .shards
            .into_iter()
            .flat_map(|map| map.into_iter())
            .filter(|&(_, count)| count >= min_count)
            .collect();
        counts.par_sort_unstable_by_key(|&(kmer, _)| kmer);
        counts
    }
}

/// Put the canonical k-mers of some sequences into one bucket per shard.
fn bucket_kmers(sequences: &[Vec<u8>], k: usize) -> Vec<Vec<u64>> {
    let mut buckets = vec![Vec::new(); NUM_SHARDS];
    for seq in sequences {
        for kmer in canonical_kmers(seq, k) {
            buckets[shard_of(kmer)].push(kmer);
        }
    }
    buckets
}

/// Count canonical k-mers over all input files (counts are summed).
pub fn count_files(paths: &[impl AsRef<Path>], k: usize) -> Result<KmerCounter> {
    let mut counter = KmerCounter::new(k)?;
    for path in paths {
        for_each_batch(path.as_ref(), BATCH_BASES, |batch| {
            counter.add_sequences(batch)
        })?;
    }
    Ok(counter)
}

/// Write counts as TSV with a header line: `kmer<TAB>count`.
pub fn write_counts(out: &mut impl Write, counts: &[(u64, u64)], k: usize) -> Result<()> {
    writeln!(out, "kmer\tcount")?;
    for &(kmer, count) in counts {
        writeln!(out, "{}\t{}", decode(kmer, k), count)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::encode;
    use proptest::prelude::*;
    use std::collections::HashMap;

    fn as_strings(counts: &[(u64, u64)], k: usize) -> Vec<(String, u64)> {
        counts.iter().map(|&(km, c)| (decode(km, k), c)).collect()
    }

    #[test]
    fn counts_small_example_by_hand() {
        // k=2 windows of ACGTA: AC, CG, GT, TA
        // canonical: AC, CG (palindrome), AC (GT -> AC), TA (palindrome)
        let mut counter = KmerCounter::new(2).unwrap();
        counter.add_sequences(&[b"ACGTA".to_vec()]);
        let counts = counter.into_sorted_counts(1);
        let expected = vec![
            ("AC".to_string(), 2),
            ("CG".to_string(), 1),
            ("TA".to_string(), 1),
        ];
        assert_eq!(as_strings(&counts, 2), expected);
    }

    #[test]
    fn min_count_filters() {
        let mut counter = KmerCounter::new(2).unwrap();
        counter.add_sequences(&[b"ACGTA".to_vec()]);
        let counts = counter.into_sorted_counts(2);
        assert_eq!(counts, vec![(encode(b"AC").unwrap(), 2)]);
    }

    #[test]
    fn writes_tsv() {
        let counts = vec![(encode(b"AAC").unwrap(), 3), (encode(b"ACG").unwrap(), 1)];
        let mut out = Vec::new();
        write_counts(&mut out, &counts, 3).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "kmer\tcount\nAAC\t3\nACG\t1\n"
        );
    }

    proptest! {
        /// Counting in several batches must equal a plain HashMap count.
        #[test]
        fn batched_count_matches_simple_count(
            seqs in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"ACGTN".to_vec()), 0..80),
                0..20),
            k in 1usize..=12,
            split in 0usize..20,
        ) {
            let mut expected: HashMap<u64, u64> = HashMap::new();
            for seq in &seqs {
                for kmer in canonical_kmers(seq, k) {
                    *expected.entry(kmer).or_insert(0) += 1;
                }
            }
            let mut expected: Vec<(u64, u64)> = expected.into_iter().collect();
            expected.sort();

            let split = split.min(seqs.len());
            let mut counter = KmerCounter::new(k).unwrap();
            counter.add_sequences(&seqs[..split]);
            counter.add_sequences(&seqs[split..]);
            prop_assert_eq!(counter.into_sorted_counts(1), expected);
        }
    }
}
