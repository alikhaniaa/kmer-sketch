//! A fixed, deterministic 64-bit hash for encoded k-mers.
//!
//! Sketches written by one run must be comparable with sketches written by
//! another run (or another machine), so the hash must not be randomly seeded.
//! We use the SplitMix64 finaliser: a few shifts, xors and multiplications
//! that mix all input bits into all output bits. It is a bijection on u64,
//! so two different k-mers never get the same hash.

/// Hash an encoded (canonical) k-mer to a well-mixed 64-bit value.
pub fn hash_kmer(kmer: u64) -> u64 {
    let mut z = kmer.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values_are_stable() {
        // The first two outputs of a SplitMix64 generator seeded with 0 are
        // the published reference values below. The Python scripts in
        // scripts/ implement the same function.
        assert_eq!(hash_kmer(0), 0xE220_A839_7B1D_CDAF);
        assert_eq!(hash_kmer(0x9E37_79B9_7F4A_7C15), 0x6E78_9E6A_A1B9_65F4);
    }

    #[test]
    fn different_inputs_give_different_hashes() {
        let hashes: std::collections::HashSet<u64> = (0..10_000u64).map(hash_kmer).collect();
        assert_eq!(hashes.len(), 10_000);
    }
}
