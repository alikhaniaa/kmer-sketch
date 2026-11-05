//! 2-bit k-mer encoding, reverse complement, canonical k-mers and a rolling
//! iterator that produces every canonical k-mer of a sequence.
//!
//! Encoding: A=0, C=1, G=2, T=3. A k-mer of length k is stored in the low
//! 2*k bits of a `u64`, with the FIRST base in the MOST significant position.
//! With this layout, comparing two encoded k-mers as integers gives the same
//! order as comparing the strings alphabetically (A < C < G < T).
//!
//! Complement in this code is simply `3 - code` (A<->T, C<->G).

use anyhow::{bail, Result};

/// Largest supported k. With k <= 31 a k-mer uses at most 62 bits, so the
/// mask `(1 << 2k) - 1` never overflows.
pub const MAX_K: usize = 31;

/// Check that `k` is in the supported range 1..=31.
pub fn check_k(k: usize) -> Result<()> {
    if k == 0 || k > MAX_K {
        bail!("k must be between 1 and {MAX_K}, got {k}");
    }
    Ok(())
}

/// Map a nucleotide byte to its 2-bit code. Lowercase is accepted.
/// Returns `None` for anything that is not A/C/G/T (N, IUPAC codes, gaps).
pub fn encode_base(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Map a 2-bit code back to an uppercase base.
pub fn decode_base(code: u64) -> u8 {
    b"ACGT"[(code & 3) as usize]
}

/// Bit mask covering the low 2*k bits.
pub fn kmer_mask(k: usize) -> u64 {
    (1u64 << (2 * k)) - 1
}

/// Encode a whole k-mer. Returns `None` if it contains a non-ACGT base.
pub fn encode(kmer: &[u8]) -> Option<u64> {
    let mut code = 0u64;
    for &base in kmer {
        code = (code << 2) | encode_base(base)?;
    }
    Some(code)
}

/// Decode an encoded k-mer of length `k` back into a string.
pub fn decode(code: u64, k: usize) -> String {
    let mut bases = Vec::with_capacity(k);
    for i in (0..k).rev() {
        bases.push(decode_base(code >> (2 * i)));
    }
    // Only A/C/G/T bytes were pushed, so this is valid UTF-8.
    String::from_utf8(bases).expect("decoded bases are ASCII")
}

/// Reverse complement of an encoded k-mer of length `k`.
///
/// We read bases off the low end of `code` (the LAST base first) and push
/// their complements onto the result, so the result starts with the
/// complement of the last base.
pub fn reverse_complement(code: u64, k: usize) -> u64 {
    let mut remaining = code;
    let mut rc = 0u64;
    for _ in 0..k {
        let base = remaining & 3;
        rc = (rc << 2) | (3 - base);
        remaining >>= 2;
    }
    rc
}

/// Canonical form of a k-mer: the smaller of the k-mer and its reverse
/// complement. A k-mer and its reverse complement have the same canonical
/// form, so both DNA strands are counted as the same k-mer.
pub fn canonical(code: u64, k: usize) -> u64 {
    code.min(reverse_complement(code, k))
}

/// Iterator over the canonical k-mers of a sequence, computed with a rolling
/// update: each new base costs O(1) instead of re-encoding k bases.
///
/// Any base that is not A/C/G/T resets the window, so k-mers containing N
/// are skipped.
pub struct CanonicalKmers<'a> {
    seq: &'a [u8],
    k: usize,
    /// Index of the next base to read from `seq`.
    pos: usize,
    /// Number of consecutive A/C/G/T bases read since the last reset.
    valid_len: usize,
    /// Forward k-mer: new bases enter at the low end.
    forward: u64,
    /// Reverse complement k-mer: new complemented bases enter at the high end.
    reverse: u64,
    mask: u64,
    /// Shift that places a 2-bit code at the most significant base position.
    top_shift: usize,
}

impl<'a> CanonicalKmers<'a> {
    /// Create the iterator. `k` must be in 1..=31 (see [`check_k`]).
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!((1..=MAX_K).contains(&k), "k must be in 1..=31");
        CanonicalKmers {
            seq,
            k,
            pos: 0,
            valid_len: 0,
            forward: 0,
            reverse: 0,
            mask: kmer_mask(k),
            top_shift: 2 * (k - 1),
        }
    }
}

impl Iterator for CanonicalKmers<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        while self.pos < self.seq.len() {
            let base = self.seq[self.pos];
            self.pos += 1;

            match encode_base(base) {
                Some(code) => {
                    // Forward strand: shift left, add the base at the end,
                    // drop the base that fell out of the window.
                    self.forward = ((self.forward << 2) | code) & self.mask;
                    // Reverse strand: shift right (drops the oldest base)
                    // and put the complement of the new base at the front.
                    self.reverse = (self.reverse >> 2) | ((3 - code) << self.top_shift);
                    self.valid_len += 1;

                    if self.valid_len >= self.k {
                        return Some(self.forward.min(self.reverse));
                    }
                }
                None => {
                    // Not ACGT: start a fresh window after this base.
                    self.valid_len = 0;
                    self.forward = 0;
                    self.reverse = 0;
                }
            }
        }
        None
    }
}

/// Convenience function: iterate over the canonical k-mers of `seq`.
pub fn canonical_kmers(seq: &[u8], k: usize) -> CanonicalKmers<'_> {
    CanonicalKmers::new(seq, k)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Slow reference: encode every window from scratch.
    fn naive_canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
        if seq.len() < k {
            return Vec::new();
        }
        seq.windows(k)
            .filter_map(encode)
            .map(|code| canonical(code, k))
            .collect()
    }

    #[test]
    fn encode_and_decode_known_kmer() {
        // A=00 C=01 G=10 T=11  ->  ACGT = 00 01 10 11 = 27
        assert_eq!(encode(b"ACGT"), Some(27));
        assert_eq!(decode(27, 4), "ACGT");
        assert_eq!(encode(b"acgt"), Some(27));
        assert_eq!(encode(b"ACNT"), None);
    }

    #[test]
    fn reverse_complement_known_kmer() {
        let k = 6;
        let code = encode(b"AACGTG").unwrap();
        assert_eq!(decode(reverse_complement(code, k), k), "CACGTT");
    }

    #[test]
    fn canonical_picks_alphabetically_smaller() {
        let k = 3;
        // TTT and AAA are reverse complements; AAA is smaller.
        let ttt = encode(b"TTT").unwrap();
        assert_eq!(decode(canonical(ttt, k), k), "AAA");
        // ACG has reverse complement CGT; ACG is smaller.
        let cgt = encode(b"CGT").unwrap();
        assert_eq!(decode(canonical(cgt, k), k), "ACG");
    }

    #[test]
    fn integer_order_matches_string_order() {
        let mut words = vec!["TGCA", "AAAA", "ACGT", "CAAA", "GTTT", "TTTT"];
        let mut codes: Vec<u64> = words
            .iter()
            .map(|w| encode(w.as_bytes()).unwrap())
            .collect();
        words.sort();
        codes.sort();
        let decoded: Vec<String> = codes.iter().map(|&c| decode(c, 4)).collect();
        assert_eq!(decoded, words);
    }

    #[test]
    fn rolling_skips_kmers_with_n() {
        let kmers: Vec<String> = canonical_kmers(b"ACGNACGTA", 3)
            .map(|c| decode(c, 3))
            .collect();
        // Valid windows: ACG, then (reset at N) ACG, CGT, GTA.
        // Canonical forms: ACG; ACG; CGT -> rc ACG is smaller; GTA (rc TAC is larger).
        assert_eq!(kmers, vec!["ACG", "ACG", "ACG", "GTA"]);
    }

    #[test]
    fn k_of_31_works() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTT";
        let rolling: Vec<u64> = canonical_kmers(seq, 31).collect();
        assert_eq!(rolling, naive_canonical_kmers(seq, 31));
        assert_eq!(rolling.len(), seq.len() - 31 + 1);
    }

    #[test]
    fn check_k_rejects_bad_values() {
        assert!(check_k(0).is_err());
        assert!(check_k(32).is_err());
        assert!(check_k(31).is_ok());
    }

    /// Strategy: a k in 1..=31 and a random k-mer value below 4^k.
    fn k_and_kmer() -> impl Strategy<Value = (usize, u64)> {
        (1usize..=MAX_K).prop_flat_map(|k| (Just(k), 0..=kmer_mask(k)))
    }

    proptest! {
        #[test]
        fn rc_of_rc_is_identity((k, code) in k_and_kmer()) {
            prop_assert_eq!(reverse_complement(reverse_complement(code, k), k), code);
        }

        #[test]
        fn decode_then_encode_is_identity((k, code) in k_and_kmer()) {
            prop_assert_eq!(encode(decode(code, k).as_bytes()), Some(code));
        }

        #[test]
        fn canonical_is_strand_independent((k, code) in k_and_kmer()) {
            let rc = reverse_complement(code, k);
            prop_assert_eq!(canonical(code, k), canonical(rc, k));
            prop_assert!(canonical(code, k) <= code);
        }

        #[test]
        fn rolling_matches_naive(
            seq in proptest::collection::vec(proptest::sample::select(b"ACGTNacgt".to_vec()), 0..300),
            k in 1usize..=MAX_K,
        ) {
            let rolling: Vec<u64> = canonical_kmers(&seq, k).collect();
            prop_assert_eq!(rolling, naive_canonical_kmers(&seq, k));
        }
    }
}
