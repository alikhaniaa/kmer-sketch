//! Pairwise comparison of sketches: Jaccard, containment and Mash distance.
//!
//! All functions work on sorted, duplicate-free hash lists (as stored in a
//! [`Sketch`]), so the intersection can be computed with a single merge pass.

use crate::sketch::Sketch;
use anyhow::{bail, Result};
use clap::ValueEnum;
use rayon::prelude::*;
use std::io::Write;

/// Which number to put in the comparison matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Metric {
    /// |A ∩ B| / |A ∪ B| (symmetric)
    Jaccard,
    /// |A ∩ B| / |A|: fraction of the ROW sample found in the COLUMN sample
    Containment,
    /// Mash distance -1/k * ln(2J / (1 + J)), capped at 1.0
    Distance,
}

/// Number of values present in both sorted lists.
pub fn intersection_size(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut shared) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            shared += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    shared
}

/// Jaccard index. Defined as 0.0 when both sets are empty.
pub fn jaccard(a: &[u64], b: &[u64]) -> f64 {
    let shared = intersection_size(a, b);
    let union = a.len() + b.len() - shared;
    if union == 0 {
        0.0
    } else {
        shared as f64 / union as f64
    }
}

/// Containment of `a` in `b`: the fraction of `a` that is also in `b`.
/// Defined as 0.0 when `a` is empty.
pub fn containment(a: &[u64], b: &[u64]) -> f64 {
    if a.is_empty() {
        0.0
    } else {
        intersection_size(a, b) as f64 / a.len() as f64
    }
}

/// Mash distance from a Jaccard index (Ondov et al. 2016). It approximates
/// the per-base mutation rate between two genomes. When nothing is shared
/// (J = 0) the formula gives infinity; like Mash we report 1.0 instead.
pub fn mash_distance(jaccard: f64, k: usize) -> f64 {
    if jaccard <= 0.0 {
        return 1.0;
    }
    if jaccard >= 1.0 {
        return 0.0; // avoids printing "-0.000000" (ln(1) = 0, negated)
    }
    let distance = -(2.0 * jaccard / (1.0 + jaccard)).ln() / k as f64;
    distance.clamp(0.0, 1.0)
}

/// Check that all sketches were built with the same k and scaled values;
/// comparing sketches built differently would give meaningless numbers.
pub fn check_compatible(sketches: &[Sketch]) -> Result<()> {
    if let Some(first) = sketches.first() {
        for sketch in sketches {
            if sketch.ksize != first.ksize || sketch.scaled != first.scaled {
                bail!(
                    "sketch '{}' has ksize={} scaled={}, but '{}' has ksize={} scaled={}",
                    sketch.name,
                    sketch.ksize,
                    sketch.scaled,
                    first.name,
                    first.ksize,
                    first.scaled
                );
            }
        }
    }
    Ok(())
}

/// Compute one metric between two sketches.
pub fn compare_pair(a: &Sketch, b: &Sketch, metric: Metric) -> f64 {
    match metric {
        Metric::Jaccard => jaccard(&a.hashes, &b.hashes),
        Metric::Containment => containment(&a.hashes, &b.hashes),
        Metric::Distance => mash_distance(jaccard(&a.hashes, &b.hashes), a.ksize),
    }
}

/// Full n x n matrix; rows are computed in parallel.
pub fn compare_all(sketches: &[Sketch], metric: Metric) -> Result<Vec<Vec<f64>>> {
    check_compatible(sketches)?;
    let matrix = sketches
        .par_iter()
        .map(|row| {
            sketches
                .iter()
                .map(|column| compare_pair(row, column, metric))
                .collect()
        })
        .collect();
    Ok(matrix)
}

/// Write the matrix as TSV. The first row and first column hold the names.
pub fn write_matrix(out: &mut impl Write, names: &[String], matrix: &[Vec<f64>]) -> Result<()> {
    write!(out, "name")?;
    for name in names {
        write!(out, "\t{name}")?;
    }
    writeln!(out)?;
    for (name, row) in names.iter().zip(matrix) {
        write!(out, "{name}")?;
        for value in row {
            write!(out, "\t{value:.6}")?;
        }
        writeln!(out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sketch(name: &str, ksize: usize, hashes: &[u64]) -> Sketch {
        Sketch {
            name: name.to_string(),
            ksize,
            scaled: 10,
            hashes: hashes.to_vec(),
        }
    }

    #[test]
    fn metrics_on_small_sets() {
        let a = [1, 2, 3];
        let b = [2, 3, 4];
        assert_eq!(intersection_size(&a, &b), 2);
        assert!((jaccard(&a, &b) - 0.5).abs() < 1e-12);
        assert!((containment(&a, &b) - 2.0 / 3.0).abs() < 1e-12);
        assert!((containment(&[2, 3], &a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn empty_sets() {
        assert_eq!(jaccard(&[], &[]), 0.0);
        assert_eq!(containment(&[], &[1]), 0.0);
        assert_eq!(jaccard(&[1], &[]), 0.0);
    }

    #[test]
    fn mash_distance_edge_cases() {
        assert_eq!(mash_distance(1.0, 21), 0.0);
        assert!(mash_distance(1.0, 21).is_sign_positive());
        assert_eq!(mash_distance(0.0, 21), 1.0);
        // J = 0.5 -> -(1/21) ln(2/3)
        let expected = -(2.0f64 / 3.0).ln() / 21.0;
        assert!((mash_distance(0.5, 21) - expected).abs() < 1e-12);
        // Distance decreases as similarity increases.
        assert!(mash_distance(0.2, 21) > mash_distance(0.8, 21));
    }

    #[test]
    fn matrix_is_symmetric_for_jaccard() {
        let sketches = vec![
            sketch("a", 21, &[1, 2, 3, 4]),
            sketch("b", 21, &[3, 4, 5]),
            sketch("c", 21, &[9]),
        ];
        let m = compare_all(&sketches, Metric::Jaccard).unwrap();
        for (i, row) in m.iter().enumerate() {
            assert_eq!(row[i], 1.0);
            for (j, value) in row.iter().enumerate() {
                assert_eq!(*value, m[j][i]);
            }
        }
        assert!((m[0][1] - 2.0 / 5.0).abs() < 1e-12);
        assert_eq!(m[0][2], 0.0);
    }

    #[test]
    fn incompatible_sketches_are_rejected() {
        let sketches = vec![sketch("a", 21, &[1]), sketch("b", 31, &[1])];
        assert!(compare_all(&sketches, Metric::Jaccard).is_err());
    }

    #[test]
    fn writes_matrix_tsv() {
        let names = vec!["a".to_string(), "b".to_string()];
        let matrix = vec![vec![1.0, 0.25], vec![0.25, 1.0]];
        let mut out = Vec::new();
        write_matrix(&mut out, &names, &matrix).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text,
            "name\ta\tb\na\t1.000000\t0.250000\nb\t0.250000\t1.000000\n"
        );
    }
}
