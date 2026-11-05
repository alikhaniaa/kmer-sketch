//! Single-linkage clustering of cells from a comparison matrix, and a report
//! that compares the clusters with the true species labels.
//!
//! Single linkage with a threshold is the same as: draw an edge between two
//! cells when they are similar enough, then every connected group of cells is
//! one cluster. We find the connected groups with a union-find structure.

use anyhow::{bail, Context, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, Write};

/// A square matrix with a name for each row/column.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub names: Vec<String>,
    pub values: Vec<Vec<f64>>,
}

/// Read a TSV matrix as written by `compare` (first row and column = names).
pub fn read_matrix(input: impl BufRead) -> Result<Matrix> {
    let mut lines = input.lines();
    let header = lines.next().context("matrix file is empty")??;
    let names: Vec<String> = header.split('\t').skip(1).map(str::to_string).collect();

    let mut values = Vec::new();
    for (row_index, line) in lines.enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let row_name = fields.next().unwrap_or_default();
        if names.get(row_index).map(String::as_str) != Some(row_name) {
            bail!(
                "row {} is named '{row_name}', expected the same order as the header",
                row_index + 1
            );
        }
        let row: Vec<f64> = fields
            .map(|field| {
                field
                    .parse::<f64>()
                    .with_context(|| format!("bad number '{field}'"))
            })
            .collect::<Result<_>>()?;
        if row.len() != names.len() {
            bail!(
                "row '{row_name}' has {} values, expected {}",
                row.len(),
                names.len()
            );
        }
        values.push(row);
    }
    if values.len() != names.len() {
        bail!(
            "matrix has {} rows but {} columns",
            values.len(),
            names.len()
        );
    }
    Ok(Matrix { names, values })
}

/// Union-find (disjoint set) structure.
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }

    /// Representative of the set containing `x`.
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            // Path halving: point x at its grandparent to keep trees flat.
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let root_a = self.find(a);
        let root_b = self.find(b);
        if root_a != root_b {
            self.parent[root_b] = root_a;
        }
    }
}

/// How to read the matrix values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Linkage {
    /// Values are similarities: link when value >= threshold.
    Similarity,
    /// Values are distances: link when value <= threshold.
    Distance,
}

/// Should cells `i` and `j` be joined?
fn is_linked(values: &[Vec<f64>], i: usize, j: usize, threshold: f64, linkage: Linkage) -> bool {
    match linkage {
        Linkage::Similarity => values[i][j].max(values[j][i]) >= threshold,
        Linkage::Distance => values[i][j].min(values[j][i]) <= threshold,
    }
}

/// Single-linkage clustering. Returns a cluster id for every row; ids are
/// 0, 1, 2, ... in order of first appearance.
///
/// The matrix may be asymmetric (containment is). For similarities we use
/// the larger of the two directions, for distances the smaller one.
pub fn single_linkage(values: &[Vec<f64>], threshold: f64, linkage: Linkage) -> Vec<usize> {
    let n = values.len();
    let mut sets = UnionFind::new(n);
    for i in 0..n {
        for j in (i + 1)..n {
            if is_linked(values, i, j, threshold, linkage) {
                sets.union(i, j);
            }
        }
    }

    // Renumber the set representatives as 0, 1, 2, ...
    let mut ids: BTreeMap<usize, usize> = BTreeMap::new();
    let mut clusters = Vec::with_capacity(n);
    for i in 0..n {
        let root = sets.find(i);
        let next_id = ids.len();
        clusters.push(*ids.entry(root).or_insert(next_id));
    }
    clusters
}

/// The true label of a cell is the part of its name before `separator`
/// (e.g. `ecoli__cell007` -> `ecoli`). Names without the separator are
/// their own label.
pub fn label_from_name(name: &str, separator: &str) -> String {
    match name.split_once(separator) {
        Some((label, _)) => label.to_string(),
        None => name.to_string(),
    }
}

/// n choose 2, as f64.
fn pairs(n: usize) -> f64 {
    let n = n as f64;
    n * (n - 1.0) / 2.0
}

/// Adjusted Rand index between two labelings of the same items.
/// 1.0 = identical partitions, around 0.0 = no better than chance.
/// If the formula is 0/0 (e.g. fewer than two items) we return 1.0 when the
/// partitions agree on every pair and 0.0 otherwise.
pub fn adjusted_rand_index(truth: &[String], clusters: &[usize]) -> f64 {
    assert_eq!(truth.len(), clusters.len());
    let mut table: BTreeMap<(&str, usize), usize> = BTreeMap::new();
    let mut label_sizes: BTreeMap<&str, usize> = BTreeMap::new();
    let mut cluster_sizes: BTreeMap<usize, usize> = BTreeMap::new();
    for (label, &cluster) in truth.iter().zip(clusters) {
        *table.entry((label.as_str(), cluster)).or_insert(0) += 1;
        *label_sizes.entry(label.as_str()).or_insert(0) += 1;
        *cluster_sizes.entry(cluster).or_insert(0) += 1;
    }

    let index: f64 = table.values().map(|&n| pairs(n)).sum();
    let label_pairs: f64 = label_sizes.values().map(|&n| pairs(n)).sum();
    let cluster_pairs: f64 = cluster_sizes.values().map(|&n| pairs(n)).sum();
    let total_pairs = pairs(truth.len());

    let expected = if total_pairs > 0.0 {
        label_pairs * cluster_pairs / total_pairs
    } else {
        0.0
    };
    let maximum = (label_pairs + cluster_pairs) / 2.0;
    if (maximum - expected).abs() < 1e-12 {
        return if (index - expected).abs() < 1e-12 {
            1.0
        } else {
            0.0
        };
    }
    (index - expected) / (maximum - expected)
}

/// Composition of one cluster.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterSummary {
    pub id: usize,
    pub size: usize,
    /// Label -> number of cells with that label in this cluster.
    pub composition: BTreeMap<String, usize>,
}

impl ClusterSummary {
    /// Most common label (ties broken alphabetically) and its count.
    pub fn majority(&self) -> (&str, usize) {
        let mut best = ("", 0);
        for (label, &count) in &self.composition {
            if count > best.1 {
                best = (label.as_str(), count);
            }
        }
        best
    }
}

/// Overall comparison of clusters with true labels.
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    pub clusters: Vec<ClusterSummary>,
    pub n_cells: usize,
    pub n_labels: usize,
    pub ari: f64,
    /// Fraction of cells whose label is the majority label of their cluster.
    pub purity: f64,
    /// Labels whose cells form exactly one cluster with nothing else in it.
    pub recovered_labels: usize,
}

/// Build the report from per-cell labels and cluster ids.
pub fn make_report(labels: &[String], clusters: &[usize]) -> Report {
    let n_clusters = clusters.iter().max().map_or(0, |&max| max + 1);
    let mut summaries: Vec<ClusterSummary> = (0..n_clusters)
        .map(|id| ClusterSummary {
            id,
            size: 0,
            composition: BTreeMap::new(),
        })
        .collect();
    for (label, &cluster) in labels.iter().zip(clusters) {
        let summary = &mut summaries[cluster];
        summary.size += 1;
        *summary.composition.entry(label.clone()).or_insert(0) += 1;
    }

    let mut label_sizes: BTreeMap<&str, usize> = BTreeMap::new();
    for label in labels {
        *label_sizes.entry(label.as_str()).or_insert(0) += 1;
    }

    let majority_total: usize = summaries.iter().map(|s| s.majority().1).sum();
    // A label is recovered if some cluster contains all of its cells and
    // nothing else.
    let recovered_labels = summaries
        .iter()
        .filter(|s| s.composition.len() == 1)
        .filter(|s| {
            let (label, count) = s.majority();
            label_sizes.get(label) == Some(&count)
        })
        .count();

    Report {
        n_cells: labels.len(),
        n_labels: label_sizes.len(),
        ari: adjusted_rand_index(labels, clusters),
        purity: if labels.is_empty() {
            0.0
        } else {
            majority_total as f64 / labels.len() as f64
        },
        recovered_labels,
        clusters: summaries,
    }
}

/// Write the report: one summary comment line, then one TSV row per cluster.
pub fn write_report(out: &mut impl Write, report: &Report) -> Result<()> {
    writeln!(
        out,
        "# cells={} clusters={} labels={} ARI={:.4} purity={:.4} species_recovered={}/{}",
        report.n_cells,
        report.clusters.len(),
        report.n_labels,
        report.ari,
        report.purity,
        report.recovered_labels,
        report.n_labels
    )?;
    writeln!(out, "cluster\tsize\tmajority\tpurity\tcomposition")?;
    for summary in &report.clusters {
        let (majority, count) = summary.majority();
        let composition: Vec<String> = summary
            .composition
            .iter()
            .map(|(label, n)| format!("{label}:{n}"))
            .collect();
        writeln!(
            out,
            "{}\t{}\t{}\t{:.4}\t{}",
            summary.id,
            summary.size,
            majority,
            count as f64 / summary.size as f64,
            composition.join(",")
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(items: &[&str]) -> Vec<String> {
        items.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn reads_matrix() {
        let text = "name\ta\tb\na\t1.0\t0.2\nb\t0.2\t1.0\n";
        let m = read_matrix(text.as_bytes()).unwrap();
        assert_eq!(m.names, strings(&["a", "b"]));
        assert_eq!(m.values, vec![vec![1.0, 0.2], vec![0.2, 1.0]]);
    }

    #[test]
    fn rejects_ragged_matrix() {
        let text = "name\ta\tb\na\t1.0\nb\t0.2\t1.0\n";
        assert!(read_matrix(text.as_bytes()).is_err());
    }

    #[test]
    fn single_linkage_chains_through_neighbours() {
        // a-b similar, b-c similar, a-c not: single linkage still joins all three.
        let values = vec![
            vec![1.0, 0.5, 0.0, 0.0],
            vec![0.5, 1.0, 0.5, 0.0],
            vec![0.0, 0.5, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        assert_eq!(
            single_linkage(&values, 0.4, Linkage::Similarity),
            vec![0, 0, 0, 1]
        );
        assert_eq!(
            single_linkage(&values, 0.6, Linkage::Similarity),
            vec![0, 1, 2, 3]
        );
    }

    #[test]
    fn distance_linkage_and_asymmetric_values() {
        let values = vec![vec![0.0, 0.9], vec![0.05, 0.0]];
        assert_eq!(single_linkage(&values, 0.1, Linkage::Distance), vec![0, 0]);
        let similarities = vec![vec![1.0, 0.1], vec![0.8, 1.0]];
        assert_eq!(
            single_linkage(&similarities, 0.5, Linkage::Similarity),
            vec![0, 0]
        );
    }

    #[test]
    fn labels_from_names() {
        assert_eq!(label_from_name("ecoli__cell007", "__"), "ecoli");
        assert_eq!(label_from_name("plain", "__"), "plain");
    }

    #[test]
    fn ari_values() {
        let truth = strings(&["x", "x", "y", "y"]);
        assert!((adjusted_rand_index(&truth, &[0, 0, 1, 1]) - 1.0).abs() < 1e-12);
        // Cluster ids are arbitrary: swapping them changes nothing.
        assert!((adjusted_rand_index(&truth, &[1, 1, 0, 0]) - 1.0).abs() < 1e-12);
        // Worked by hand: index=1, expected=1, max=2.5 -> ARI = 0.
        assert!(adjusted_rand_index(&truth, &[0, 0, 0, 1]).abs() < 1e-12);
    }

    #[test]
    fn report_counts_recovered_species() {
        let labels = strings(&["x", "x", "y", "y", "z"]);
        // x is recovered; y is split; z sits alone and is recovered.
        let report = make_report(&labels, &[0, 0, 1, 2, 3]);
        assert_eq!(report.recovered_labels, 2);
        assert_eq!(report.clusters.len(), 4);
        assert!((report.purity - 1.0).abs() < 1e-12);

        // Everything merged: nothing recovered (the cluster is mixed).
        let merged = make_report(&labels, &[0, 0, 0, 0, 0]);
        assert_eq!(merged.recovered_labels, 0);
        assert!((merged.purity - 0.4).abs() < 1e-12);
    }
}
