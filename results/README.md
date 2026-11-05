# Results

Raw outputs of the scripts in `scripts/`. Every number in the top-level
README comes from one of these files. Machine: `machine.txt`. Genomes and
lengths: `genomes.tsv`.

## Main pipeline (`pipeline/`, from `scripts/run_pipeline.sh`)

40 cells (10 per species), 150 bp reads, coverage drawn uniformly in 2x-30x
(realised 2.09x-28.10x, `pipeline/cells.tsv`), 0.1% substitution errors,
seed 42; 18,125,004 reads, 5,985,928,540 bytes of FASTQ (`pipeline/sizes.tsv`).
k=21, scaled=1000; median sketch 5,632 hashes (`sweep/separation.tsv`, same
setting and seed).

| metric | thresholds | clusters | ARI | species recovered |
|---|---|---|---|---|
| Jaccard | 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5 | 4 | 1.0000 | 4/4 |
| Jaccard | 0.7 | 19 | 0.4947 | 0/4 (cells split, every cluster pure) |
| containment | 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9 | 4 | 1.0000 | 4/4 |
| Mash distance | 0.05, 0.1, 0.2 | 4 | 1.0000 | 4/4 |
| Mash distance | 0.01 | 16 | 0.5900 | 0/4 |

Source: `pipeline/cluster_summary.tsv`. Separation (`pipeline/matrix_stats.tsv`):
largest cross-species Jaccard 0.000611, smallest best same-species Jaccard
0.607961. Wall times with 8 threads (`pipeline/timings.tsv`): simulate 33.13 s,
sketch 2.85 s, each compare 0.02 s.

## Sweep (`sweep/`, from `scripts/sweep.sh` and `scripts/sweep_fine.sh`)

40 cells per setting, seed 42, k=21. "Working thresholds" are the thresholds
on the tried grid (Jaccard 0.0005 to 0.7, containment 0.001 to 0.9) that give
ARI = 1. The threshold is chosen after looking at the answer, so these are
best cases, not what an unsupervised run would pick. Each setting is one
simulation (seed 42); at the lowest coverages a cell keeps only ~40 hashes,
so a different seed can change these results.

| setting | scaled | median hashes/cell | max cross-species J | min same-species nearest J | best ARI | working Jaccard thresholds |
|---|---|---|---|---|---|---|
| 2-30x, 0.1% errors | 1000 | 5,632 | 0.0006 | 0.608 | 1.0 | 0.001-0.5 |
| 2-30x, 1% errors | 1000 | 14,774 | 0.0005 | 0.158 | 1.0 | 0.001-0.1 |
| 2-30x, 5% errors | 1000 | 40,136 | 0.0005 | 0.071 | 1.0 | 0.001-0.05 |
| 0.5-2x | 1000 | 2,602 | 0.0009 | 0.294 | 1.0 | 0.001-0.2 |
| 0.1-0.5x | 1000 | 911 | 0.0016 | 0.073 | 1.0 | 0.002-0.05 |
| 0.02-0.1x | 1000 | 217 | 0.0036 | 0.017 | 1.0 | 0.005-0.01 |
| 0.005-0.02x | 1000 | 43 | 0.0000 | 0.009 | 1.0 | 0.0005-0.005 |
| 0.005-0.02x | 100 | 439 | 0.0027 | 0.005 | 1.0 | 0.005 only |
| + E. coli strain (1% subst.), 2-30x | 1000 | 5,983 | 0.518 | 0.563 | 1.0 | 0.52-0.56 (fine scan) |
| + E. coli strain (1% subst.), 0.1-0.5x | 1000 | 919 | 0.150 | 0.078 | 0.764 (J), 0.898 (containment, fine scan) | none |

Sources: `sweep/sweep_best.tsv` (best ARI and thresholds on the coarse grid),
`sweep/separation.tsv` (hash counts and separation), `sweep/strains_fine.tsv`
and `sweep/strains_window.tsv` (0.01-step scan for the strain settings).
The artificial strain has 45,846 substitutions (rate 0.009877,
`sweep/ecolimut_substitutions.tsv`).

- Errors and low coverage shrink same-species similarity a lot, but the four
  species share almost nothing, so a low threshold still separates them.
- With the close strain at 0.1-0.5x, cross-strain similarity (up to 0.150) is
  higher than the weakest same-strain link (0.078), so no threshold works:
  either the two E. coli strains merge (Jaccard, ARI 0.764, 4 clusters) or the
  clusters stay pure but fragmented (containment 0.29: 9 clusters, ARI 0.898).

## Correctness (`count_check.tsv`, `sketch_check/`)

| check | result | source |
|---|---|---|
| `count -k 21` on the E. coli genome vs `scripts/py_count.py` | byte-identical, 4,543,891 distinct canonical 21-mers, 1/2/4/8 threads | `count_check.tsv` |
| `count -k 21` on 50x simulated E. coli reads (1,547,217 reads of 150 bp, 0.1% errors) vs Python | byte-identical, 8,697,834 distinct canonical 21-mers, 1/2/4/8 threads | `count_check.tsv` |
| Rust sketch hashes vs Python FracMinHash (same hash and threshold) | identical for 6 cells (scaled 1000) and 5 genomes (scaled 100 and 1000) | `sketch_check/*_hashes*.tsv` |
| sketch vs exact Jaccard, 15 pairs of 6 cells (4 E. coli at 2.1x-17.2x, 2 S. aureus), scaled 1000 | largest error 0.0136 (same-species pairs, exact J 0.59-0.73); cross-species exact J below 0.00021 | `sketch_check/cells_pairs.tsv`, `cells_used.tsv` |
| sketch vs exact Jaccard, 10 pairs of E. coli, 3 mutated copies and B. subtilis, scaled 100 | largest error 0.0010 (related pairs, exact J 0.16-0.96) | `sketch_check/genomes_pairs_100.tsv` |
| same, scaled 1000 | largest error 0.0121 | `sketch_check/genomes_pairs_1000.tsv` |
| Mash distance vs realised substitution rate | true 0.000999 / 0.009877 / 0.050007; scaled 100: 0.000993 / 0.009979 / 0.051248; scaled 1000: 0.001013 / 0.009951 / 0.049056 | `sketch_check/mash_vs_true.tsv`, `mutants.tsv` |
| simulator refactor (record picking without per-read allocation) | same-seed outputs byte-identical to the earlier runs | `simulate_determinism.txt` |
| fmt, clippy `-D warnings`, tests | all pass | `quality.txt` |

## Speed of `count` (`bench_summary.tsv`, raw runs in `bench_raw.tsv`)

k=21, median of 3 runs, including reading the input and writing the sorted
TSV to a file. Peak memory is the maximum resident set size of the process.

| input | tool | threads | median wall time | peak memory |
|---|---|---|---|---|
| E. coli genome (4.64 Mbp, 4.54 M distinct k-mers) | py_count.py | 1 | 10.23 s | 464 MB |
| | kmer-sketch | 1 | 1.15 s | 219 MB |
| | kmer-sketch | 2 | 0.98 s | 205 MB |
| | kmer-sketch | 4 | 0.87 s | 208 MB |
| | kmer-sketch | 8 | 0.87 s | 233 MB |
| 50x reads (232 Mbp, 506 MB FASTQ, 8.70 M distinct k-mers) | py_count.py | 1 | 125.07 s | 898 MB |
| | kmer-sketch | 1 | 7.17 s | 425 MB |
| | kmer-sketch | 2 | 4.69 s | 426 MB |
| | kmer-sketch | 4 | 3.68 s | 430 MB |
| | kmer-sketch | 8 | 3.61 s | 434 MB |

The benchmark was run before the later input-validation changes. The only
change on the counting path is an `is_file()` check and the empty-file case
when a file is opened; it was not re-timed.
