# kmer-sketch

A small Rust command-line tool for k-mer work on bacterial reads, written as a
learning project. It counts canonical k-mers, builds FracMinHash sketches,
compares sketches (Jaccard, containment, Mash distance), simulates
single-cell read sets from reference genomes, and clusters the simulated
cells to see whether species can be recovered from k-mer sketches alone.

The single-cell part is a toy version of one step of a microbial single-cell
pipeline: every "cell" is a FASTQ file of reads from one genome, with its own
coverage, and the goal is to group cells by species from k-mer sketches
without assembly.

## Build and run

Needs Rust 1.88 or newer (developed with rustc 1.99.1).

```
cargo build --release
cargo test
./target/release/kmer-sketch --help
```

| command    | what it does |
|------------|--------------|
| `count`    | canonical k-mer counts (k <= 31) from FASTA/FASTQ, plain or `.gz`; sorted TSV; `--min-count` |
| `sketch`   | one FracMinHash sketch per input file (`--ksize`, default 21; `--scaled`, default 1000) |
| `compare`  | pairwise matrix of `jaccard`, `containment` or Mash `distance` between sketches |
| `simulate` | reads from genomes as "cells": random substrings from both strands, per-cell coverage drawn uniformly in a range, substitution errors, fixed seed; one FASTQ per cell |
| `cluster`  | single-linkage clustering of cells on a matrix with a threshold; reports cluster composition against the species label in the cell names |

`-t/--threads N` sets the number of worker threads (default: all cores).
Parallelism (rayon) is over cells in `simulate`, over files in `sketch`, over
matrix rows in `compare`, and over chunks of reads and hash-map shards in
`count`. An empty input file (also an empty `.gz`) is read as zero records;
its sketch has 0 hashes (with a warning) and gives Jaccard 0 and Mash
distance 1 against everything, itself included.

## Example commands

Outputs are copied from `results/examples.txt` (made by `scripts/examples.sh`;
directory names shortened, `kmer-sketch` stands for `target/release/kmer-sketch`).

```
$ kmer-sketch count -k 21 --min-count 20 -o ecoli_min20.tsv data/genomes/ecoli.fna
$ head -n 4 ecoli_min20.tsv
kmer	count
AAACGCCTTATCCGGCCTACA	21
AACGCCTTATCCGGCCTACAA	45
AACGCCTTATCCGGCCTACGG	25

$ kmer-sketch sketch -k 21 --scaled 1000 -o genome_sketches data/genomes/*.fna
wrote 4 sketches to genome_sketches
$ head -n 9 genome_sketches/ecoli.sketch
# kmer-sketch FracMinHash v1
name	ecoli
ksize	21
scaled	1000
hash	splitmix64
count	4583
hashes
348139822534
3829463084548
$ kmer-sketch compare --metric distance genome_sketches/*.sketch
name	bsubtilis	ecoli	paeruginosa	saureus
bsubtilis	0.000000	1.000000	0.407357	1.000000
ecoli	1.000000	0.000000	0.323871	1.000000
paeruginosa	0.407357	0.323871	0.000000	1.000000
saureus	1.000000	1.000000	1.000000	0.000000

$ kmer-sketch simulate --cells 40 --min-coverage 2 --max-coverage 30 \
      --error-rate 0.001 --read-length 150 --seed 42 -o cells data/genomes/*.fna
$ kmer-sketch sketch -k 21 --scaled 1000 -o sketches cells/*.fastq
$ kmer-sketch compare --metric jaccard -o jaccard.tsv sketches/*.sketch
$ kmer-sketch cluster --threshold 0.1 jaccard.tsv
# cells=40 clusters=4 labels=4 ARI=1.0000 purity=1.0000 species_recovered=4/4
cluster	size	majority	purity	composition
0	10	bsubtilis	1.0000	bsubtilis:10
1	10	ecoli	1.0000	ecoli:10
2	10	paeruginosa	1.0000	paeruginosa:10
3	10	saureus	1.0000	saureus:10
```

A Mash distance of 1.0 means no shared hashes. In a containment matrix each
row gives the fraction of the row sample found in the column sample.

## Reproducing the results

```
bash scripts/fetch_genomes.sh   # 4 genomes from NCBI efetch into data/genomes/ (not committed)
bash scripts/run_pipeline.sh    # 40 cells -> sketch -> compare -> cluster   -> results/pipeline/
bash scripts/check_sketch.sh    # sketches vs exact Python k-mer sets         -> results/sketch_check/
bash scripts/sweep.sh           # coverage / error / close-strain sweep       -> results/sweep/
bash scripts/sweep_fine.sh      # fine threshold scan for the strain settings -> results/sweep/
bash scripts/bench.sh           # count: Rust vs pure Python, timing + diff   -> results/bench_*.tsv
bash scripts/examples.sh        # the example commands above                  -> results/examples.txt
bash scripts/quality.sh         # fmt, clippy, tests                          -> results/quality.txt
```

The Python helpers (`python3 scripts/py_count.py`, `python3 scripts/py_sketch_check.py`,
...) use only the standard library. Genomes: E. coli K-12 MG1655 (NC_000913.3),
B. subtilis 168 (NC_000964.3), S. aureus NCTC 8325 (NC_007795.1),
P. aeruginosa PAO1 (NC_002516.2); lengths in `results/genomes.tsv`.

## Results

Machine (`results/machine.txt`): Intel Core i7-6700 @ 3.40 GHz, 4 cores / 8
logical CPUs, 23.5 GiB RAM, Linux under WSL2, rustc 1.99.1, Python 3.14.8.
Full tables, with the file each number comes from, are in `results/README.md`.

**Species recovery.** 40 cells (10 per species), coverage drawn in 2x-30x,
0.1% substitution errors, seed 42, k=21, scaled=1000: 4 clusters, ARI 1.0,
4/4 species recovered for every Jaccard threshold from 0.001 to 0.5
(`results/pipeline/cluster_summary.tsv`). This is easy because the species are
far apart: the largest Jaccard between cells of different species is 0.0006,
and every cell has a same-species neighbour with Jaccard at least 0.61.

**Where it breaks** (`results/sweep/`). With 1% or 5% read errors and with
coverage down to 0.005-0.02x the four species are still recovered, but the
working threshold window narrows (at 0.005-0.02x and scaled=100 only one
tried value works). Adding an artificial close strain (E. coli with 1% random
substitutions) is the real failure: at 2-30x it separates only for Jaccard
thresholds 0.52-0.56, and at 0.1-0.5x no threshold works (best ARI 0.764 with
Jaccard, 0.898 with containment). Thresholds were picked after seeing the
answer, so these are best cases, not what an unsupervised run would pick.
Each setting is one simulation (seed 42); at the lowest coverages a cell
keeps only ~40 hashes, so a different seed can change these results.

**Correctness** (`results/count_check.tsv`, `results/sketch_check/`).
`count` output is byte-identical to an independent pure-Python counter on the
E. coli genome (4,543,891 distinct 21-mers) and on a 50x read set (8,697,834)
at 1, 2, 4 and 8 threads. Rust sketches contain exactly the same hashes as a
Python FracMinHash for 6 cells and 5 genomes. Sketch Jaccard is within 0.0136
of exact Jaccard for cell pairs at scaled=1000, and within 0.0010 for E. coli
vs mutated copies at scaled=100. Mash distance to copies with 0.0999%,
0.9877% and 5.0007% substitutions is 0.000993, 0.009979 and 0.051248
(scaled=100). All CI checks pass (`results/quality.txt`).

**Speed of `count`** (k=21, median of 3 runs, including writing the sorted
TSV; `results/bench_summary.tsv`):

| input | py_count.py | kmer-sketch, 1 thread | kmer-sketch, 8 threads |
|---|---|---|---|
| E. coli genome (4.64 Mbp) | 10.23 s | 1.15 s | 0.87 s |
| 50x simulated reads (232 Mbp) | 125.07 s | 7.17 s | 3.61 s |

For the genome (one FASTA record) k-mer extraction runs on one thread,
because `count` parallelises over reads; only the shard-map updates and the
sort are parallel. On the read set, scaling stops at about 4 threads: reading
the file and writing the sorted output are serial, and the CPU has 4
physical cores. Pure Python is a weak baseline; no optimised k-mer counter
was run for comparison.

## Design notes

- **Canonical k-mers.** A read can come from either DNA strand, so a k-mer
  and its reverse complement are the same biological sequence. Counting the
  smaller of the two (the canonical form) makes counts and sketches
  independent of strand.
- **2-bit encoding.** A, C, G, T are stored as 0, 1, 2, 3, so a k-mer with
  k <= 31 fits in one `u64` (62 bits; k <= 31 also keeps the mask
  `(1 << 2k) - 1` from overflowing). Integers are cheaper to hash and compare
  than strings. The first base sits in the most significant bits, so sorting
  the integers gives alphabetical order, which is how `count` sorts its
  output. Complement is `3 - code`.
- **Rolling update.** When the window moves by one base, the forward k-mer is
  shifted left and the new base added at the bottom; the reverse complement
  is shifted right and the complement of the new base added at the top. Each
  base costs O(1). Any non-ACGT base (for example N) resets the window, so
  k-mers containing N are skipped.
- **Counting in parallel.** k-mers are split into 64 shards by hash. Threads
  first sort k-mers from chunks of reads into per-shard buckets, then each
  shard's hash map (`rustc-hash` FxHashMap) is filled by a single thread, so
  no locks are needed. Reads are processed in batches of about 16 Mbases.
- **FracMinHash.** Keep every k-mer whose hash is `<= u64::MAX / scaled`.
  Unlike a fixed-size MinHash, the sketch grows with the number of distinct
  k-mers, so it can estimate containment between sets of very different size
  (a 2x cell against a 30x cell, or a cell against a genome), and any two
  sketches with the same k and scaled can be compared directly. The hash is
  SplitMix64, fixed and unseeded, so sketches from different runs are
  comparable. The sketch files are this tool's own format and are not
  interchangeable with other tools' sketches.
- **Sketch file format.** Plain text: a `#` comment line, `key<TAB>value`
  lines for `name`, `ksize`, `scaled`, `hash`, `count`, then a `hashes` line
  and the sorted hash values in decimal, one per line (see `src/sketch.rs`).
- **Mash distance** is `-1/k * ln(2J / (1 + J))`, reported as 1.0 when J = 0.
- **Clustering.** Single linkage with a threshold is connected components of
  the graph "similarity >= threshold" (or "distance <= threshold"), found with
  union-find. For the asymmetric containment matrix the larger of the two
  directions is used. The true label is the part of the cell name before
  `__` (`ecoli__cell001` -> `ecoli`). The report gives the adjusted Rand
  index, purity, and how many species form exactly one pure cluster.

## Limitations

- All reads are simulated: uniform random positions, substitutions only (no
  indels, no quality-dependent errors, no GC or amplification bias, no
  chimeras, no contamination between cells, constant quality string).
  Real single-cell bacterial data is much messier.
- No real barcode handling or demultiplexing: each cell is already its own file.
- No assembly, no taxonomic classification against a database, no abundance
  tracking in sketches.
- The clustering threshold is picked by hand; the sweep shows the right value
  depends on coverage and on how related the organisms are.
- The "close strain" is an artificial genome with uniform random
  substitutions, not a real strain (real strains differ by gene content,
  recombination and mobile elements).
- `count` keeps all distinct k-mers in memory, which is fine for bacterial
  genomes and read sets of this size but not for large metagenomes.
- This is a teaching-sized tool. It is not a replacement for sourmash,
  jellyfish, Mash or single-cell genomics pipelines.

## License

MIT, see `LICENSE`.
