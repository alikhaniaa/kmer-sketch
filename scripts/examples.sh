#!/usr/bin/env bash
# Run the example commands shown in the README and save their output to
# results/examples.txt. Needs data/genomes (fetch_genomes.sh) and the
# pipeline outputs (run_pipeline.sh).
set -euo pipefail

bin=target/release/kmer-sketch
work=data/examples
mkdir -p "$work"
cargo build --release --quiet

show() {
    echo "\$ $*"
    eval "$@"
    echo
}

{
    show "$bin count -k 21 --min-count 20 -o $work/ecoli_min20.tsv data/genomes/ecoli.fna"
    show "head -n 4 $work/ecoli_min20.tsv"
    show "wc -l < $work/ecoli_min20.tsv"
    show "$bin sketch -k 21 --scaled 1000 -o $work/genome_sketches data/genomes/*.fna"
    show "head -n 9 $work/genome_sketches/ecoli.sketch"
    show "$bin compare --metric jaccard $work/genome_sketches/*.sketch"
    show "$bin compare --metric distance $work/genome_sketches/*.sketch"
    show "$bin compare --metric containment $work/genome_sketches/ecoli.sketch $work/genome_sketches/paeruginosa.sketch"
    show "head -n 5 data/pipeline/cells/cells.tsv"
    show "$bin cluster --threshold 0.1 results/pipeline/jaccard.tsv"
} > results/examples.txt 2>&1
cat results/examples.txt
