#!/usr/bin/env python3
"""Check `kmer-sketch sketch` and `compare` against exact k-mer sets.

For every input file this script
  1. builds the EXACT set of canonical k-mers (pure Python, strings),
  2. builds its own FracMinHash sketch from that set with the same hash
     (SplitMix64 of the 2-bit k-mer value) and threshold (2^64-1) // scaled,
  3. checks that the Rust sketch file holds exactly the same hashes.
For every pair of inputs it reports the exact Jaccard index, the Python
sketch estimate and the Rust estimate (read from the `compare` matrix).

Usage:
  py_sketch_check.py -k 21 --scaled 1000 --rust-sketch-dir DIR \
      --rust-matrix jaccard.tsv --pairs-out pairs.tsv --hashes-out hashes.tsv \
      input1.fastq input2.fastq ...
"""
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from py_count import canonical_kmers, read_sequences  # noqa: E402

MASK64 = (1 << 64) - 1
TO_DIGITS = str.maketrans("ACGT", "0123")


def splitmix64(x):
    z = (x + 0x9E3779B97F4A7C15) & MASK64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return z ^ (z >> 31)


def sample_name(path):
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in (".fastq", ".fq", ".fasta", ".fa", ".fna"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def exact_kmer_set(path, k):
    kmers = set()
    for seq in read_sequences(path):
        kmers.update(canonical_kmers(seq, k))
    return kmers


def python_sketch(kmers, scaled):
    max_hash = MASK64 // scaled
    kept = set()
    for kmer in kmers:
        h = splitmix64(int(kmer.translate(TO_DIGITS), 4))
        if h <= max_hash:
            kept.add(h)
    return kept


def read_rust_sketch(path):
    hashes = []
    in_hashes = False
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if in_hashes:
                hashes.append(int(line))
            elif line == "hashes":
                in_hashes = True
    return set(hashes)


def read_matrix(path):
    with open(path) as handle:
        names = handle.readline().rstrip("\n").split("\t")[1:]
        values = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            for name, value in zip(names, fields[1:]):
                values[(fields[0], name)] = float(value)
    return values


def jaccard(a, b):
    union = len(a | b)
    return len(a & b) / union if union else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--scaled", type=int, default=1000)
    parser.add_argument("--rust-sketch-dir", required=True)
    parser.add_argument("--rust-matrix", required=True)
    parser.add_argument("--pairs-out", required=True)
    parser.add_argument("--hashes-out", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    names = [sample_name(p) for p in args.inputs]
    exact = {}
    sketches = {}
    with open(args.hashes_out, "w") as out:
        out.write("name\texact_kmers\tpython_hashes\trust_hashes\tidentical\n")
        for name, path in zip(names, args.inputs):
            exact[name] = exact_kmer_set(path, args.k)
            sketches[name] = python_sketch(exact[name], args.scaled)
            rust = read_rust_sketch(os.path.join(args.rust_sketch_dir, name + ".sketch"))
            same = rust == sketches[name]
            out.write(f"{name}\t{len(exact[name])}\t{len(sketches[name])}\t{len(rust)}\t{same}\n")
            print(f"{name}: {len(exact[name])} k-mers, hashes identical: {same}", file=sys.stderr)

    rust_matrix = read_matrix(args.rust_matrix)
    with open(args.pairs_out, "w") as out:
        out.write("a\tb\texact_jaccard\tpython_sketch_jaccard\trust_sketch_jaccard\tabs_error\n")
        for a, b in itertools.combinations(names, 2):
            j_exact = jaccard(exact[a], exact[b])
            j_python = jaccard(sketches[a], sketches[b])
            j_rust = rust_matrix[(a, b)]
            out.write(
                f"{a}\t{b}\t{j_exact:.6f}\t{j_python:.6f}\t{j_rust:.6f}\t{abs(j_rust - j_exact):.6f}\n"
            )


if __name__ == "__main__":
    main()
