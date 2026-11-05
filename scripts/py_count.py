#!/usr/bin/env python3
"""Independent pure-Python canonical k-mer counter.

Used to check `kmer-sketch count` and as a speed baseline. It works on
plain strings only (no 2-bit encoding): the canonical k-mer is the
alphabetically smaller of the k-mer and its reverse complement. Any
character other than A/C/G/T (after upper-casing) splits the sequence, so
k-mers containing N are skipped. The output format matches the Rust tool:
a header line `kmer<TAB>count`, then k-mers in alphabetical order.

Usage: py_count.py -k 21 [--min-count 1] input.fa[.gz] [more inputs] > counts.tsv
"""
import argparse
import gzip
import re
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")
NOT_ACGT = re.compile("[^ACGT]+")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_sequences(path):
    """Yield sequences from a FASTA (multi-line ok) or FASTQ (4-line) file."""
    with open_text(path) as handle:
        first = handle.readline()
        if first.startswith(">"):
            parts = []
            for line in handle:
                if line.startswith(">"):
                    yield "".join(parts)
                    parts = []
                else:
                    parts.append(line.strip())
            yield "".join(parts)
        elif first.startswith("@"):
            while True:
                seq = handle.readline()
                handle.readline()  # '+' line
                handle.readline()  # quality line
                if not seq:
                    break
                yield seq.strip()
                if not handle.readline():  # next header
                    break
        elif first:
            raise ValueError(f"{path}: not FASTA or FASTQ")


def canonical_kmers(seq, k):
    """Yield the canonical k-mers of one sequence (same rules as below)."""
    for segment in NOT_ACGT.split(seq.upper()):
        n = len(segment)
        reverse = segment.translate(COMPLEMENT)[::-1]
        for i in range(n - k + 1):
            forward_kmer = segment[i : i + k]
            reverse_kmer = reverse[n - k - i : n - i]
            yield forward_kmer if forward_kmer <= reverse_kmer else reverse_kmer


def count_canonical(paths, k):
    # Same logic as canonical_kmers(), written inline because this loop is
    # the speed baseline and a generator adds overhead per k-mer.
    counts = {}
    for path in paths:
        for seq in read_sequences(path):
            for segment in NOT_ACGT.split(seq.upper()):
                n = len(segment)
                if n < k:
                    continue
                reverse = segment.translate(COMPLEMENT)[::-1]
                for i in range(n - k + 1):
                    forward_kmer = segment[i : i + k]
                    # reverse complement of segment[i:i+k] in the reversed string
                    reverse_kmer = reverse[n - k - i : n - i]
                    kmer = forward_kmer if forward_kmer <= reverse_kmer else reverse_kmer
                    counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-k", type=int, default=21)
    parser.add_argument("--min-count", type=int, default=1)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()

    counts = count_canonical(args.inputs, args.k)
    out = sys.stdout
    out.write("kmer\tcount\n")
    for kmer in sorted(counts):
        count = counts[kmer]
        if count >= args.min_count:
            out.write(f"{kmer}\t{count}\n")


if __name__ == "__main__":
    main()
