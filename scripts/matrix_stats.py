#!/usr/bin/env python3
"""Summarise a similarity matrix written by `kmer-sketch compare`.

Prints one TSV line with:
  cross_max      largest similarity between cells of DIFFERENT species
  within_nn_min  over all cells, the smallest "best same-species similarity"
                 (how weakly the worst-connected cell is tied to its species)
  within_mean    mean similarity between different cells of the same species
For single linkage at threshold t to recover every species exactly, the
condition cross_max < t <= within_nn_min is necessary (no link between
species, and no cell left without a same-species link) but not sufficient:
a species can still split into several groups that are not linked to each
other.
Usage: matrix_stats.py matrix.tsv [--header]
"""
import sys


def label(name):
    return name.split("__", 1)[0]


def main():
    path = sys.argv[1]
    with open(path) as handle:
        names = handle.readline().rstrip("\n").split("\t")[1:]
        rows = [list(map(float, line.rstrip("\n").split("\t")[1:])) for line in handle if line.strip()]
    n = len(names)
    cross_max = 0.0
    within = []
    nearest = []
    for i in range(n):
        best_same = 0.0
        for j in range(n):
            if i == j:
                continue
            s = max(rows[i][j], rows[j][i])
            if label(names[i]) == label(names[j]):
                within.append(s)
                best_same = max(best_same, s)
            else:
                cross_max = max(cross_max, s)
        nearest.append(best_same)
    within_mean = sum(within) / len(within) if within else 0.0
    if "--header" in sys.argv:
        print("cross_max\twithin_nn_min\twithin_mean")
    print(f"{cross_max:.6f}\t{min(nearest):.6f}\t{within_mean:.6f}")


if __name__ == "__main__":
    main()
