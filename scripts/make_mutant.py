#!/usr/bin/env python3
"""Write a copy of a FASTA genome with random substitutions.

Each A/C/G/T position is replaced, with probability RATE, by one of the
three other bases. Used to make artificial "strains" with a known
divergence, for checking Mash distances and for harder clustering tests.

Usage: make_mutant.py input.fna output.fna RATE SEED
Prints the number of substituted positions and the realised rate.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from py_count import read_sequences  # noqa: E402

OTHERS = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}


def main():
    source, target, rate, seed = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    rng = random.Random(seed)
    changed = 0
    total = 0
    with open(target, "w") as out:
        for number, seq in enumerate(read_sequences(source)):
            bases = list(seq.upper())
            for i, base in enumerate(bases):
                if base in OTHERS:
                    total += 1
                    if rng.random() < rate:
                        bases[i] = rng.choice(OTHERS[base])
                        changed += 1
            out.write(f">mutant_{number} rate={rate} seed={seed}\n")
            text = "".join(bases)
            for start in range(0, len(text), 80):
                out.write(text[start : start + 80] + "\n")
    print(f"{changed}\t{changed / total:.6f}")


if __name__ == "__main__":
    main()
