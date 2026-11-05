#!/usr/bin/env bash
# Check sketches and Jaccard estimates against exact k-mer sets computed in
# pure Python (scripts/py_sketch_check.py), and check Mash distances against
# known substitution rates. Needs the pipeline cells from run_pipeline.sh.
# Raw outputs go to results/sketch_check/.
set -euo pipefail

bin=target/release/kmer-sketch
work=data/sketch_check
out=results/sketch_check
mkdir -p "$work" "$out"
cargo build --release --quiet

# A. Six simulated cells (four E. coli, two S. aureus) at the default scaled=1000.
cells=(ecoli__cell001 ecoli__cell005 ecoli__cell009 ecoli__cell013 saureus__cell003 saureus__cell007)
cell_files=()
for c in "${cells[@]}"; do cell_files+=("data/pipeline/cells/$c.fastq"); done
grep -E "^($(IFS='|'; echo "${cells[*]}"))\." data/pipeline/cells/cells.tsv > "$out/cells_used.tsv"
"$bin" sketch --scaled 1000 -o "$work/cells_1000" "${cell_files[@]}"
"$bin" compare -o "$out/cells_rust_jaccard.tsv" "$work"/cells_1000/*.sketch
python3 scripts/py_sketch_check.py -k 21 --scaled 1000 --rust-sketch-dir "$work/cells_1000" \
    --rust-matrix "$out/cells_rust_jaccard.tsv" \
    --pairs-out "$out/cells_pairs.tsv" --hashes-out "$out/cells_hashes.tsv" "${cell_files[@]}"

# B. E. coli against copies with known substitution rates, plus B. subtilis.
genomes="$work/genomes"
mkdir -p "$genomes"
cp data/genomes/ecoli.fna data/genomes/bsubtilis.fna "$genomes/"
echo -e "genome\tsubstitutions\trealised_rate" > "$out/mutants.tsv"
for rate in 0.001 0.01 0.05; do
    name="ecolimut$rate"
    echo -e "$name\t$(python3 scripts/make_mutant.py data/genomes/ecoli.fna "$genomes/$name.fna" "$rate" 1)" >> "$out/mutants.tsv"
done
genome_files=("$genomes"/*.fna)
for scaled in 100 1000; do
    "$bin" sketch --scaled "$scaled" -o "$work/genomes_$scaled" "${genome_files[@]}"
    "$bin" compare -o "$out/genomes_rust_jaccard_$scaled.tsv" "$work/genomes_$scaled"/*.sketch
    "$bin" compare --metric distance -o "$out/genomes_rust_distance_$scaled.tsv" "$work/genomes_$scaled"/*.sketch
    python3 scripts/py_sketch_check.py -k 21 --scaled "$scaled" --rust-sketch-dir "$work/genomes_$scaled" \
        --rust-matrix "$out/genomes_rust_jaccard_$scaled.tsv" \
        --pairs-out "$out/genomes_pairs_$scaled.tsv" --hashes-out "$out/genomes_hashes_$scaled.tsv" \
        "${genome_files[@]}"
done

# Mash distance from E. coli to each mutant vs the realised substitution rate.
python3 - "$out" > "$out/mash_vs_true.tsv" <<'EOF'
import csv, sys
out = sys.argv[1]
rates = {r["genome"]: float(r["realised_rate"]) for r in csv.DictReader(open(f"{out}/mutants.tsv"), delimiter="\t")}
print("genome\ttrue_rate\tmash_distance_scaled100\tmash_distance_scaled1000")
def row(path):
    lines = [l.rstrip("\n").split("\t") for l in open(path)]
    names = lines[0][1:]
    ecoli = next(l for l in lines[1:] if l[0] == "ecoli")
    return dict(zip(names, map(float, ecoli[1:])))
d100 = row(f"{out}/genomes_rust_distance_100.tsv")
d1000 = row(f"{out}/genomes_rust_distance_1000.tsv")
for name, rate in rates.items():
    print(f"{name}\t{rate:.6f}\t{d100[name]:.6f}\t{d1000[name]:.6f}")
EOF
cat "$out"/cells_hashes.tsv "$out"/cells_pairs.tsv "$out"/genomes_pairs_100.tsv "$out/mash_vs_true.tsv"
