#!/usr/bin/env bash
# Robustness sweep: where does species recovery break down?
#
# For each setting (coverage range, error rate, genome set) we simulate 40
# cells, sketch them at scaled=1000 and scaled=100, compare with Jaccard and
# containment, and cluster at many thresholds. One row per
# (setting, scaled, metric, threshold) goes to results/sweep/sweep.tsv, and
# results/sweep/separation.tsv records how far apart the species are.
#
# Genome sets:
#   four     the 4 reference genomes (very distant species)
#   strains  the 4 genomes plus an E. coli copy with 1% substitutions
#            (data/genomes_strains/ecolimut.fna, a stand-in for a close strain)
set -euo pipefail

bin=target/release/kmer-sketch
work=data/sweep
out=results/sweep
mkdir -p "$work" "$out"
cargo build --release --quiet

# Genome set with an artificial close relative of E. coli.
strains=data/genomes_strains
if [[ ! -s "$strains/ecolimut.fna" ]]; then
    mkdir -p "$strains"
    cp data/genomes/*.fna "$strains/"
    python3 scripts/make_mutant.py data/genomes/ecoli.fna "$strains/ecolimut.fna" 0.01 1 \
        > "$out/ecolimut_substitutions.tsv"
fi

# name  min_cov  max_cov  error_rate  genome_dir
settings=(
    "cov2-30_err0.001     2     30    0.001 data/genomes"
    "cov2-30_err0.01      2     30    0.01  data/genomes"
    "cov2-30_err0.05      2     30    0.05  data/genomes"
    "cov0.5-2_err0.001    0.5   2     0.001 data/genomes"
    "cov0.1-0.5_err0.001  0.1   0.5   0.001 data/genomes"
    "cov0.02-0.1_err0.001 0.02  0.1   0.001 data/genomes"
    "cov0.005-0.02_err0.001 0.005 0.02 0.001 data/genomes"
    "strains_cov2-30_err0.001 2  30    0.001 $strains"
    "strains_cov0.1-0.5_err0.001 0.1 0.5 0.001 $strains"
)
jaccard_thresholds="0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 0.7"
containment_thresholds="0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 0.7 0.9"

table="$out/sweep.tsv"
separation="$out/separation.tsv"
echo -e "setting\tscaled\tmetric\tthreshold\tclusters\tARI\tpurity\tspecies_recovered" > "$table"
echo -e "setting\tscaled\tmetric\tcross_max\twithin_nn_min\twithin_mean\tmedian_sketch_size" > "$separation"

for line in "${settings[@]}"; do
    read -r name min_cov max_cov error genome_dir <<< "$line"
    echo "== $name"
    cells="$work/$name/cells"
    rm -rf "$work/$name"
    "$bin" simulate --cells 40 --min-coverage "$min_cov" --max-coverage "$max_cov" \
        --error-rate "$error" --seed 42 -o "$cells" "$genome_dir"/*.fna
    cp "$cells/cells.tsv" "$out/${name}_cells.tsv"

    for scaled in 1000 100; do
        sketches="$work/$name/sketches_$scaled"
        "$bin" sketch --scaled "$scaled" -o "$sketches" "$cells"/*.fastq 2> /dev/null
        median_size="$(for f in "$sketches"/*.sketch; do grep '^count' "$f" | cut -f2; done | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}')"
        for metric in jaccard containment; do
            matrix="$work/$name/${metric}_$scaled.tsv"
            "$bin" compare --metric "$metric" -o "$matrix" "$sketches"/*.sketch
            echo -e "$name\t$scaled\t$metric\t$(python3 scripts/matrix_stats.py "$matrix")\t$median_size" >> "$separation"
            if [[ "$metric" == jaccard ]]; then thresholds="$jaccard_thresholds"; else thresholds="$containment_thresholds"; fi
            for t in $thresholds; do
                "$bin" cluster --threshold "$t" "$matrix" | head -n 1 | awk \
                    -v s="$name" -v sc="$scaled" -v m="$metric" -v t="$t" '{
                    for (i = 2; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
                    print s "\t" sc "\t" m "\t" t "\t" v["clusters"] "\t" v["ARI"] "\t" v["purity"] "\t" v["species_recovered"] }' >> "$table"
            done
        done
    done
    rm -rf "$cells" # the FASTQ files are large; keep only sketches and matrices
done

# Best result per (setting, scaled, metric), and the threshold range that reaches it.
python3 - "$table" > "$out/sweep_best.tsv" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
groups = {}
for r in rows:
    groups.setdefault((r["setting"], r["scaled"], r["metric"]), []).append(r)
print("setting\tscaled\tmetric\tbest_ARI\tthresholds_at_best\tclusters_at_best\trecovered_at_best")
for (setting, scaled, metric), rs in groups.items():
    best = max(float(r["ARI"]) for r in rs)
    at_best = [r for r in rs if float(r["ARI"]) == best]
    ts = ",".join(r["threshold"] for r in at_best)
    print(f"{setting}\t{scaled}\t{metric}\t{best:.4f}\t{ts}\t{at_best[0]['clusters']}\t{at_best[0]['species_recovered']}")
EOF
cat "$out/sweep_best.tsv" "$separation"
