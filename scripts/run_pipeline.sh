#!/usr/bin/env bash
# Main example: simulate 40 cells from the 4 genomes (10 per species,
# coverage 2x-30x, 0.1% substitution errors), sketch them, compare them and
# cluster them at several thresholds. Raw outputs go to results/pipeline/.
set -euo pipefail

bin=target/release/kmer-sketch
work=data/pipeline
out=results/pipeline
mkdir -p "$work" "$out"
cargo build --release --quiet

timings="$out/timings.tsv"
echo -e "step\twall_s" > "$timings"
timed() {
    local step="$1"
    shift
    local start end
    start="$(date +%s.%N)"
    "$@"
    end="$(date +%s.%N)"
    awk -v s="$step" -v a="$start" -v b="$end" 'BEGIN { printf "%s\t%.2f\n", s, b - a }' >> "$timings"
}

rm -rf "$work/cells" "$work/sketches"
timed simulate "$bin" simulate --cells 40 --min-coverage 2 --max-coverage 30 \
    --error-rate 0.001 --read-length 150 --seed 42 -o "$work/cells" data/genomes/*.fna
cp "$work/cells/cells.tsv" "$out/cells.tsv"
du -cb "$work"/cells/*.fastq | tail -n 1 | awk '{printf "total_fastq_bytes\t%d\n", $1}' > "$out/sizes.tsv"

timed sketch "$bin" sketch -k 21 --scaled 1000 -o "$work/sketches" "$work"/cells/*.fastq
for metric in jaccard containment distance; do
    timed "compare_$metric" "$bin" compare --metric "$metric" -o "$out/$metric.tsv" "$work"/sketches/*.sketch
done

summary="$out/cluster_summary.tsv"
echo -e "metric\tthreshold\tclusters\tARI\tpurity\tspecies_recovered" > "$summary"
cluster() {
    local metric="$1" threshold="$2"
    shift 2
    local report="$out/cluster_${metric}_${threshold}.txt"
    "$bin" cluster --threshold "$threshold" "$@" \
        --assignments "$out/assign_${metric}_${threshold}.tsv" "$out/$metric.tsv" > "$report"
    head -n 1 "$report" | awk -v m="$metric" -v t="$threshold" '{
        for (i = 2; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
        print m "\t" t "\t" v["clusters"] "\t" v["ARI"] "\t" v["purity"] "\t" v["species_recovered"] }' >> "$summary"
}
for t in 0.001 0.01 0.05 0.1 0.2 0.3 0.5 0.7; do cluster jaccard "$t"; done
for t in 0.01 0.05 0.1 0.3 0.5 0.7 0.9; do cluster containment "$t"; done
for t in 0.01 0.05 0.1 0.2; do cluster distance "$t" --distance; done

{
    echo -e "metric\tcross_max\twithin_nn_min\twithin_mean"
    for metric in jaccard containment; do
        echo -e "$metric\t$(python3 scripts/matrix_stats.py "$out/$metric.tsv")"
    done
} > "$out/matrix_stats.tsv"

# Keep only the assignment and report files of the default setting.
find "$out" -name 'assign_*' ! -name 'assign_jaccard_0.1.tsv' -delete
find "$out" -name 'cluster_*.txt' ! -name 'cluster_jaccard_0.1.txt' \
    ! -name 'cluster_containment_0.5.txt' -delete

cat "$summary" "$out/matrix_stats.tsv" "$timings"
