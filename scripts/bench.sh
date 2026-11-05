#!/usr/bin/env bash
# Benchmark `kmer-sketch count` against the pure-Python counter, and check
# that both produce byte-identical output.
#
# Inputs: the E. coli genome, and a simulated 50x read set of it.
# Each setting is run 3 times; the summary reports the median wall time.
# Writes results/bench_raw.tsv, results/bench_summary.tsv,
# results/count_check.tsv and results/machine.txt.
set -euo pipefail

bin=target/release/kmer-sketch
genome=data/genomes/ecoli.fna
work=data/bench
results=results
k=21
mkdir -p "$work" "$results"

cargo build --release --quiet

# Machine description.
{
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    echo "physical_cores: $(grep -m1 'cpu cores' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //')"
    echo "logical_cpus: $(nproc)"
    echo "mem_total: $(grep MemTotal /proc/meminfo | awk '{printf "%.1f GiB", $2/1048576}')"
    echo "kernel: $(uname -r)"
    echo "rustc: $(rustc --version)"
    echo "python: $(python3 --version)"
} > "$results/machine.txt"

# 50x read set: one simulated cell with coverage exactly 50.
"$bin" simulate --cells 1 --min-coverage 50 --max-coverage 50 --error-rate 0.001 \
    --seed 7 -o "$work/reads50x" "$genome"
reads="$work/reads50x/ecoli__cell000.fastq"

raw="$results/bench_raw.tsv"
echo -e "input\ttool\tthreads\trun\twall_s\tpeak_rss_mb" > "$raw"

run() {
    local input_name="$1" tool="$2" threads="$3" rep="$4"
    shift 4
    local out="$work/${input_name}.${tool}.t${threads}.tsv"
    local timing
    timing="$(python3 scripts/time_cmd.py "$out" "$@")"
    echo -e "${input_name}\t${tool}\t${threads}\t${rep}\t${timing}" | tee -a "$raw"
}

for input_name in genome reads50x; do
    if [[ "$input_name" == genome ]]; then input="$genome"; else input="$reads"; fi
    for rep in 1 2 3; do
        for threads in 1 2 4 8; do
            run "$input_name" rust "$threads" "$rep" "$bin" -t "$threads" count -k "$k" "$input"
        done
        run "$input_name" python 1 "$rep" python3 scripts/py_count.py -k "$k" "$input"
    done
done

# Median of the 3 runs per setting.
python3 - "$raw" > "$results/bench_summary.tsv" <<'EOF'
import csv, statistics, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
groups = {}
for r in rows:
    groups.setdefault((r["input"], r["tool"], r["threads"]), []).append(r)
print("input\ttool\tthreads\tmedian_wall_s\tmin_wall_s\tmax_wall_s\tmedian_peak_rss_mb")
for (inp, tool, threads), rs in groups.items():
    walls = [float(r["wall_s"]) for r in rs]
    rss = [float(r["peak_rss_mb"]) for r in rs]
    print(f"{inp}\t{tool}\t{threads}\t{statistics.median(walls):.2f}\t{min(walls):.2f}\t{max(walls):.2f}\t{statistics.median(rss):.0f}")
EOF

# Correctness: every Rust output must equal the Python output byte for byte.
check="$results/count_check.tsv"
echo -e "input\toutput\tlines\tmd5\tidentical_to_python" > "$check"
for input_name in genome reads50x; do
    python_md5="$(md5sum < "$work/${input_name}.python.t1.tsv" | cut -d' ' -f1)"
    for out in "$work/${input_name}".*.tsv; do
        md5="$(md5sum < "$out" | cut -d' ' -f1)"
        lines="$(wc -l < "$out")"
        same=no
        [[ "$md5" == "$python_md5" ]] && same=yes
        echo -e "${input_name}\t$(basename "$out")\t${lines}\t${md5}\t${same}" >> "$check"
    done
done
du -b "$reads" | awk '{printf "reads50x_fastq_bytes\t%d\n", $1}' >> "$results/machine.txt"
cat "$results/bench_summary.tsv" "$check"
