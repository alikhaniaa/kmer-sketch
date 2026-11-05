#!/usr/bin/env bash
# Fine threshold scan (step 0.01) for the two "strains" settings of sweep.sh,
# where the coarse threshold grid misses the narrow working window.
# Reuses the matrices that sweep.sh left in data/sweep/.
# Output: results/sweep/strains_fine.tsv
set -euo pipefail

bin=target/release/kmer-sketch
out=results/sweep/strains_fine.tsv
echo -e "setting\tscaled\tmetric\tthreshold\tclusters\tARI\tpurity\tspecies_recovered" > "$out"
for name in strains_cov2-30_err0.001 strains_cov0.1-0.5_err0.001; do
    for scaled in 1000 100; do
        for metric in jaccard containment; do
            matrix="data/sweep/$name/${metric}_$scaled.tsv"
            for t in $(seq 0.01 0.01 0.99); do
                "$bin" cluster --threshold "$t" "$matrix" | head -n 1 | awk \
                    -v s="$name" -v sc="$scaled" -v m="$metric" -v t="$t" '{
                    for (i = 2; i <= NF; i++) { split($i, kv, "="); v[kv[1]] = kv[2] }
                    print s "\t" sc "\t" m "\t" t "\t" v["clusters"] "\t" v["ARI"] "\t" v["purity"] "\t" v["species_recovered"] }' >> "$out"
            done
        done
    done
done
# Threshold window with ARI = 1 (all 5 labels recovered), per setting.
awk -F'\t' 'NR > 1 && $6 == "1.0000" { k = $1 "\t" $2 "\t" $3; if (!(k in lo)) lo[k] = $4; hi[k] = $4 }
    END { print "setting\tscaled\tmetric\tfirst_t_ARI1\tlast_t_ARI1"; for (k in lo) print k "\t" lo[k] "\t" hi[k] }' "$out" \
    | sort > results/sweep/strains_window.tsv
awk -F'\t' 'NR > 1 { k = $1 "\t" $2 "\t" $3; if ($6 + 0 > best[k] + 0) best[k] = $6 }
    END { for (k in best) print k "\tbest_ARI=" best[k] }' "$out" | sort
cat results/sweep/strains_window.tsv
