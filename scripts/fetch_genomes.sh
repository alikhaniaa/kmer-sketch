#!/usr/bin/env bash
# Download the four bacterial reference genomes used in the examples from
# NCBI nuccore with the E-utilities efetch endpoint (FASTA format).
#
# Usage: scripts/fetch_genomes.sh [output_dir]    (default: data/genomes)
#
# The genomes are not part of the repository (data/ is in .gitignore).
set -euo pipefail

out_dir="${1:-data/genomes}"
mkdir -p "$out_dir"
efetch_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

fetch() {
    local accession="$1"
    local name="$2"
    local path="$out_dir/$name.fna"
    if [[ -s "$path" ]]; then
        echo "$path already exists, skipping" >&2
        return
    fi
    echo "fetching $accession -> $path" >&2
    curl --fail --silent --show-error --retry 3 \
        "$efetch_url?db=nuccore&id=$accession&rettype=fasta&retmode=text" \
        -o "$path.part"
    if [[ "$(head -c 1 "$path.part")" != ">" ]]; then
        echo "error: $accession did not return FASTA" >&2
        exit 1
    fi
    mv "$path.part" "$path"
    sleep 1 # stay well below the NCBI limit of 3 requests per second
}

fetch NC_000913.3 ecoli        # Escherichia coli K-12 MG1655
fetch NC_000964.3 bsubtilis    # Bacillus subtilis 168
fetch NC_007795.1 saureus      # Staphylococcus aureus NCTC 8325
fetch NC_002516.2 paeruginosa  # Pseudomonas aeruginosa PAO1

# Print name, length and FASTA header of each genome (stdout, TSV).
for file in "$out_dir"/*.fna; do
    header="$(head -n 1 "$file")"
    length="$(grep -v '^>' "$file" | tr -d '\n\r' | wc -c)"
    printf '%s\t%s\t%s\n' "$(basename "$file" .fna)" "$length" "${header#>}"
done
