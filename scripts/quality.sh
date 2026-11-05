#!/usr/bin/env bash
# Run the same checks as CI (format, lint, tests) and save a short log to
# results/quality.txt.
set -euo pipefail
mkdir -p results
{
    rustc --version
    cargo fmt --all -- --check && echo "cargo fmt --check: ok"
    cargo clippy --locked --all-targets -- -D warnings 2>&1 | tail -n 1
    echo "cargo clippy --locked --all-targets -- -D warnings: ok"
    cargo test --locked 2>&1 | grep -E "Running|test result|^test "
} > results/quality.txt
cat results/quality.txt
