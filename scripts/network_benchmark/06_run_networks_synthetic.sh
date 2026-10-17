#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Redes sobre un dataset sintético: corto (06_run_corto_synthetic.R) y
# ARACNe3 (vía scripts/run_aracne3.sh, con los parámetros de
# network_inference.aracne3 de config.yaml).
#
# Uso: bash scripts/network_benchmark/06_run_networks_synthetic.sh <dataset_dir> [threads]
#   ARACNE3_DIR / MRA_CONFIG: igual que en 03_run_aracne3_real.sh.
# Salida: <dataset_dir>/corto_{current,lenient_r<r>}.rds
#         <dataset_dir>/aracne3/consolidated-net_<dataset>.tsv
# ============================================================

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Uso: bash scripts/network_benchmark/06_run_networks_synthetic.sh <dataset_dir> [threads]"
    exit 1
fi

DS="$(realpath "$1")"
NAME="$(basename "$DS")"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG="${MRA_CONFIG:-$REPO_ROOT/config/config.yaml}"

cfg() {
    Rscript -e 'a <- commandArgs(TRUE); x <- yaml::yaml.load(readLines(a[1], warn = FALSE, encoding = "UTF-8")); for (k in strsplit(a[2], ".", fixed = TRUE)[[1]]) x <- x[[k]]; cat(x)' \
        "$CONFIG" "$1"
}

ARACNE3_DIR="${ARACNE3_DIR:-$(cfg network_inference.aracne3.dir)}"
THREADS="${2:-18}"

cd "$REPO_ROOT"
Rscript scripts/network_benchmark/06_run_corto_synthetic.R "$DS"

if [[ ! -s "$DS/aracne3/consolidated-net_${NAME}.tsv" ]]; then
    bash scripts/run_aracne3.sh "$ARACNE3_DIR" "$DS/matrix.txt" "$DS/tfs.txt" "$DS/aracne3" "$NAME" "$THREADS" \
        "$(cfg network_inference.aracne3.subnets)" "$(cfg network_inference.aracne3.seed)" \
        "$(cfg network_inference.aracne3.alpha)"
fi
