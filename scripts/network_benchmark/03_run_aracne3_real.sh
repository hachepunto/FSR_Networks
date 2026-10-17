#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# ARACNe3 sobre la matriz tumoral de un analysis_id, para el módulo
# network_benchmark (referencia de comparación en los analysis_id con
# n_red >= 95). Delega en scripts/run_aracne3.sh con los parámetros de
# network_inference.aracne3 de config.yaml, así que hereda su
# verificación del commit compilado.
#
# Uso: bash scripts/network_benchmark/03_run_aracne3_real.sh <analysis_id> [threads]
#   threads: default 18 (no network_inference.aracne3.threads = 50: run_all.sh
#     corre dos a la vez y el tope de la máquina es (cores/2)-1)
#   ARACNE3_DIR (variable de entorno, opcional): build de ARACNe3 a usar
#     en vez de network_inference.aracne3.dir (debe estar en el mismo commit).
#   MRA_CONFIG (variable de entorno, opcional): config alternativo.
#
# Salida: <out_dir>/real/aracne3/<analysis_id>/consolidated-net_<analysis_id>.tsv
# ============================================================

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Uso: bash scripts/network_benchmark/03_run_aracne3_real.sh <analysis_id> [threads]"
    exit 1
fi

ANALYSIS_ID="$1"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG="${MRA_CONFIG:-$REPO_ROOT/config/config.yaml}"

cfg() {
    Rscript -e 'a <- commandArgs(TRUE); x <- yaml::yaml.load(readLines(a[1], warn = FALSE, encoding = "UTF-8")); for (k in strsplit(a[2], ".", fixed = TRUE)[[1]]) x <- x[[k]]; cat(if (is.null(x)) "" else x)' \
        "$CONFIG" "$1"
}

WORK_ROOT="$(cfg paths.work_root)"
OUT_DIR="$(cfg network_benchmark.out_dir)"; [[ -n "$OUT_DIR" ]] || OUT_DIR="$WORK_ROOT/network_benchmark"
OUT_DIR="$OUT_DIR/real/aracne3/$ANALYSIS_ID"
ARACNE3_DIR="${ARACNE3_DIR:-$(cfg network_inference.aracne3.dir)}"
THREADS="${2:-18}"

MATRIX="$WORK_ROOT/networks/input/${ANALYSIS_ID}_matrix.txt"
TFS="$WORK_ROOT/networks/input/$(cfg network_inference.tfs.file)"

if [[ -s "$OUT_DIR/consolidated-net_${ANALYSIS_ID}.tsv" ]]; then
    echo ">> $ANALYSIS_ID: red ya existe, se omite"
    exit 0
fi

START=$(date +%s)
bash "$REPO_ROOT/scripts/run_aracne3.sh" "$ARACNE3_DIR" "$MATRIX" "$TFS" "$OUT_DIR" "$ANALYSIS_ID" "$THREADS" \
    "$(cfg network_inference.aracne3.subnets)" "$(cfg network_inference.aracne3.seed)" \
    "$(cfg network_inference.aracne3.alpha)"

# Costo de la corrida, para la tabla de tiempos del suplementario.
printf "threads\twall_seconds\n%s\t%s\n" "$THREADS" "$(( $(date +%s) - START ))" > "$OUT_DIR/run_info.tsv"
