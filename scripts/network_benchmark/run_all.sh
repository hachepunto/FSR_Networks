#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Orquestador del módulo network_benchmark: corre los pasos 01–09, 12/12b y 13/13b
# en orden, con el paralelismo indicado. Los pasos de inferencia (01, 02, 03,
# 06, 12, 13) saltan lo que ya existe; las evaluaciones (04, 04b, 07, 12b, 13b) y el
# resumen (08, 09) se recalculan siempre (son baratos). Se puede relanzar
# tras una interrupción.
# El paso 13 necesita las matrices <CANCER>_pooled: correr antes
# preprocess_expression.R con network_benchmark.pooled.write_matrices: true y
# pooled.expr_root apuntando a su work_root (otra carpeta, para no tocar la
# del pipeline).
# Los pasos 10–11 (apéndice técnico) se corren aparte, porque escriben en
# una carpeta de manuscrito que no es parte del repo:
#   Rscript scripts/network_benchmark/10_appendix_statistics.R <appendix_dir>
#   Rscript scripts/network_benchmark/11_appendix_numbers_tex.R <appendix_dir>
#
# Uso: bash scripts/network_benchmark/run_all.sh [corto_jobs] [aracne3_jobs] [aracne3_threads]
#   corto_jobs       corridas de corto en paralelo (cada una usa
#                    network_benchmark.threads hilos; default 4)
#   aracne3_jobs     corridas de ARACNe3 en paralelo (default 2)
#   aracne3_threads  hilos por corrida de ARACNe3 (default 18)
#   Tope de la máquina: corto_jobs*threads + aracne3_jobs*aracne3_threads
#   no debe pasar de (cores/2)-1.
#
# Variables de entorno respetadas: MRA_CONFIG (config alternativo),
# ARACNE3_DIR (build de ARACNe3 distinto del de config.yaml).
#
# Requisitos previos (no los hace este script):
#   - entorno config/env_network_benchmark.yaml activo, con
#     00_install_r_packages.R ya corrido;
#   - preprocess_expression.R y run_deseq2.R corridos para los 34
#     analysis_id (firmas en <work_root>/deseq2/<analysis_id>/stat_vector.rds);
#   - ARACNe3 compilado (scripts/setup_aracne3.sh).
# ============================================================

CORTO_JOBS="${1:-4}"
ARACNE3_JOBS="${2:-2}"
ARACNE3_THREADS="${3:-18}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
export LC_ALL=C.UTF-8 OPENBLAS_NUM_THREADS=1
unset OMP_NUM_THREADS   # ARACNe3 lo respetaría y correría con un solo hilo

CONFIG="${MRA_CONFIG:-config/config.yaml}"
cfg() {
    Rscript -e 'a <- commandArgs(TRUE); x <- yaml::yaml.load(readLines(a[1], warn = FALSE, encoding = "UTF-8")); for (k in strsplit(a[2], ".", fixed = TRUE)[[1]]) x <- x[[k]]; cat(if (is.null(x)) "" else x)' \
        "$CONFIG" "$1"
}
WORK_ROOT="$(cfg paths.work_root)"
OUT_DIR="$(cfg network_benchmark.out_dir)"; [[ -n "$OUT_DIR" ]] || OUT_DIR="$WORK_ROOT/network_benchmark"
LOGS="$OUT_DIR/logs"; mkdir -p "$LOGS"

# analysis_id y método desde config.yaml, igual que los scripts R
IDS="$(Rscript -e 'source("scripts/network_benchmark/lib_benchmark.R"); at <- analysis_table(bench_config()); cat(paste(at$analysis_id, at$network_method), sep = "\n")')"

for id in $(echo "$IDS" | awk '{print $1}'); do
    [[ -s "$WORK_ROOT/deseq2/$id/stat_vector.rds" ]] || { echo "ERROR: falta la firma DESeq2 de $id (scripts/run_deseq2.R $id)"; exit 1; }
done

step() { echo ">> $(date '+%F %T') $*"; }

step "01 CollecTRI"
Rscript scripts/network_benchmark/01_download_collectri.R > "$LOGS/01_collectri.log" 2>&1

step "02 corto en datos reales ($CORTO_JOBS en paralelo)"
echo "$IDS" | awk '{print $1}' | xargs -P "$CORTO_JOBS" -I{} sh -c \
    'Rscript scripts/network_benchmark/02_run_corto_real.R {} > "'"$LOGS"'/02_corto_{}.log" 2>&1 || { echo "ERROR en 02 {}"; exit 255; }'

step "03 ARACNe3 en datos reales ($ARACNE3_JOBS x $ARACNE3_THREADS hilos)"
echo "$IDS" | awk '$2=="ARACNe3"{print $1}' | xargs -P "$ARACNE3_JOBS" -I{} sh -c \
    'bash scripts/network_benchmark/03_run_aracne3_real.sh {} '"$ARACNE3_THREADS"' > "'"$LOGS"'/03_aracne3_{}.log" 2>&1 || { echo "ERROR en 03 {}"; exit 255; }'

step "04 evaluación real"
echo "$IDS" | awk '{print $1}' | xargs -P "$CORTO_JOBS" -I{} sh -c \
    'Rscript scripts/network_benchmark/04_evaluate_real.R {} > "'"$LOGS"'/04_eval_{}.log" 2>&1 || { echo "ERROR en 04 {}"; exit 255; }'
Rscript scripts/network_benchmark/04b_cross_cohort_real.R > "$LOGS/04b_cross_cohort.log" 2>&1

step "12 submuestreo de cohortes grandes ($CORTO_JOBS en paralelo)"
Rscript -e 'source("scripts/network_benchmark/lib_benchmark.R"); ss <- bench_config()$network_benchmark$subsample; for (id in unlist(ss$cohorts)) for (n in unlist(ss$n_grid)) for (r in seq_len(ss$replicates)) cat(id, n, r, "\n")' \
    | xargs -P "$CORTO_JOBS" -L 1 sh -c \
    'Rscript scripts/network_benchmark/12_subsample_real.R $0 $1 $2 > "'"$LOGS"'/12_sub_$0_n$1_rep$2.log" 2>&1 || { echo "ERROR en 12 $0 $1 $2"; exit 255; }'
Rscript scripts/network_benchmark/12b_evaluate_subsample.R > "$LOGS/12b_subsample.log" 2>&1

step "13 red del cáncer completo para cohortes chicas ($CORTO_JOBS en paralelo)"
Rscript -e 'source("scripts/network_benchmark/lib_benchmark.R"); cfg <- bench_config(); pp <- cfg$network_benchmark$pooled; for (c in unique(analysis_table(cfg)$cancer[analysis_table(cfg)$n_red < network_min_n(cfg)])) cat("full", c, "\n"); for (r in seq_len(pp$replicates)) for (id in unlist(pp$cohorts)) cat("surrogate", id, r, "\n")' \
    | xargs -P "$CORTO_JOBS" -L 1 sh -c \
    'Rscript scripts/network_benchmark/13_pooled_networks.R $0 $1 $2 > "'"$LOGS"'/13_$0_$1_$2.log" 2>&1 || { echo "ERROR en 13 $0 $1 $2"; exit 255; }'
Rscript scripts/network_benchmark/13b_evaluate_pooled.R > "$LOGS/13b_pooled.log" 2>&1

step "05 datos sintéticos"
Rscript scripts/network_benchmark/05_simulate_synthetic.R > "$LOGS/05_simulate.log" 2>&1

step "06-07 redes y evaluación sintéticas"
ls -d "$OUT_DIR"/synthetic/*_n*_rep* | xargs -P "$ARACNE3_JOBS" -I{} sh -c \
    'n=$(basename {}); bash scripts/network_benchmark/06_run_networks_synthetic.sh {} '"$ARACNE3_THREADS"' > "'"$LOGS"'/06_net_$n.log" 2>&1 && Rscript scripts/network_benchmark/07_evaluate_synthetic.R {} > "'"$LOGS"'/07_eval_$n.log" 2>&1 || { echo "ERROR en 06/07 $n"; exit 255; }'

step "08 tablas y figuras"
Rscript scripts/network_benchmark/08_summarize.R > "$LOGS/08_summarize.log" 2>&1

step "09 texto del suplementario"
Rscript scripts/network_benchmark/09_render_text.R > "$LOGS/09_render.log" 2>&1

step "listo: results/final/{tables,figures}/network_benchmark_* y docs/supplementary_network_method.md"
