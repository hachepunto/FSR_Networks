# corto sobre un dataset sintético, con las mismas configuraciones que
# 02_run_corto_real.R (current: p fijo sobre log2; lenient_r<r>: rangos y
# |r| >= r para cada r de corto.r_min_grid, antes del recorte top-k).
#
# Uso: Rscript scripts/network_benchmark/06_run_corto_synthetic.R <dataset_dir>
# Salida: <dataset_dir>/corto_{current,lenient_r<r>}.rds

suppressMessages(library(corto))
source("scripts/network_benchmark/lib_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/06_run_corto_synthetic.R <dataset_dir>")
ds <- args[1]

config <- bench_config()
nb <- config$network_benchmark

expr <- readRDS(file.path(ds, "expr.rds"))
tfs <- intersect(readLines(file.path(ds, "tfs.txt")), rownames(expr))
n <- ncol(expr)

ranks <- rank_rows(expr)
runs <- list(current = list(mat = expr, p = nb$previous_rule$fixed_p))
for (r in unlist(nb$corto$r_min_grid)) runs[[lenient_name(r)]] <- list(mat = ranks, p = r2p(r, n))

for (name in names(runs)) {
  out_file <- file.path(ds, paste0("corto_", name, ".rds"))
  if (file.exists(out_file)) next
  reg <- run_corto_or_empty(runs[[name]]$mat, tfs, nb$corto$nbootstraps, runs[[name]]$p, nb$threads)
  saveRDS(reg, out_file)
  message(basename(ds), " ", name, ": ", length(reg), " reguladores")
}
