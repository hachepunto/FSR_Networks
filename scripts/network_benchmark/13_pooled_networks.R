# Redes del cáncer completo ("pooled") para cohortes con n < min_n. Dos modos:
#
#   surrogate <analysis_id> <replicate>
#     Cohorte grande X reducida a pooled.small_n tumores al azar (X "chica").
#     Escribe dos regulones de FSR sobre esos mismos tumores:
#       _pooled: otros subtipos del cáncer + los small_n de X (matriz pooled)
#       _own:    solo los small_n de X (matriz de X)
#     13b los compara contra la red de X con todas sus muestras.
#
#   full <CANCER>
#     Red sobre la matriz pooled completa del cáncer: la que usaría el pipeline.
#
# La firma no se toca aquí (13b usa la de DESeq2 de X completa): se aísla el
# efecto de la RED. Idempotente.
#
# Uso: Rscript scripts/network_benchmark/13_pooled_networks.R surrogate <analysis_id> <replicate>
#      Rscript scripts/network_benchmark/13_pooled_networks.R full <CANCER>
# Salida: <out_dir>/real/pooled/

suppressMessages(library(corto))
source("scripts/network_benchmark/lib_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
if (!(mode %in% c("surrogate", "full")) || length(args) != c(surrogate = 3, full = 2)[[mode]])
  stop("Uso: 13_pooled_networks.R surrogate <analysis_id> <replicate> | full <CANCER>")

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
pp <- nb$pooled
pool_dir <- file.path(dirs$out, "real", "pooled"); dir.create(pool_dir, recursive = TRUE, showWarnings = FALSE)

# FSR sobre una matriz (genes x muestras); quita genes sin varianza en
# el subconjunto (con rangos constantes la correlación no está definida).
fsr_network <- function(mat, tfs) {
  mat <- mat[matrixStats::rowSds(mat) > 0, , drop = FALSE]
  t0 <- Sys.time()
  reg <- topk_regulon(run_corto_or_empty(rank_rows(mat), intersect(tfs, rownames(mat)), nb$corto$nbootstraps,
                                         r2p(nb$corto$r_min, ncol(mat)), nb$threads), nb$top_k_main)
  attr(reg, "params") <- list(n = ncol(mat), n_genes = nrow(mat),
                              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  reg
}

tfs <- read_tf_list(config)

if (mode == "full") {
  cancer <- args[2]
  out_file <- file.path(pool_dir, sprintf("%s_pooled_full.rds", cancer))
  if (file.exists(out_file)) { message(basename(out_file), ": ya existe, se omite"); quit(save = "no") }
  reg <- fsr_network(read_pooled_expr(config, cancer), tfs)
  saveRDS(reg, out_file)
  message(sprintf("%s pooled (n = %d): %d reguladores (%.1f min)", cancer, attr(reg, "params")$n,
                  length(reg), attr(reg, "params")$minutes))
  quit(save = "no")
}

analysis_id <- args[2]; rep_id <- as.integer(args[3])
stem <- file.path(pool_dir, sprintf("%s_rep%d", analysis_id, rep_id))
if (file.exists(paste0(stem, "_pooled.rds")) && file.exists(paste0(stem, "_own.rds"))) {
  message(basename(stem), ": ya existe, se omite"); quit(save = "no")
}

at <- analysis_table(config)
cancer <- at$cancer[at$analysis_id == analysis_id]
own <- read_tumor_expr(config, analysis_id, cancer)
pooled <- read_pooled_expr(config, cancer)
stopifnot("Las muestras de X deben estar en la matriz pooled" = all(colnames(own) %in% colnames(pooled)))
if (pp$small_n >= ncol(own)) stop(analysis_id, ": small_n no es menor que la cohorte")

# Semilla determinista por (cohorte, réplica)
set.seed(pp$seed + 1e6 * match(analysis_id, unlist(pp$cohorts)) + rep_id)
keep_x <- sort(sample(colnames(own), pp$small_n))
pool_cols <- c(setdiff(colnames(pooled), colnames(own)), keep_x)

reg_pooled <- fsr_network(pooled[, pool_cols, drop = FALSE], tfs)
reg_own <- fsr_network(own[, keep_x, drop = FALSE], tfs)
for (r in c("reg_pooled", "reg_own")) {
  reg <- get(r)
  attr(reg, "params") <- c(attr(reg, "params"), list(analysis_id = analysis_id, replicate = rep_id))
  saveRDS(reg, paste0(stem, sub("reg", "", r), ".rds"))
}
writeLines(keep_x, paste0(stem, "_samples.txt"))
message(sprintf("%s rep %d: pooled n = %d (%d reg., %.1f min); propia n = %d (%d reg., %.1f min)",
                analysis_id, rep_id, length(pool_cols), length(reg_pooled), attr(reg_pooled, "params")$minutes,
                pp$small_n, length(reg_own), attr(reg_own, "params")$minutes))
