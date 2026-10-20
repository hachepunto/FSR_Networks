# Evalúa los regulones submuestreados (12) contra la referencia con todas las
# muestras (el regulón A_r<r_min> top-k de 02/04): NaRnEA con la firma fija
# de DESeq2, concordancia de reguladores maestros (rho de NES, top-20/50) y
# estabilidad de aristas (Jaccard con la referencia, por regulador y global),
# más cuántos reguladores quedan evaluables.
#
# Uso: Rscript scripts/network_benchmark/12b_evaluate_subsample.R
# Salida: <out_dir>/real/evaluation/subsample.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
ss <- nb$subsample
sub_dir <- file.path(dirs$out, "real", "subsample")

edge_keys <- function(reg) unlist(lapply(names(reg), function(tf) paste(tf, names(reg[[tf]]$tfmode))), use.names = FALSE)

at <- analysis_table(config)
rows <- list()
for (id in unlist(ss$cohorts)) {
  ref_reg <- topk_regulon(readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds"))), nb$top_k_main)
  ges <- read_deseq2_stat(config, id)
  ref_res <- run_narnea(to_narnea(ref_reg, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed)
  ref_keys <- edge_keys(ref_reg)
  n_full <- ncol(read_tumor_expr(config, id, at$cancer[at$analysis_id == id]))
  for (n_sub in unlist(ss$n_grid)) for (r in seq_len(ss$replicates)) {
    f <- file.path(sub_dir, sprintf("%s_n%d_rep%d.rds", id, n_sub, r))
    if (!file.exists(f)) next
    reg <- readRDS(f)
    res <- run_narnea(to_narnea(reg, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed)
    cc <- mr_concordance(ref_res, res)
    keys <- edge_keys(reg)
    # Jaccard por regulador (sobre los comunes) y global
    per_tf <- vapply(intersect(names(reg), names(ref_reg)), function(tf) {
      a <- names(reg[[tf]]$tfmode); b <- names(ref_reg[[tf]]$tfmode)
      length(intersect(a, b)) / length(union(a, b)) }, numeric(1))
    rows[[length(rows) + 1]] <- data.frame(
      analysis_id = id, n_full = n_full, n = n_sub, replicate = r,
      n_regulators = length(reg), n_evaluable = if (is.null(res)) 0 else nrow(res),
      n_shared = cc["n_shared"], rho = cc["rho"], as.list(cc[grepl("^top", names(cc))]),
      edge_jaccard_global = length(intersect(keys, ref_keys)) / length(union(keys, ref_keys)),
      edge_jaccard_per_tf_median = median(per_tf),
      minutes = attr(reg, "params")$minutes)
  }
}
out <- do.call(rbind, rows); rownames(out) <- NULL
saveRDS(out, file.path(dirs$real_eval, "subsample.rds"))
agg <- aggregate(cbind(top20, top50, rho, edge_jaccard_per_tf_median, n_evaluable) ~ n, out, median)
message("Mediana sobre cohortes y réplicas:"); print(agg, digits = 3)
