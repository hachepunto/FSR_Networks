# Compara FSR con DPI por blanco (target, el actual) y por triángulos (14).
#
# Sintéticos (de evaluation.rds de 07, que evalúa A_tri_k<main> si existe):
#   aristas (precisión, recall, FP indirectos por 100) y fracción de MR
#   verdaderos en el top-10, por régimen y n; diferencia pareada por dataset.
# Reales (dpi_triangle.cohorts): NaRnEA con la firma de DESeq2; top-20 de
#   target vs triangle y de cada uno vs ARACNe3-top100; tamaño de las redes.
#
# Criterio de adopción (fijado antes de correr): triangle reemplaza a target
# si mejora la fracción de MR verdaderos en el top-10 en sintéticos por más
# de 0.05 (≈ 1 EE) sin bajar la precisión de aristas.
#
# Uso: Rscript scripts/network_benchmark/14b_evaluate_dpi_triangle.R
# Salida: <out_dir>/real/evaluation/dpi_triangle.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
k <- nb$top_k_main
A <- paste0("A_k", k); TRI <- paste0("A_tri_k", k)

# ── Sintéticos ──
ev <- lapply(list.files(dirs$synthetic, "^evaluation\\.rds$", recursive = TRUE, full.names = TRUE), readRDS)
ev <- ev[vapply(ev, function(e) TRI %in% e$edges$regulon, logical(1))]
stopifnot("Ningún sintético evaluado con triangle (correr 14 y luego 07)" = length(ev) > 0)
ed <- do.call(rbind, lapply(ev, function(e) e$edges[e$edges$regulon %in% c(A, TRI), ]))
mr <- do.call(rbind, lapply(ev, function(e) e$master_regulators[e$master_regulators$regulon %in% c(A, TRI) &
                                                                  e$master_regulators$min_targets == nb$narnea$min_targets, ]))
pair <- function(d, v) {
  w <- reshape(d[, c("dataset", "regime", "n", "regulon", v)], idvar = c("dataset", "regime", "n"),
               timevar = "regulon", direction = "wide")
  w$diff <- w[[paste0(v, ".", TRI)]] - w[[paste0(v, ".", A)]]
  w
}
mr_w <- pair(mr, "precision_top10"); prec_w <- pair(ed, "precision"); rec_w <- pair(ed, "recall")
ind_w <- pair(ed, "fp_indirect_per_100_edges"); ne_w <- pair(ed, "n_edges")
stopifnot(identical(mr_w$dataset, prec_w$dataset), identical(mr_w$dataset, ind_w$dataset))
diffs <- data.frame(regime = mr_w$regime, n = mr_w$n, mr_diff = mr_w$diff, precision_diff = prec_w$diff,
                    recall_diff = rec_w$diff, indirect100_diff = ind_w$diff, n_edges_diff = ne_w$diff)
syn <- aggregate(. ~ regime + n, diffs, mean)
syn_level <- aggregate(cbind(precision_top10) ~ regulon + regime + n, mr, mean)

# ── Reales ──
at <- analysis_table(config)
real <- list()
for (id in unlist(nb$dpi_triangle$cohorts)) {
  f <- file.path(dirs$out, "real", "corto_triangle", paste0(id, "_", lenient_name(nb$corto$r_min), "_triangle.rds"))
  if (!file.exists(f)) { message(id, ": falta la red triangle; se omite"); next }
  tri_full <- readRDS(f)
  tri <- topk_regulon(tri_full, k)
  tgt <- topk_regulon(readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds"))), k)
  ges <- read_deseq2_stat(config, id)
  nar <- function(reg) run_narnea(to_narnea(reg, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed)
  r_tgt <- nar(tgt); r_tri <- nar(tri)
  r_B <- if (at$network_method[at$analysis_id == id] == "ARACNe3")
    nar(topk_regulon(aracne3_to_regulon(read_aracne3_net(dirs, id), config$network_inference$aracne3$subnets), k)) else NULL
  n_edges <- function(reg) sum(vapply(reg, function(x) length(x$tfmode), integer(1)))
  real[[length(real) + 1]] <- data.frame(
    analysis_id = id, n = at$n_red[at$analysis_id == id],
    edges_full_target = n_edges(readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds")))),
    edges_full_triangle = n_edges(tri_full), minutes_triangle = attr(tri_full, "params")$minutes,
    top20_target_vs_triangle = mr_concordance(r_tgt, r_tri)["top20"], rho_target_vs_triangle = mr_concordance(r_tgt, r_tri)["rho"],
    top20_target_vs_B = mr_concordance(r_tgt, r_B)["top20"], top20_triangle_vs_B = mr_concordance(r_tri, r_B)["top20"])
}
real <- do.call(rbind, real); rownames(real) <- NULL

saveRDS(list(synthetic_paired = list(mr = mr_w, precision = prec_w, recall = rec_w, indirect = ind_w),
             synthetic_summary = syn, synthetic_levels = syn_level, real = real),
        file.path(dirs$real_eval, "dpi_triangle.rds"))

message("Sintéticos, triangle - target (media sobre réplicas):"); print(syn, digits = 3)
message("Fracción de MR verdaderos en el top-10:"); print(reshape(syn_level, idvar = c("regime", "n"), timevar = "regulon", direction = "wide"), digits = 3)
message("Diferencia global en MR top-10: media ", round(mean(mr_w$diff), 3), "; precisión de aristas: media ", round(mean(prec_w$diff), 3))
message("\nReales:"); print(real, digits = 3)
