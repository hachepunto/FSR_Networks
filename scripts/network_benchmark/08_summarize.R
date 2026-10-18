# Reúne las evaluaciones reales (04) y sintéticas (07) en las tablas y
# figuras del material suplementario. Etiquetas en inglés (manuscrito).
#
# Uso: Rscript scripts/network_benchmark/08_summarize.R
# Salidas (results/final/):
#   tables/network_benchmark_real_structure.tsv
#   tables/network_benchmark_real_concordance.tsv
#   tables/network_benchmark_real_nonmonotonic.tsv
#   tables/network_benchmark_runtime.tsv
#   tables/network_benchmark_synthetic_edges.tsv
#   tables/network_benchmark_synthetic_master_regulators.tsv
#   tables/network_benchmark_real_cross_cohort{,_decomposition}.tsv
#   tables/network_benchmark_real_subsample{,_summary}.tsv
#   tables/network_benchmark_synthetic_{average_precision,nonmonotonic,truth_regulon_nonmonotone}.tsv
#   tables/network_benchmark_key_numbers.tsv   (cada cifra citada en el texto)
#   tables/network_benchmark_collectri_{edges,provenance}.tsv (copia archivada de la referencia usada)
#   tables/network_benchmark_real_cross_cohort.tsv
#   tables/network_benchmark_synthetic_nonmonotonic.tsv
#   figures/network_benchmark_fig{1..7}_*.pdf
#
# Aborta si faltan evaluaciones (reales o sintéticas): el resumen debe ser
# sobre el conjunto completo, no sobre lo que haya alcanzado a correr.

suppressMessages({ library(ggplot2); library(patchwork) })
source("scripts/network_benchmark/lib_benchmark.R")

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
k_main <- nb$top_k_main
A <- paste0("A_k", k_main); B <- paste0("B_k", k_main); C <- paste0("aracne3_c", nb$aracne3_count_filter)
at <- analysis_table(config)

# Copia archivada de CollecTRI tal como se usó (la API sirve la versión vigente).
ct_src <- file.path(dirs$reference, "collectri.tsv")
if (file.exists(ct_src)) {
  file.copy(ct_src, file.path(dirs$tables, "network_benchmark_collectri_edges.tsv"), overwrite = TRUE)
  file.copy(file.path(dirs$reference, "collectri_provenance.tsv"), file.path(dirs$tables, "network_benchmark_collectri_provenance.tsv"), overwrite = TRUE)
}

write_tsv <- function(df, name) {
  write.table(df, file.path(dirs$tables, paste0("network_benchmark_", name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
save_fig <- function(p, name, w, h) ggsave(file.path(dirs$figures, paste0("network_benchmark_", name, ".pdf")), p, width = w, height = h)
theme_set(theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom"))

key <- list()
# Enteros (conteos) exactos; el resto a 3 cifras significativas.
add_key <- function(id, value, description) {
  v <- if (!is.na(value) && value == round(value)) value else signif(value, 3)
  key[[length(key) + 1]] <<- data.frame(id = id, value = v, description = description)
}

r_grid <- unlist(nb$corto$r_min_grid)
Ar <- sprintf("A_r%.1f", r_grid)
regulon_labels <- c(current = "corto, fixed p = 1e-8 (log2TPM)",
                    setNames(sprintf("corto, ranks, |r| >= %.1f, top-%d", r_grid, k_main), Ar),
                    setNames(sprintf("Scheme A: corto, ranks, |r| >= %.1f, top-%d", nb$corto$r_min, unlist(nb$top_k)), paste0("A_k", unlist(nb$top_k))),
                    aracne3_full = "ARACNe3, full consolidated",
                    setNames(paste0("ARACNe3, support >= ", nb$aracne3_count_filter), C),
                    setNames(paste0("Scheme B: ARACNe3, top-", k_main), B),
                    B_matched = "ARACNe3 by support, same edge count as scheme A",
                    truth = "True network (upper bound)")

# ============================================================
# DATOS REALES
# ============================================================

ev_files <- file.path(dirs$real_eval, paste0(at$analysis_id, ".rds"))
stopifnot("Faltan evaluaciones reales (04_evaluate_real.R)" = all(file.exists(ev_files)))
ev <- lapply(ev_files, readRDS)
message(length(ev), " analysis_id reales evaluados")
cross_file <- file.path(dirs$real_eval, "cross_cohort.rds")
stopifnot("Falta cross_cohort.rds (04b_cross_cohort_real.R)" = file.exists(cross_file))
cross <- readRDS(cross_file)

if (length(ev)) {
  regulon_structure_tbl <- do.call(rbind, lapply(ev, function(x) cbind(x$structure, n = x$n, network_method = x$network_method)))
  concordance <- do.call(rbind, lapply(ev, function(x) cbind(x$concordance, n = x$n)))
  nonmono <- do.call(rbind, lapply(ev, function(x) x$nonmonotonic))
  write_tsv(regulon_structure_tbl, "real_structure"); write_tsv(concordance, "real_concordance")
  if (!is.null(nonmono)) write_tsv(nonmono, "real_nonmonotonic")

  # ── Figura 1: estructura (current vs A) ──
  s1 <- regulon_structure_tbl[regulon_structure_tbl$regulon %in% c("current", A), ]
  s1$regulon <- factor(regulon_labels[s1$regulon], levels = regulon_labels[c("current", A)])
  f1a <- ggplot(s1, aes(n, regulators_ge_min, colour = regulon)) + geom_point() + scale_x_log10() +
    labs(x = "Tumour samples (n, log scale)", y = paste0("Regulators with >= ", nb$narnea$min_targets, " targets"), colour = NULL)
  f1b <- ggplot(s1, aes(n, median_targets, colour = regulon)) + geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x = "Tumour samples (n, log scale)", y = "Median targets per regulator (log scale)", colour = NULL)
  f1c <- ggplot(s1, aes(n, frac_negative, colour = regulon)) + geom_point() + scale_x_log10() +
    labs(x = "Tumour samples (n, log scale)", y = "Fraction of negative-mode edges", colour = NULL)
  s1c <- regulon_structure_tbl[regulon_structure_tbl$regulon %in% c("current", A, "aracne3_full", C, B), ]
  s1c$regulon <- factor(regulon_labels[s1c$regulon], levels = regulon_labels[c("current", A, "aracne3_full", C, B)])
  f1d <- ggplot(s1c, aes(regulon, collectri_enrichment)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) + geom_hline(yintercept = 1, linetype = 2) +
    labs(x = NULL, y = "CollecTRI enrichment over random") + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig((f1a | f1b) / (f1c | f1d) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A"),
           "fig1_regulon_structure", 10, 8)

  cur <- regulon_structure_tbl[regulon_structure_tbl$regulon == "current", ]; a_s <- regulon_structure_tbl[regulon_structure_tbl$regulon == A, ]
  min_n_prev <- nb$previous_rule$aracne3_min_n
  small <- cur$n < min_n_prev
  thr <- do.call(rbind, lapply(ev, function(x) data.frame(n = x$n, r = x$corto_params$current$r_threshold)))
  add_key("real_n_analysis", nrow(cur), "real analysis_ids evaluated")
  add_key("current_small_n_median_regulators_ge_min", median(cur$regulators_ge_min[small]), "current corto, n<previous threshold: median regulators with >= min targets")
  add_key("A_small_n_median_regulators_ge_min", median(a_s$regulators_ge_min[a_s$n < 95]), "scheme A, n<previous threshold: median regulators with >= min targets")
  add_key("current_median_targets_range_low", min(cur$median_targets), "current corto: min over analysis_ids of median targets per regulator")
  add_key("current_median_targets_range_high", max(cur$median_targets), "current corto: max over analysis_ids of median targets per regulator")
  add_key("A_median_targets_range_low", min(a_s$median_targets), "scheme A: min over analysis_ids of median targets per regulator")
  add_key("A_median_targets_range_high", max(a_s$median_targets), "scheme A: max over analysis_ids of median targets per regulator")
  add_key("current_cor_n_vs_median_targets", cor(cur$n, cur$median_targets, method = "spearman"), "current corto: Spearman(n, median targets)")
  add_key("A_cor_n_vs_median_targets", cor(a_s$n, a_s$median_targets, method = "spearman"), "scheme A: Spearman(n, median targets)")
  add_key("n_min", min(thr$n), "smallest cohort size")
  add_key("n_max", max(thr$n), "largest cohort size")
  add_key("current_r_threshold_at_n_min", thr$r[which.min(thr$n)], "current corto: effective |r| threshold at the smallest n")
  add_key("current_r_threshold_at_n_max", thr$r[which.max(thr$n)], "current corto: effective |r| threshold at the largest n")
  add_key("current_small_n_median_frac_negative", median(cur$frac_negative[small]), "current corto, n<previous threshold: median fraction of negative-mode edges")
  add_key("current_small_n_count_frac_negative_lt_0.01", sum(cur$frac_negative[small] < 0.01), "current corto, n<previous threshold: cohorts with < 1% negative-mode edges")
  add_key("n_small_cohorts", sum(small), "cohorts below the previous ARACNe3 threshold")
  add_key("A_n_cohorts_lt30", sum(a_s$n < 30), "cohorts with n < 30")
  for (r in c("current", A, "aracne3_full", C, B)) {
    v <- regulon_structure_tbl$collectri_enrichment[regulon_structure_tbl$regulon == r]
    if (length(v)) add_key(paste0("collectri_enrichment_median_", r), median(v), paste("median CollecTRI enrichment, all cohorts where defined:", r))
  }
  # Sobre las MISMAS cohortes (las de ARACNe3) y al mismo k, que es lo comparable:
  big <- unique(regulon_structure_tbl$analysis_id[regulon_structure_tbl$regulon == B])
  for (r in c("current", A, "aracne3_full", C, B)) {
    v <- regulon_structure_tbl$collectri_enrichment[regulon_structure_tbl$regulon == r & regulon_structure_tbl$analysis_id %in% big]
    add_key(paste0("collectri_enrichment_median_big_", r), median(v), paste("median CollecTRI enrichment over the ARACNe3 cohorts:", r))
  }
  af <- regulon_structure_tbl[regulon_structure_tbl$regulon == "aracne3_full" & regulon_structure_tbl$n <= 130, ]
  add_key("aracne3_full_n_le130_edges_min", min(af$n_edges), "full ARACNe3 network, cohorts with 95 <= n <= 130: min edges")
  add_key("aracne3_full_n_le130_edges_max", max(af$n_edges), "full ARACNe3 network, cohorts with 95 <= n <= 130: max edges")
  bf <- regulon_structure_tbl[regulon_structure_tbl$regulon == B & regulon_structure_tbl$n <= 130, ]
  add_key("B_n_le130_regulators_ge_min_min", min(bf$regulators_ge_min), "scheme B, cohorts with 95 <= n <= 130: min evaluable regulators")
  add_key("B_n_le130_regulators_ge_min_max", max(bf$regulators_ge_min), "scheme B, cohorts with 95 <= n <= 130: max evaluable regulators")
  smallc <- unique(regulon_structure_tbl$analysis_id[!regulon_structure_tbl$analysis_id %in% big])
  for (r in c("current", A)) {
    v <- regulon_structure_tbl$collectri_enrichment[regulon_structure_tbl$regulon == r & regulon_structure_tbl$analysis_id %in% smallc]
    add_key(paste0("collectri_enrichment_median_small_", r), median(v), paste("median CollecTRI enrichment over the cohorts below the previous threshold:", r))
  }
  add_key("A_collectri_enrichment_cor_n", cor(a_s$n, a_s$collectri_enrichment, method = "spearman"), "scheme A: Spearman(n, CollecTRI enrichment)")
  add_key("A_collectri_enrichment_min_n_lt30", min(a_s$collectri_enrichment[a_s$n < 30]), "scheme A, n<30: min CollecTRI enrichment")
  add_key("A_collectri_enrichment_max_n_lt30", max(a_s$collectri_enrichment[a_s$n < 30]), "scheme A, n<30: max CollecTRI enrichment")
  for (k in setdiff(unlist(nb$top_k), k_main)) add_key(paste0("collectri_enrichment_median_A_k", k), median(regulon_structure_tbl$collectri_enrichment[regulon_structure_tbl$regulon == paste0("A_k", k)]), paste0("median CollecTRI enrichment, scheme A with k = ", k))
  # "current rule" en su dominio real (n < 95)
  add_key("current_small_n_median_targets_low", min(cur$median_targets[small]), "corto fixed p, n<previous threshold: min median targets")
  add_key("current_small_n_median_targets_high", max(cur$median_targets[small]), "corto fixed p, n<previous threshold: max median targets")
  add_key("current_small_n_cor_n_vs_median_targets", cor(cur$n[small], cur$median_targets[small], method = "spearman"), "corto fixed p, n<previous threshold: Spearman(n, median targets)")
  thr_small <- thr[thr$n < min_n_prev, ]
  add_key("current_r_threshold_at_n95", max(thr_small$r[thr_small$n == max(thr_small$n)]), "corto fixed p: effective |r| threshold at the largest n < 95")
  add_key("n_max_small", max(thr_small$n), "largest cohort size below the previous ARACNe3 threshold")

  # ── Figura 2: A vs B frente a la variabilidad interna de ARACNe3 ──
  c2 <- concordance[(concordance$x == A & concordance$y == B) | (concordance$x == B & concordance$y == C), ]
  if (nrow(c2)) {
    c2$comparison <- ifelse(c2$x == A, "Scheme A vs scheme B (between methods)",
                            paste0("Scheme B vs ARACNe3 support >= ", nb$aracne3_count_filter, " (within ARACNe3)"))
    long <- do.call(rbind, lapply(c("rho", "top20", "top50"), function(m)
      data.frame(analysis_id = c2$analysis_id, n = c2$n, comparison = c2$comparison, metric = m, value = c2[[m]])))
    long$metric <- factor(long$metric, c("rho", "top20", "top50"),
                          c("Spearman rho of NES", "Top-20 overlap (|NES|)", "Top-50 overlap (|NES|)"))
    f2 <- ggplot(long, aes(n, value, colour = comparison)) + geom_line(aes(group = analysis_id), colour = "grey80") +
      geom_point(size = 2) + facet_wrap(~metric) + scale_x_log10() +
      labs(x = "Tumour samples (n, log scale)", y = "Agreement of master-regulator results", colour = NULL)
    save_fig(f2, "fig2_scheme_A_vs_B", 10, 4.5)
    ab <- c2[c2$x == A, ]; bc <- c2[c2$x == B, ]
    for (m in c("rho", "top20", "top50")) {
      add_key(paste0("AvsB_median_", m), median(ab[[m]], na.rm = TRUE), paste("scheme A vs B, median over ARACNe3 analysis_ids:", m))
      add_key(paste0("BvsC_median_", m), median(bc[[m]], na.rm = TRUE), paste("scheme B vs support-filtered ARACNe3, median:", m))
      add_key(paste0("AvsB_min_", m), min(ab[[m]], na.rm = TRUE), paste("scheme A vs B, minimum:", m))
    }
    add_key("AvsB_n_aracne3", nrow(ab), "ARACNe3 analysis_ids in the A vs B comparison")
    bf <- concordance[concordance$x == B & concordance$y == "aracne3_full", ]
    for (m in c("rho", "top20", "top50")) add_key(paste0("BvsFull_median_", m), median(bf[[m]], na.rm = TRUE), paste("scheme B vs full ARACNe3 network, median:", m))
    add_key("AvsB_rho_lower_than_BvsC_count", sum(ab$rho[match(bc$analysis_id, ab$analysis_id)] < bc$rho), "cohorts where rho(A,B) < rho(B, support-filtered ARACNe3)")
    d20 <- ab$top20[match(bc$analysis_id, ab$analysis_id)] - bc$top20
    add_key("AvsB_minus_BvsC_top20_min", min(d20), "paired difference top-20 overlap, (A vs B) - (B vs C): min")
    add_key("AvsB_minus_BvsC_top20_max", max(d20), "paired difference top-20 overlap, (A vs B) - (B vs C): max")
    add_key("AvsB_rho_at_smallest_n", ab$rho[which.min(ab$n)], "scheme A vs B: Spearman rho in the smallest ARACNe3 cohort")
    add_key("AvsB_smallest_n", min(ab$n), "smallest ARACNe3 cohort size")
  }

  # Sensibilidad (k, aw) y relación con el esquema actual
  sens <- concordance[concordance$x == A & grepl("^A_k|_awlik$|^current$", concordance$y), ]
  for (y in unique(sens$y)) add_key(paste0("A_vs_", y, "_median_rho"), median(sens$rho[sens$y == y], na.rm = TRUE),
                                    paste("scheme A vs", y, ": median Spearman rho of NES"))
  ct_c <- concordance[concordance$y == "CollecTRI", ]
  for (x in unique(ct_c$x)) add_key(paste0(x, "_vs_CollecTRI_median_rho"), median(ct_c$rho[ct_c$x == x], na.rm = TRUE),
                                    paste(x, "vs CollecTRI: median Spearman rho of NES"))

  # ── Línea nula y efecto entre subtipos (04b) ──
  write_tsv(rbind(cbind(comparison = "within_cancer", cross$within[, c("x", "y", "n_shared", "rho", "top20", "top50")]),
                  cbind(comparison = "null_other_cancer_network", setNames(cross$null[, c("analysis_id", "network_from", "n_shared", "rho", "top20", "top50")], c("x", "y", "n_shared", "rho", "top20", "top50")))),
            "real_cross_cohort")
  for (m in c("rho", "top20", "top50")) {
    add_key(paste0("within_cancer_median_", m), median(cross$within[[m]], na.rm = TRUE), paste("subtypes of the same cancer, scheme A, median:", m))
    add_key(paste0("null_median_", m), median(cross$null[[m]], na.rm = TRUE), paste("null (own signature, scheme A network from another cancer), median:", m))
    add_key(paste0("null_max_", m), max(cross$null[[m]], na.rm = TRUE), paste("null, max:", m))
  }
  add_key("within_cancer_n_pairs", nrow(cross$within), "pairs of subtypes of the same cancer")
  add_key("null_n_pairs", nrow(cross$null), "null pairs")
  add_key("method_effect_top20_turnover", 1 - median(concordance$top20[concordance$x == A & concordance$y == B], na.rm = TRUE), "top-20 turnover when changing method (1 - median overlap A vs B)")
  add_key("subtype_effect_top20_turnover", 1 - median(cross$within$top20, na.rm = TRUE), "top-20 turnover between subtypes of the same cancer (1 - median overlap)")
  if (!is.null(cross$decomp)) {
    for (cmp in c("within_net", "within_sig")) {
      d <- cross$decomp[cross$decomp$comparison == cmp, ]
      for (m in c("rho", "top20", "top50")) add_key(paste0(cmp, "_median_", m), median(d[[m]], na.rm = TRUE),
        paste(if (cmp == "within_net") "own signature, scheme-A network of another subtype of the same cancer" else "own network, signature of another subtype of the same cancer", ", median:", m))
      add_key(paste0(cmp, "_n_pairs"), nrow(d), paste(cmp, "directed pairs"))
    }
    add_key("network_effect_between_subtypes_top20_turnover", 1 - median(cross$decomp$top20[cross$decomp$comparison == "within_net"], na.rm = TRUE), "top-20 turnover when swapping only the network between subtypes of the same cancer")
    add_key("signature_effect_between_subtypes_top20_turnover", 1 - median(cross$decomp$top20[cross$decomp$comparison == "within_sig"], na.rm = TRUE), "top-20 turnover when swapping only the signature between subtypes of the same cancer")
    write_tsv(cross$decomp, "real_cross_cohort_decomposition")
    wt <- wilcox.test(cross$within$top20, cross$null$top20)$p.value
    add_key("within_vs_null_top20_wilcoxon_p", wt, "Wilcoxon p, top-20 overlap within-cancer vs null")
    add_key("within_cancer_top20_q25", quantile(cross$within$top20, 0.25, na.rm = TRUE), "within-cancer top-20 overlap, 25th percentile")
    add_key("within_cancer_top20_q75", quantile(cross$within$top20, 0.75, na.rm = TRUE), "within-cancer top-20 overlap, 75th percentile")
    add_key("null_top20_q25", quantile(cross$null$top20, 0.25, na.rm = TRUE), "null top-20 overlap, 25th percentile")
    add_key("null_top20_q75", quantile(cross$null$top20, 0.75, na.rm = TRUE), "null top-20 overlap, 75th percentile")
  }

  # ── Figura 3: no monotonía en ARACNe3 ──
  if (!is.null(nonmono)) {
    nm_real <- nonmono[!is.na(nonmono$min_count), ]
    for (m in sort(unique(nm_real$min_count))) {
      v <- nm_real$frac_abs_scc_lt_0.2[nm_real$min_count == m]
      add_key(paste0("nonmono_frac_low_scc_min_c", m), min(v), paste0("support >= ", m, ": min fraction of edges with |Spearman| < 0.2"))
      add_key(paste0("nonmono_frac_low_scc_max_c", m), max(v), paste0("support >= ", m, ": max fraction of edges with |Spearman| < 0.2"))
    }
    vB <- nonmono$frac_abs_scc_lt_0.2[nonmono$regulon == B]
    add_key("nonmono_frac_low_scc_median_B", median(vB), "scheme B regulons: median fraction of edges with |Spearman| < 0.2")
    add_key("nonmono_frac_low_scc_min_B", min(vB), "scheme B regulons: min fraction with |Spearman| < 0.2")
    add_key("nonmono_frac_low_scc_max_B", max(vB), "scheme B regulons: max fraction with |Spearman| < 0.2")
    nonmono <- nm_real
    for (m in sort(unique(nonmono$min_count))) {
      d <- nonmono[nonmono$min_count == m, ]
      add_key(paste0("nonmono_frac_low_scc_c", m), median(d$frac_abs_scc_lt_0.2), paste0("support >= ", m, ": median fraction of edges with |Spearman| < 0.2"))
      add_key(paste0("nonmono_hi_mi_low_scc_total_c", m), sum(d$n_top_quartile_mi_abs_scc_lt_0.2), paste0("support >= ", m, ": total top-quartile-MI edges with |Spearman| < 0.2"))
      add_key(paste0("nonmono_hi_mi_total_c", m), sum(d$n_top_quartile_mi), paste0("support >= ", m, ": total top-quartile-MI edges"))
      add_key(paste0("nonmono_cor_mi_scc_c", m), median(d$spearman_mi_vs_abs_scc), paste0("support >= ", m, ": median Spearman(MI, |Spearman|)"))
    }
  }

  # ── Figura 6A-B: sensibilidad al umbral |r| en datos reales ──
  s6 <- regulon_structure_tbl[regulon_structure_tbl$regulon %in% Ar, ]
  s6$threshold <- factor(sub("A_r", "|r| >= ", s6$regulon), levels = sub("A_r", "|r| >= ", Ar))
  f6a <- ggplot(s6, aes(n, median_targets, colour = threshold)) + geom_point() + scale_x_log10() +
    labs(x = "Tumour samples (n, log scale)", y = "Median targets per regulator", colour = NULL)
  c6 <- concordance[concordance$x %in% Ar & concordance$y == B, ]
  c6$threshold <- factor(sub("A_r", "|r| >= ", c6$x), levels = sub("A_r", "|r| >= ", Ar))
  f6b <- ggplot(c6, aes(n, top20, colour = threshold)) + geom_point(position = position_dodge(width = 0.03)) + scale_x_log10() +
    labs(x = "Tumour samples (n, log scale)", y = "Top-20 overlap with scheme B", colour = NULL)
  for (a in Ar) {
    sa <- s6[s6$regulon == a, ]
    add_key(paste0(a, "_median_targets_min"), min(sa$median_targets), paste(a, ": min over analysis_ids of median targets per regulator"))
    add_key(paste0(a, "_regulators_ge_min_min"), min(sa$regulators_ge_min), paste(a, ": min over analysis_ids of regulators with >= min targets"))
    add_key(paste0(a, "_n_ids_median_targets_eq_k"), sum(sa$median_targets >= k_main), paste(a, ": analysis_ids whose median regulon size reaches k"))
    ca <- concordance[concordance$x == a & concordance$y == B, ]
    if (nrow(ca)) for (m in c("rho", "top20", "top50"))
      add_key(paste0(a, "_vsB_median_", m), median(ca[[m]], na.rm = TRUE), paste(a, "vs scheme B, median:", m))
    add_key(paste0(a, "_collectri_enrichment_median"), median(regulon_structure_tbl$collectri_enrichment[regulon_structure_tbl$regulon == a]),
            paste(a, ": median CollecTRI enrichment"))
  }
  fig6_real <- list(f6a, f6b)

  # ── Costo computacional ──
  # ARACNe3: tiempo de pared y hilos de run_info.tsv (03_run_aracne3_real.sh);
  # si falta, suma de los tiempos que reporta el log de ARACNe3.
  runtime <- do.call(rbind, lapply(ev, function(x) {
    d <- file.path(dirs$aracne3, x$analysis_id)
    info <- file.path(d, "run_info.tsv")
    lg <- file.path(d, paste0("log_", x$analysis_id, ".txt"))
    a_min <- NA; a_thr <- NA
    if (file.exists(info)) {
      ri <- read.delim(info); a_min <- ri$wall_seconds / 60; a_thr <- ri$threads
    } else if (file.exists(lg)) {
      l <- grep("time: [0-9]+s", readLines(lg), value = TRUE)
      a_min <- sum(as.numeric(sub(".*time: ([0-9]+)s.*", "\\1", l))) / 60
    }
    data.frame(analysis_id = x$analysis_id, n = x$n,
               corto_lenient_minutes = x$corto_params[[lenient_name(nb$corto$r_min)]]$minutes, corto_threads = nb$threads,
               aracne3_minutes = a_min, aracne3_threads = a_thr)
  }))
  write_tsv(runtime, "runtime")
  rt <- runtime[!is.na(runtime$aracne3_minutes), ]
  if (nrow(rt)) {
    add_key("runtime_aracne3_total_hours", sum(rt$aracne3_minutes) / 60, "ARACNe3: total wall-clock hours")
    add_key("runtime_aracne3_max_hours", max(rt$aracne3_minutes) / 60, "ARACNe3: max wall-clock hours for one analysis_id")
    add_key("runtime_corto_same_ids_total_hours", sum(rt$corto_lenient_minutes) / 60, "corto (lenient): total wall-clock hours, same analysis_ids")
    add_key("runtime_aracne3_core_hours", sum(rt$aracne3_minutes * rt$aracne3_threads) / 60, "ARACNe3: core-hours (wall x threads), same analysis_ids")
    add_key("runtime_corto_core_hours", sum(rt$corto_lenient_minutes * rt$corto_threads) / 60, "corto (lenient): core-hours (wall x threads), same analysis_ids")
    add_key("runtime_aracne3_threads", unique(rt$aracne3_threads)[1], "ARACNe3 threads per run")
    add_key("runtime_corto_threads", nb$threads, "corto threads per run")
  }
}

# ============================================================
# SUBMUESTREO DE COHORTES REALES (12, 12b)
# ============================================================
sub_file <- file.path(dirs$real_eval, "subsample.rds")
stopifnot("Falta subsample.rds (12_subsample_real.R + 12b_evaluate_subsample.R)" = file.exists(sub_file))
{
  sub <- readRDS(sub_file)
  write_tsv(sub, "real_subsample")
  sub_sum <- do.call(rbind, lapply(split(sub, sub$n), function(d) data.frame(
    n = d$n[1], n_runs = nrow(d),
    top20_median = median(d$top20, na.rm = TRUE), top20_q25 = quantile(d$top20, 0.25, na.rm = TRUE), top20_q75 = quantile(d$top20, 0.75, na.rm = TRUE),
    top50_median = median(d$top50, na.rm = TRUE), rho_median = median(d$rho, na.rm = TRUE),
    edge_jaccard_per_tf_median = median(d$edge_jaccard_per_tf_median), n_evaluable_median = median(d$n_evaluable))))
  write_tsv(sub_sum, "real_subsample_summary")
  for (i in seq_len(nrow(sub_sum))) for (v in c("top20_median", "top50_median", "rho_median", "edge_jaccard_per_tf_median", "n_evaluable_median"))
    add_key(sprintf("subsample_n%d_%s", sub_sum$n[i], v), sub_sum[[v]][i], sprintf("subsampled real cohorts, n = %d: %s over cohorts x replicates", sub_sum$n[i], v))
  add_key("subsample_n_cohorts", length(unique(sub$analysis_id)), "cohorts subsampled")
  add_key("subsample_replicates", max(sub$replicate), "replicates per n")
  add_key("subsample_n_full_min", min(sub$n_full), "smallest full cohort subsampled")
  add_key("subsample_n_full_max", max(sub$n_full), "largest full cohort subsampled")
  # n mínimo cuyo top-20 mediano alcanza el yardstick A vs B
  ab_med <- median(concordance$top20[concordance$x == A & concordance$y == B], na.rm = TRUE)
  reach <- sub_sum$n[sub_sum$top20_median >= ab_med]
  add_key("subsample_min_n_reaching_AvsB_top20", if (length(reach)) min(reach) else NA, "smallest n whose median top-20 overlap with the full cohort reaches the A-vs-B median")
  sub$cohort_lab <- paste0(sub$analysis_id, " (n = ", sub$n_full, ")")
  yard <- data.frame(y = c(ab_med, median(cross$within$top20, na.rm = TRUE), median(cross$null$top20, na.rm = TRUE)),
                     lab = c("Scheme A vs B (same cohort)", "Subtypes of the same cancer", "Null (network from another cancer)"))
  f7a <- ggplot(sub, aes(n, top20, colour = cohort_lab)) + stat_summary(fun = median, geom = "line") + geom_point(size = 0.8, alpha = 0.6, position = position_jitter(width = 0.02, height = 0)) +
    geom_hline(data = yard, aes(yintercept = y, linetype = lab), colour = "grey40") + scale_x_log10(breaks = unlist(nb$subsample$n_grid)) +
    labs(x = "Subsampled tumours (n, log scale)", y = "Top-20 MR overlap with the full cohort", colour = NULL, linetype = NULL)
  f7b <- ggplot(sub, aes(n, edge_jaccard_per_tf_median, colour = cohort_lab)) + stat_summary(fun = median, geom = "line") + geom_point(size = 0.8, alpha = 0.6) +
    scale_x_log10(breaks = unlist(nb$subsample$n_grid)) + labs(x = "Subsampled tumours (n, log scale)", y = "Median per-regulator edge Jaccard vs full cohort", colour = NULL)
  f7c <- ggplot(sub, aes(n, n_evaluable, colour = cohort_lab)) + stat_summary(fun = median, geom = "line") + geom_point(size = 0.8, alpha = 0.6) +
    scale_x_log10(breaks = unlist(nb$subsample$n_grid)) + labs(x = "Subsampled tumours (n, log scale)", y = paste0("Regulators with >= ", nb$narnea$min_targets, " targets"), colour = NULL)
  f7a <- f7a + guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 3))
  save_fig((f7a / (f7b | f7c)) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A"), "fig7_subsampling_real", 11, 10)
}

# ============================================================
# DATOS SINTÉTICOS
# ============================================================

syn_files <- list.files(dirs$synthetic, pattern = "^evaluation\\.rds$", recursive = TRUE, full.names = TRUE)
n_syn_expected <- length(nb$synthetic$noise_sd) * length(unlist(nb$synthetic$n_network)) * nb$synthetic$replicates
stopifnot("Faltan evaluaciones sintéticas (07_evaluate_synthetic.R)" = length(syn_files) == n_syn_expected)
message(length(syn_files), " datasets sintéticos evaluados")
if (length(syn_files)) {
  syn <- lapply(syn_files, readRDS)
  se <- do.call(rbind, lapply(syn, `[[`, "edges"))
  sm <- do.call(rbind, lapply(syn, `[[`, "master_regulators"))
  sap <- do.call(rbind, lapply(syn, `[[`, "average_precision"))
  summarise <- function(d, vars, by) {
    g <- interaction(d[by], drop = TRUE, sep = "|")
    out <- do.call(rbind, lapply(split(d, g), function(x) {
      cbind(x[1, by, drop = FALSE], n_replicates = nrow(x),
            as.data.frame(as.list(unlist(lapply(vars, function(v)
              c(setNames(mean(x[[v]], na.rm = TRUE), paste0(v, "_mean")), setNames(sd(x[[v]], na.rm = TRUE), paste0(v, "_sd"))))))))
    }))
    rownames(out) <- NULL
    out[order(out[[by[1]]], out[[by[2]]], out$n), ]
  }
  edge_vars <- c("n_edges", "precision", "recall", "recall_linear", "recall_monotone_nonlinear",
                 "recall_non_monotone", "sign_accuracy_monotone", "n_fp", "n_fp_indirect",
                 "fp_indirect_per_100_edges", "frac_fp_indirect")
  se_sum <- summarise(se, edge_vars, c("regime", "regulon", "n"))
  sm$regime_mt <- paste0(sm$regime, "|", sm$min_targets)
  sm_sum <- summarise(sm, c("auroc", "auroc_all", "precision_top10", "mr_sign_correct", "mr_evaluated", "n_regulators"),
                      c("regime_mt", "regulon", "n"))
  sm_sum <- cbind(regime = sub("\\|.*", "", sm_sum$regime_mt), min_targets = as.integer(sub(".*\\|", "", sm_sum$regime_mt)),
                  sm_sum[, setdiff(colnames(sm_sum), "regime_mt")])
  sap_sum <- summarise(sap, "average_precision", c("regime", "method", "n"))
  write_tsv(se_sum, "synthetic_edges")
  write_tsv(sap_sum, "synthetic_average_precision")
  snm <- do.call(rbind, lapply(syn, `[[`, "nonmonotonic"))
  snm_sum <- summarise(snm, c("n_edges", "frac_abs_scc_lt_0.2", "n_true_non_monotone_recovered",
                              "frac_true_non_monotone_abs_scc_lt_0.2", "frac_edges_true_non_monotone"), c("regime", "regulon", "n"))
  write_tsv(snm_sum, "synthetic_nonmonotonic")
  # Figura 3: fracción |Spearman| < 0.2 en ARACNe3 real vs sintético (control positivo)
  nm_real_plot <- nonmono; nm_real_plot$source <- "Real tumours"
  nm_syn_plot <- snm[!is.na(snm$min_count) & snm$n %in% c(100, 200, max(unlist(nb$synthetic$n_network))), ]
  nm_syn_plot$source <- paste0("Synthetic, n = ", nm_syn_plot$n)
  nm_plot <- rbind(nm_real_plot[, c("source", "min_count", "frac_abs_scc_lt_0.2")], nm_syn_plot[, c("source", "min_count", "frac_abs_scc_lt_0.2")])
  nm_plot$support <- factor(paste0(">= ", nm_plot$min_count, " / 100"), levels = paste0(">= ", c(1, nb$aracne3_count_filter, 50), " / 100"))
  f3 <- ggplot(nm_plot, aes(support, frac_abs_scc_lt_0.2, fill = source)) + geom_boxplot(position = position_dodge(0.8), outlier.size = 0.6) +
    labs(x = "ARACNe3 edge support (subnetworks)", y = "Fraction of edges with |Spearman| < 0.2",
         fill = "Synthetic networks contain 20 % non-monotone true edges") + guides(fill = guide_legend(nrow = 2, title.position = "top"))
  save_fig(f3, "fig3_nonmonotonic_real_vs_synthetic", 8, 5)
  for (m in c(100, 200)) for (rg in unique(snm_sum$regime)) {
    d <- snm_sum[snm_sum$regime == rg & snm_sum$n == m & snm_sum$regulon == C, ]
    if (nrow(d)) add_key(sprintf("syn_nonmono_%s_n%d_c%d_frac_low_scc", rg, m, nb$aracne3_count_filter), d$frac_abs_scc_lt_0.2_mean, sprintf("synthetic %s, n=%d, support >= %d: mean fraction of edges with |Spearman| < 0.2", rg, m, nb$aracne3_count_filter))
  }
  for (rg in unique(snm_sum$regime)) {
    d <- snm_sum[snm_sum$regime == rg & snm_sum$n == max(unlist(nb$synthetic$n_network)), ]
    for (i in seq_len(nrow(d))) {
      tag <- sub("^aracne3_", "", d$regulon[i])   # "aracne3_c10" -> "c10"; "B_k100" queda igual
      add_key(sprintf("syn_nonmono_%s_%s_frac_low_scc", rg, tag), d$frac_abs_scc_lt_0.2_mean[i], sprintf("synthetic %s, n max, %s: mean fraction of edges with |Spearman| < 0.2", rg, d$regulon[i]))
      add_key(sprintf("syn_nonmono_%s_%s_frac_true_nm", rg, tag), d$frac_edges_true_non_monotone_mean[i], sprintf("synthetic %s, n max, %s: mean fraction of edges that are true non-monotone", rg, d$regulon[i]))
      add_key(sprintf("syn_nonmono_%s_%s_n_true_nm", rg, tag), d$n_true_non_monotone_recovered_mean[i], sprintf("synthetic %s, n max, %s: mean true non-monotone edges recovered", rg, d$regulon[i]))
      if (!is.na(d$frac_true_non_monotone_abs_scc_lt_0.2_mean[i])) add_key(sprintf("syn_nonmono_%s_%s_frac_true_nm_low_scc", rg, tag), d$frac_true_non_monotone_abs_scc_lt_0.2_mean[i], sprintf("synthetic %s, n max, %s: mean fraction of recovered true non-monotone edges with |Spearman| < 0.2", rg, d$regulon[i]))
    }
  }
  write_tsv(sm_sum, "synthetic_master_regulators")
  se_mt <- sm[sm$regulon == A & sm$min_targets == nb$narnea$min_targets, ]
  se_vals <- tapply(se_mt$precision_top10, list(se_mt$regime, se_mt$n), function(v) sd(v) / sqrt(length(v)))
  add_key("syn_prec10_se_min", min(se_vals, na.rm = TRUE), "synthetic: min over (regime, n) of the s.e. of the top-10 fraction, scheme A")
  add_key("syn_prec10_se_max", max(se_vals, na.rm = TRUE), "synthetic: max over (regime, n) of the s.e. of the top-10 fraction, scheme A")
  # diferencia máxima entre esquemas a n máximo, mt = min_targets
  nmax <- max(unlist(nb$synthetic$n_network))
  for (rg in unique(sm_sum$regime)) {
    d <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == nb$narnea$min_targets & sm_sum$n == nmax & sm_sum$regulon %in% c(A, B, "current"), ]
    add_key(sprintf("syn_%s_nmax_prec10_range_three_schemes", rg), max(d$precision_top10_mean) - min(d$precision_top10_mean), sprintf("synthetic %s, n = %d: range of top-10 fraction across A, B and fixed p", rg, nmax))
    tr <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == nb$narnea$min_targets & sm_sum$regulon == "truth", ]
    aa <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == nb$narnea$min_targets & sm_sum$regulon == A, ]
    gap <- max(tr$precision_top10_mean[tr$n >= 100] - aa$precision_top10_mean[match(tr$n[tr$n >= 100], aa$n)])
    add_key(sprintf("syn_%s_A_gap_to_truth_max_n_ge100", rg), gap, sprintf("synthetic %s: max over n >= 100 of (truth - scheme A) top-10 fraction", rg))
  }
  # prueba con el regulón verdadero: aristas no monótonas con modo 0 / sin ellas / con signo perfecto
  tn <- do.call(rbind, lapply(syn, `[[`, "nonmono_nes_truth"))
  if (!is.null(tn)) {
    tn_sum <- aggregate(cbind(rho_vs_truth_full, top20_vs_truth_full, precision_top10) ~ regime + n + variant, tn, mean)
    write_tsv(tn_sum, "synthetic_truth_regulon_nonmonotone")
    for (i in seq_len(nrow(tn_sum))) add_key(sprintf("syn_truth_%s_n%d_%s_prec10", tn_sum$regime[i], tn_sum$n[i], tn_sum$variant[i]), tn_sum$precision_top10[i],
      sprintf("synthetic %s, n=%d, true regulon variant %s: mean top-10 fraction", tn_sum$regime[i], tn_sum$n[i], tn_sum$variant[i]))
  }

  show <- c("current", A, "aracne3_full", B)
  se$regulon_lab <- factor(regulon_labels[se$regulon], levels = regulon_labels[c("current", A, "aracne3_full", C, B)])
  rl <- do.call(rbind, lapply(c("recall_linear", "recall_monotone_nonlinear", "recall_non_monotone", "precision"), function(v)
    data.frame(se[se$regulon %in% show, c("regime", "n", "regulon_lab")], metric = v, value = se[[v]][se$regulon %in% show])))
  rl$metric <- factor(rl$metric, c("precision", "recall_linear", "recall_monotone_nonlinear", "recall_non_monotone"),
                      c("Precision", "Recall: linear", "Recall: monotone non-linear", "Recall: non-monotone"))
  f4 <- ggplot(rl, aes(n, value, colour = regulon_lab)) + stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.2) + facet_grid(regime ~ metric) + scale_x_log10() +
    labs(x = "Samples used for network inference (n, log scale)", y = "Mean over replicates (+- s.e.)", colour = NULL)
  save_fig(f4, "fig4_synthetic_edges", 12, 6)

  sm$regulon_lab <- factor(regulon_labels[sm$regulon], levels = regulon_labels[c("current", A, "aracne3_full", C, B, "truth")])
  sel <- sm$regulon %in% c(show, "truth")
  ml <- do.call(rbind, lapply(c("mr_evaluated", "auroc_all", "precision_top10"), function(v)
    data.frame(sm[sel, c("regime", "min_targets", "n", "regulon_lab")], metric = v, value = sm[[v]][sel])))
  ml$metric <- factor(ml$metric, c("mr_evaluated", "auroc_all", "precision_top10"),
                      c("True MRs evaluable (of 10)", "AUROC, all TFs (non-evaluable = undetected)", "Fraction of true MRs in top-10 |NES|"))
  ml$panel <- paste0(ml$regime, " signal, min. targets = ", ml$min_targets)
  f5 <- ggplot(ml, aes(n, value, colour = regulon_lab)) + stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.2) + facet_grid(metric ~ panel, scales = "free_y") + scale_x_log10() +
    labs(x = "Samples used for network inference (n, log scale)", y = "Mean over replicates (+- s.e.)", colour = NULL)
  save_fig(f5, "fig5_synthetic_master_regulators", 13, 8)

  for (rg in unique(sm_sum$regime)) for (mt in unique(sm_sum$min_targets)) for (r in c("current", A, "aracne3_full", B, "truth")) {
    d <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == mt & sm_sum$regulon == r, ]
    for (i in seq_len(nrow(d))) {
      add_key(sprintf("syn_%s_mt%d_%s_n%d_auroc", rg, mt, r, d$n[i]), d$auroc_all_mean[i], sprintf("synthetic %s, min targets %d, %s, n=%d: mean AUROC over all TFs", rg, mt, r, d$n[i]))
      add_key(sprintf("syn_%s_mt%d_%s_n%d_prec10", rg, mt, r, d$n[i]), d$precision_top10_mean[i], sprintf("synthetic %s, min targets %d, %s, n=%d: mean precision in top-10", rg, mt, r, d$n[i]))
      add_key(sprintf("syn_%s_mt%d_%s_n%d_mreval", rg, mt, r, d$n[i]), d$mr_evaluated_mean[i], sprintf("synthetic %s, min targets %d, %s, n=%d: mean true MRs evaluable", rg, mt, r, d$n[i]))
    }
  }
  for (rg in unique(se_sum$regime)) for (r in c("current", A, "aracne3_full", B)) {
    e <- se_sum[se_sum$regime == rg & se_sum$regulon == r, ]
    for (i in seq_len(nrow(e))) for (v in c("precision", "recall_linear", "recall_monotone_nonlinear", "recall_non_monotone", "n_fp_indirect", "fp_indirect_per_100_edges", "frac_fp_indirect"))
      add_key(sprintf("syn_%s_%s_n%d_%s", rg, r, e$n[i], v), e[[paste0(v, "_mean")]][i], sprintf("synthetic %s, %s, n=%d: mean %s", rg, r, e$n[i], v))
  }
  add_key("syn_prec10_random_baseline", nb$synthetic$n_master_regulators / nb$synthetic$n_tfs, "fraction of true MRs in a random top-10 (n_MR / n_TFs)")
  for (rg in unique(se_sum$regime)) for (r in c("B_matched")) {
    e <- se_sum[se_sum$regime == rg & se_sum$regulon == r, ]
    for (i in seq_len(nrow(e))) for (v in c("precision", "recall_linear", "n_fp_indirect", "fp_indirect_per_100_edges", "n_edges"))
      add_key(sprintf("syn_%s_%s_n%d_%s", rg, r, e$n[i], v), e[[paste0(v, "_mean")]][i], sprintf("synthetic %s, %s, n=%d: mean %s", rg, r, e$n[i], v))
  }
  for (rg in unique(sm_sum$regime)) { d <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == nb$narnea$min_targets & sm_sum$regulon == "truth", ]
    for (i in seq_len(nrow(d))) add_key(sprintf("syn_truth_%s_n%d_truth_full_prec10", rg, d$n[i]), d$precision_top10_mean[i], sprintf("synthetic %s, n=%d, true regulon with non-monotone edges at mode 0: mean top-10 fraction", rg, d$n[i])) }
  for (rg in unique(se_sum$regime)) for (r in c(A, B, "aracne3_full")) { e <- se_sum[se_sum$regime == rg & se_sum$regulon == r, ]
    for (i in seq_len(nrow(e))) add_key(sprintf("syn_%s_%s_n%d_n_edges", rg, r, e$n[i]), e$n_edges_mean[i], sprintf("synthetic %s, %s, n=%d: mean edges", rg, r, e$n[i])) }
  m6 <- sm[sm$regulon %in% c(Ar, B, "current", "truth") & sm$min_targets == nb$narnea$min_targets, ]
  m6$scheme <- factor(regulon_labels[m6$regulon], levels = regulon_labels[c(Ar, B, "current", "truth")])
  f6c <- ggplot(m6, aes(n, precision_top10, colour = scheme)) + stat_summary(fun = mean, geom = "line") +
    stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.2) + facet_wrap(~regime) + scale_x_log10() +
    labs(x = "Samples used for network inference (n, log scale)", y = "Fraction of true MRs in top-10 |NES|", colour = NULL)
  if (exists("fig6_real")) {
    save_fig((fig6_real[[1]] | fig6_real[[2]]) / f6c + plot_annotation(tag_levels = "A"), "fig6_threshold_sensitivity", 11, 9)
  }
  for (rg in unique(sm_sum$regime)) for (r in Ar) {
    d <- sm_sum[sm_sum$regime == rg & sm_sum$min_targets == nb$narnea$min_targets & sm_sum$regulon == r, ]
    for (i in seq_len(nrow(d))) {
      add_key(sprintf("syn_%s_mt%d_%s_n%d_prec10", rg, nb$narnea$min_targets, r, d$n[i]), d$precision_top10_mean[i],
              sprintf("synthetic %s, min targets %d, %s, n=%d: mean fraction of true MRs in top-10", rg, nb$narnea$min_targets, r, d$n[i]))
      add_key(sprintf("syn_%s_mt%d_%s_n%d_auroc", rg, nb$narnea$min_targets, r, d$n[i]), d$auroc_all_mean[i],
              sprintf("synthetic %s, min targets %d, %s, n=%d: mean AUROC over all TFs", rg, nb$narnea$min_targets, r, d$n[i]))
    }
    e <- se_sum[se_sum$regime == rg & se_sum$regulon == r, ]
    for (i in seq_len(nrow(e))) for (v in c("precision", "recall_linear", "recall_non_monotone"))
      add_key(sprintf("syn_%s_%s_n%d_%s", rg, r, e$n[i], v), e[[paste0(v, "_mean")]][i], sprintf("synthetic %s, %s, n=%d: mean %s", rg, r, e$n[i], v))
  }
  for (i in seq_len(nrow(sap_sum))) add_key(sprintf("syn_%s_%s_n%d_AP", sap_sum$regime[i], sap_sum$method[i], sap_sum$n[i]),
                                            sap_sum$average_precision_mean[i], "average precision of the full edge ranking")
}

if (length(key)) {
  key_df <- do.call(rbind, key)
  stopifnot("ids duplicados en key_numbers" = !any(duplicated(key_df$id)))
  write_tsv(key_df, "key_numbers")
}
message("Tablas: ", dirs$tables, "\nFiguras: ", dirs$figures)
