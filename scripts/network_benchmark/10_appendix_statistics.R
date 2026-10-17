# Cálculos del apéndice técnico (manuscript/appendix_network_method/, para un
# lector de física estadística). Son análisis propios que NO están en el
# suplementario: eficiencia de estimadores, validez empírica de la DPI sobre
# correlación, aristas indirectas frente a correlación parcial regularizada, y
# contribución de las aristas no monótonas al NES de NaRnEA. Reutilizan los
# datos sintéticos del módulo (05) y una red real de ARACNe3.
#
# Uso: Rscript scripts/network_benchmark/10_appendix_statistics.R <out_dir>
#   out_dir: carpeta donde dejar figures/ y appendix_numbers.tsv
#   (p.ej. la carpeta de Drive manuscript/appendix_network_method)
# Entradas: intermedios del módulo según config (bench_dirs), más
#   network_benchmark_key_numbers.tsv ya generado.

suppressMessages({ library(ggplot2); library(patchwork); library(glasso); library(NaRnEA) })
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/10_appendix_statistics.R <out_dir>")
out <- args[1]
fig_dir <- file.path(out, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

config <- bench_config()
dirs <- bench_dirs(config, create = FALSE)
nb <- config$network_benchmark
set.seed(20261017)

theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom"))
save_fig <- function(p, name, w, h) ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
nums <- list()
add_num <- function(id, value, description) nums[[length(nums) + 1]] <<- data.frame(id = id, value = signif(value, 4), description = description)

# ============================================================
# A1. Eficiencia de estimadores: Pearson, Spearman y tres estimadores de MI
#     (binning ingenuo de 8 cajas; Miller-Madow con B = sqrt(n/5) cajas;
#     kNN de Kraskov-Stögbauer-Grassberger con k = 3), bajo tres modelos de
#     ruido calibrados con la curtosis real de la expresión: gaussiano,
#     t con 5 g.l. (curtosis en exceso 6, por encima del p99 real) y
#     gaussiano con 3 % de outliers (el patrón dominante en RNA-seq).
# ============================================================
# Para cada n y modelo: R réplicas de pares (x, y) con rho = 0.4. Se mide
# (i) la varianza y el sesgo de cada estimador respecto a su propio valor
# poblacional (Spearman estima (6/pi) asin(rho/2), no rho), (ii) la potencia
# para detectar dependencia contra un nulo permutado (alpha = 0.05), que es la
# cantidad comparable entre medidas de escala distinta.
mi_bins <- function(x, y, b, miller_madow = FALSE) {
  tx <- cut(rank(x, ties.method = "random"), b, labels = FALSE); ty <- cut(rank(y, ties.method = "random"), b, labels = FALSE)
  p <- table(factor(tx, 1:b), factor(ty, 1:b)) / length(x); px <- rowSums(p); py <- colSums(p)
  nz <- p > 0; mi <- sum(p[nz] * log(p[nz] / outer(px, py)[nz]))
  if (miller_madow) mi <- mi - (sum(nz) - sum(px > 0) - sum(py > 0) + 1) / (2 * length(x))
  mi
}
mi_ksg <- function(x, y, k = 3) {
  n <- length(x); x <- rank(x, ties.method = "random"); y <- rank(y, ties.method = "random")
  dx <- abs(outer(x, x, "-")); dy <- abs(outer(y, y, "-")); d <- pmax(dx, dy); diag(d) <- Inf
  eps <- apply(d, 1, function(r) sort(r)[k])
  nx <- rowSums(dx < eps) - 1; ny <- rowSums(dy < eps) - 1
  digamma(k) + digamma(n) - mean(digamma(nx + 1) + digamma(ny + 1))
}
simulate_pair <- function(n, rho, noise) {
  z <- matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  switch(noise,
    gaussian = z,
    t5 = z / sqrt(rchisq(n, 5) / 5),
    outliers = { k <- sample(n, max(1, round(0.03 * n))); z[k, ] <- z[k, ] + matrix(rnorm(2 * length(k), 0, 6), ncol = 2); z })
}
est_all <- function(z) {
  x <- z[, 1]; y <- z[, 2]; n <- length(x)
  c(pearson = cor(x, y), spearman = cor(x, y, method = "spearman"),
    mi_bins8 = mi_bins(x, y, 8), mi_mm = mi_bins(x, y, max(2, round(sqrt(n / 5))), TRUE), mi_ksg = mi_ksg(x, y))
}
ns <- c(27, 50, 100, 200, 400); R <- 200; rho0 <- 0.4
rho_s0 <- (6 / pi) * asin(rho0 / 2)   # valor poblacional de Spearman bajo normalidad
eff <- do.call(rbind, lapply(c("gaussian", "t5", "outliers"), function(noise) do.call(rbind, lapply(ns, function(n) {
  est <- t(replicate(R, est_all(simulate_pair(n, rho0, noise))))
  null <- t(replicate(R, est_all(simulate_pair(n, 0, noise))))
  crit <- apply(null, 2, quantile, 0.95)
  target <- c(pearson = rho0, spearman = rho_s0, mi_bins8 = NA, mi_mm = NA, mi_ksg = NA)
  data.frame(noise = noise, n = n, estimator = colnames(est),
             sd = apply(est, 2, sd), bias = colMeans(est) - target[colnames(est)],
             power = sapply(colnames(est), function(e) mean(est[, e] > crit[e])))
}))))
eff$estimator <- factor(eff$estimator, c("pearson", "spearman", "mi_bins8", "mi_mm", "mi_ksg"),
                        c("Pearson", "Spearman", "MI, 8 bins (naive)", "MI, Miller-Madow, sqrt(n/5) bins", "MI, KSG k = 3"))
eff$noise <- factor(eff$noise, c("gaussian", "t5", "outliers"), c("Gaussian", "t, 5 d.f. (excess kurtosis 6)", "3 % outliers"))
fA1 <- ggplot(eff[eff$estimator %in% c("Pearson", "Spearman"), ], aes(n, sd, colour = estimator)) + geom_line() + geom_point() + facet_wrap(~noise) + scale_x_log10() +
  scale_colour_manual(values = c(Pearson = "#F8766D", Spearman = "#B79F00"), guide = "none") +
  labs(x = "n (log scale)", y = "s.d. of the estimate (true rho = 0.4)", colour = NULL)
fA1b <- ggplot(eff, aes(n, power, colour = estimator)) + geom_line() + geom_point() + facet_wrap(~noise) + scale_x_log10() +
  labs(x = "n (log scale)", y = expression(Power~to~detect~rho==0.4~"(" * alpha == 0.05 * ")"), colour = NULL)
# curtosis real por gen (COAD CMS2), para calibrar el modelo de ruido
at <- analysis_table(config); real_id <- "COAD_cms2"
expr <- read_tumor_expr(config, real_id, at$cancer[at$analysis_id == real_id])
kurt <- apply(expr, 1, function(v) { s <- sd(v); if (s == 0) NA else mean((v - mean(v))^4) / s^4 - 3 })
kq <- quantile(kurt, c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99), na.rm = TRUE)
for (i in seq_along(kq)) add_num(sprintf("A1_real_kurtosis_q%s", sub("%", "", names(kq)[i])), kq[i], "A1: per-gene excess kurtosis of log2TPM, COAD CMS2, quantile")
add_num("A1_real_kurtosis_frac_gt3", mean(kurt > 3, na.rm = TRUE), "A1: fraction of genes with excess kurtosis > 3")
fA1c <- ggplot(data.frame(k = kurt[!is.na(kurt) & kurt < 10]), aes(k)) + geom_histogram(bins = 80, fill = "grey60") +
  geom_vline(xintercept = c(0, 6), linetype = 2) + labs(x = "Per-gene excess kurtosis, log2(TPM+1), COAD CMS2 (truncated at 10)", y = "Genes")
save_fig((fA1 / fA1b / fA1c) + plot_layout(guides = "collect", heights = c(1, 1.2, 0.7)) + plot_annotation(tag_levels = "A"), "figA1_estimator_efficiency", 9, 10)
for (i in seq_len(nrow(eff))) {
  tag <- sprintf("A1_%s_n%d_%s", c("gaussian", "t5", "outliers")[as.integer(eff$noise[i])], eff$n[i], c("pearson", "spearman", "mibins", "mimm", "miksg")[as.integer(eff$estimator[i])])
  add_num(paste0(tag, "_sd"), eff$sd[i], "A1 s.d."); add_num(paste0(tag, "_power"), eff$power[i], "A1 power")
  if (!is.na(eff$bias[i])) add_num(paste0(tag, "_bias"), eff$bias[i], "A1 bias vs own population value")
}
vr <- eff$sd[eff$noise == "Gaussian" & eff$estimator == "Pearson"]^2 / eff$sd[eff$noise == "Gaussian" & eff$estimator == "Spearman"]^2
add_num("A1_gaussian_var_ratio_pearson_over_spearman_mean", mean(vr), "A1: Var(r_P)/Var(r_S) under Gaussianity at rho = 0.4, mean over n (not the Pitman ARE at rho = 0)")
add_num("A1_pitman_ARE_spearman_vs_pearson", (3 / pi)^2, "A1: Pitman asymptotic relative efficiency of Spearman vs Pearson for testing rho = 0 under normality, (3/pi)^2")

# ============================================================
# A2. ¿Se cumple la DPI sobre |rho| en cadenas X -> Y -> Z?
#     (a) Teoría: lineal-gaussiano, rho_XZ = rho_XY * rho_YZ.
#     (b) Simulación con f no lineales (monótonas y no monótonas).
#     (c) Datos reales: tripletes TF -> TF -> target en una red de ARACNe3.
# ============================================================
f_list <- list(linear = function(x) x, sigmoid = function(x) plogis(3 * x), square = function(x) x^2, abs = function(x) abs(x))
dpi_sim <- do.call(rbind, lapply(names(f_list), function(f1) do.call(rbind, lapply(names(f_list), function(f2) {
  viol <- replicate(200, {
    x <- rnorm(400); y <- f_list[[f1]](x) + rnorm(400, sd = 0.5); z <- f_list[[f2]](y) + rnorm(400, sd = 0.5)
    r <- abs(c(xy = cor(x, y, method = "spearman"), yz = cor(y, z, method = "spearman"), xz = cor(x, z, method = "spearman")))
    c(viol = r["xz"] > min(r["xy"], r["yz"]), margin = r["xz"] - min(r["xy"], r["yz"]))
  })
  data.frame(f_xy = f1, f_yz = f2, frac_violations = mean(viol["viol.xz", ]), mean_margin = mean(viol["margin.xz", ]))
}))))
for (i in seq_len(nrow(dpi_sim))) add_num(sprintf("A2_sim_%s_%s_viol", dpi_sim$f_xy[i], dpi_sim$f_yz[i]), dpi_sim$frac_violations[i], "A2: fraction of chains violating |rho_XZ| <= min(|rho_XY|,|rho_YZ|)")

# Real: tripletes en ARACNe3 (COAD_cms2), aristas con soporte >= 10. Hasta 5
# targets por par (a, b): los tripletes NO son independientes; se reporta
# también el número de pares.
net <- read_aracne3_net(dirs, real_id); net <- net[net$count.values >= nb$aracne3_count_filter, ]
tfs <- intersect(read_tf_list(config), rownames(expr))
tf_tf <- net[net$regulator.values %in% tfs & net$target.values %in% tfs, ]
set.seed(1); tf_tf <- tf_tf[sample(nrow(tf_tf), min(3000, nrow(tf_tf))), ]
rk <- rank_rows(expr)
trip <- do.call(rbind, lapply(seq_len(nrow(tf_tf)), function(i) {
  a <- tf_tf$regulator.values[i]; b <- tf_tf$target.values[i]
  tg <- net$target.values[net$regulator.values == b & !net$target.values %in% c(a, b)]
  if (!length(tg)) return(NULL)
  tg <- sample(tg, min(5, length(tg)))
  rab <- abs(cor(rk[a, ], rk[b, ]))
  do.call(rbind, lapply(tg, function(cc) {
    rbc <- abs(cor(rk[b, ], rk[cc, ])); rac <- abs(cor(rk[a, ], rk[cc, ]))
    data.frame(a = a, b = b, c = cc, r_ab = rab, r_bc = rbc, r_ac = rac, direct_ac = any(net$regulator.values == a & net$target.values == cc))
  }))
}))
trip$viol <- trip$r_ac > pmin(trip$r_ab, trip$r_bc)
trip$a_beats_b <- trip$r_ac > trip$r_bc   # lo que la regla winner-take-all de corto evalúa de verdad
add_num("A2_real_n_triplets", nrow(trip), "A2: TF->TF->target triplets sampled from COAD_cms2 ARACNe3 (support>=10)")
add_num("A2_real_n_tf_pairs", length(unique(paste(trip$a, trip$b))), "A2: distinct TF-TF pairs (a,b) behind the triplets")
add_num("A2_real_n_tfs_in_matrix", length(tfs), "A2: regulators present in the COAD_cms2 matrix")
add_num("A2_real_frac_a_beats_b_all", mean(trip$a_beats_b), "A2: fraction of triplets with |rho_ac| > |rho_bc| (a would win c over b under corto's rule), all")
add_num("A2_real_frac_a_beats_b_no_direct", mean(trip$a_beats_b[!trip$direct_ac]), "A2: same, triplets with no direct a->c edge in ARACNe3")
add_num("A2_real_frac_a_beats_b_direct", mean(trip$a_beats_b[trip$direct_ac]), "A2: same, triplets where ARACNe3 also has a->c")
add_num("A2_real_frac_violations_all", mean(trip$viol), "A2: fraction violating the correlation DPI, all triplets")
add_num("A2_real_frac_violations_no_direct", mean(trip$viol[!trip$direct_ac]), "A2: fraction violating, triplets with no direct a->c edge in ARACNe3")
add_num("A2_real_frac_violations_direct", mean(trip$viol[trip$direct_ac]), "A2: fraction violating, triplets where ARACNe3 also has a->c")
fA2 <- ggplot(trip, aes(pmin(r_ab, r_bc), r_ac, colour = direct_ac)) + geom_point(alpha = 0.3, size = 0.7) + geom_abline(linetype = 2) +
  coord_equal() + labs(x = expression(min*"(|"*rho[ab]*"|, |"*rho[bc]*"|)"), y = expression("|"*rho[ac]*"|"), colour = "ARACNe3 has a→c edge") +
  ggtitle("COAD CMS2: TF(a) -> TF(b) -> target(c) triplets, support >= 10")
save_fig(fA2, "figA2_dpi_real_triplets", 6, 6.5)

# ============================================================
# A3. Aristas indirectas: corto+DPI, ARACNe3, y correlación parcial (glasso)
#     en un dataset sintético (strong, n = 400 y n = 100), sobre un bloque
#     TF + sus targets verdaderos + targets de TFs del mismo módulo.
# ============================================================
pc_eval <- function(ds_name) {
  d <- file.path(dirs$synthetic, ds_name)
  tr <- readRDS(file.path(d, "truth.rds")); x <- readRDS(file.path(d, "expr.rds")); n <- ncol(x)
  te <- tr$edges; module <- tr$tf_module
  # submatriz: 40 TFs (2 módulos completos) + 600 targets regulados por ellos
  mods <- sort(unique(module))[1:2]; tfs_sub <- names(module)[module %in% mods]
  tg_sub <- unique(te$target[te$tf %in% tfs_sub]); tg_sub <- head(tg_sub, 600)
  genes <- c(tfs_sub, tg_sub); xs <- t(rank_rows(x[genes, ])); S <- cor(xs)
  truth_key <- paste(te$tf, te$target)
  cand <- expand.grid(tf = tfs_sub, target = tg_sub, stringsAsFactors = FALSE); cand$key <- paste(cand$tf, cand$target)
  cand$true <- cand$key %in% truth_key
  indirect_of <- function(tf, tg) { tr_ <- te$tf[te$target == tg]; length(tr_) > 0 && any(module[tr_] == module[tf]) && !(paste(tf, tg) %in% truth_key) }
  cand$indirect <- mapply(indirect_of, cand$tf, cand$target)
  # (i) correlación marginal |rho|; (ii) glasso en dos lambdas; (iii) corto A_k100; (iv) ARACNe3 B_k100
  cand$rho <- S[cbind(match(cand$tf, genes), match(cand$target, genes))]
  res <- list()
  for (lam in c(0.05, 0.1, 0.2)) {
    g <- glasso(S, rho = lam, penalize.diagonal = FALSE); P <- -cov2cor(g$wi)
    cand[[paste0("pc_", lam)]] <- abs(P[cbind(match(cand$tf, genes), match(cand$target, genes))])
  }
  regA <- topk_regulon(readRDS(file.path(d, paste0("corto_", lenient_name(nb$corto$r_min), ".rds"))), nb$top_k_main)
  eA <- regulon_edges(regA); cand$A <- cand$key %in% paste(eA$tf, eA$target)
  nt <- read.delim(file.path(d, "aracne3", paste0("consolidated-net_", ds_name, ".tsv")))
  regB_full <- aracne3_to_regulon(nt, config$network_inference$aracne3$subnets)
  eB <- regulon_edges(topk_regulon(regB_full, nb$top_k_main)); cand$B <- cand$key %in% paste(eB$tf, eB$target)
  # Para comparar a igual número de aristas seleccionadas: top-m por |rho|, por
  # |pc| y por soporte de ARACNe3 (B emparejada), con m = aristas de A en el bloque
  m <- sum(cand$A)
  eBf <- regulon_edges(regB_full); eBf <- eBf[eBf$tf %in% tfs_sub & eBf$target %in% tg_sub, ]
  eBf$tb <- unlist(lapply(regB_full[unique(eBf$tf)], tiebreak_of), use.names = FALSE)[match(paste(eBf$tf, eBf$target), unlist(lapply(unique(eBf$tf), function(t) paste(t, names(regB_full[[t]]$tfmode)))))]
  Bm <- eBf[order(eBf$likelihood, eBf$tb, decreasing = TRUE), ][seq_len(min(m, nrow(eBf))), ]
  cand$B_matched <- cand$key %in% paste(Bm$tf, Bm$target)
  summ <- function(sel, name) data.frame(dataset = ds_name, n = n, method = name, n_selected = sum(sel),
                                         precision = mean(cand$true[sel]), recall = sum(cand$true & sel) / sum(cand$true),
                                         indirect_per_100 = 100 * sum(cand$indirect & sel) / sum(sel))
  rbind(summ(cand$A, "corto, ranks, |r|>=0.2, top-100 (scheme A)"),
        summ(cand$B, "ARACNe3 top-100 (scheme B)"),
        summ(cand$B_matched, sprintf("ARACNe3 by support, top-%d (B matched)", m)),
        summ(rank(-cand$rho) <= m, sprintf("marginal |rho|, top-%d", m)),
        summ(rank(-cand$pc_0.05) <= m, sprintf("glasso partial corr. (lambda 0.05), top-%d", m)),
        summ(rank(-cand$pc_0.1) <= m, sprintf("glasso partial corr. (lambda 0.1), top-%d", m)),
        summ(rank(-cand$pc_0.2) <= m, sprintf("glasso partial corr. (lambda 0.2), top-%d", m)))
}
pc <- rbind(pc_eval("strong_n400_rep1"), pc_eval("strong_n100_rep1"))
write.table(pc, file.path(out, "appendix_A3_partial_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(pc))) {
  tag <- gsub("[^A-Za-z0-9]+", "_", pc$method[i])
  add_num(sprintf("A3_n%d_%s_precision", pc$n[i], tag), pc$precision[i], "A3 precision"); add_num(sprintf("A3_n%d_%s_recall", pc$n[i], tag), pc$recall[i], "A3 recall")
  add_num(sprintf("A3_n%d_%s_indirect100", pc$n[i], tag), pc$indirect_per_100[i], "A3 indirect FP per 100 edges")
}
pl <- pc; pl$method <- factor(pl$method, unique(pc$method))
fA3 <- ggplot(pl, aes(recall, precision, colour = method, shape = factor(n))) + geom_point(size = 3) +
  labs(x = "Recall (true edges in block)", y = "Precision", colour = NULL, shape = "n") + guides(colour = guide_legend(ncol = 2))
fA3b <- ggplot(pl, aes(method, indirect_per_100, fill = factor(n))) + geom_col(position = "dodge") + coord_flip() +
  labs(x = NULL, y = "Indirect false positives per 100 selected edges", fill = "n")
save_fig(fA3 / fA3b + plot_annotation(tag_levels = "A"), "figA3_indirect_edges_partial_correlation", 9, 9)

# ============================================================
# A4. Contribución de las aristas no monótonas al NES de NaRnEA, con la red
#     VERDADERA sintética (07_evaluate_synthetic.R, agregado sobre réplicas):
#     regulón verdadero completo (am = 0 en las no monótonas), sin ellas, y
#     con ellas pero con signo "perfecto" (|am| = 1), que ninguna correlación
#     puede dar. Sustituye a la prueba anterior de quitar |am| < 0.2 de un
#     regulón inferido, que era tautológica en la parte dirigida.
# ============================================================
syn_files <- list.files(dirs$synthetic, "^evaluation\\.rds$", recursive = TRUE, full.names = TRUE)
a4 <- do.call(rbind, lapply(syn_files, function(f) readRDS(f)$nonmono_nes_truth))
a4s <- aggregate(cbind(rho_vs_truth_full, top20_vs_truth_full, precision_top10) ~ regime + n + variant, a4, mean)
write.table(a4s, file.path(out, "appendix_A4_truth_regulon_nonmono.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_len(nrow(a4s))) for (v in c("rho_vs_truth_full", "top20_vs_truth_full", "precision_top10"))
  add_num(sprintf("A4_%s_n%d_%s_%s", a4s$regime[i], a4s$n[i], a4s$variant[i], v), a4s[[v]][i], "A4: truth regulon variant, mean over replicates")
mrfull <- do.call(rbind, lapply(syn_files, function(f) { x <- readRDS(f); m <- x$master_regulators; m <- m[m$regulon == "truth" & m$min_targets == nb$narnea$min_targets, ]; cbind(x$meta, precision_top10 = m$precision_top10) }))
mf <- aggregate(precision_top10 ~ regime + n, mrfull, mean)
for (i in seq_len(nrow(mf))) add_num(sprintf("A4_%s_n%d_truth_full_precision_top10", mf$regime[i], mf$n[i]), mf$precision_top10[i], "A4: truth regulon, full, mean fraction of true MRs in top-10")
a4p <- rbind(data.frame(regime = mf$regime, n = mf$n, variant = "truth_full", precision_top10 = mf$precision_top10),
             a4s[, c("regime", "n", "variant", "precision_top10")])
a4p$variant <- factor(a4p$variant, c("truth_full", "truth_no_nonmono", "truth_nonmono_signed"),
                      c("True network, non-monotone edges with mode 0", "True network without non-monotone edges", "True network, non-monotone edges with perfect sign"))
fA4 <- ggplot(a4p, aes(n, precision_top10, colour = variant)) + geom_line() + geom_point() + facet_wrap(~regime) + scale_x_log10() +
  labs(x = "Samples (n, log scale)", y = "Fraction of true MRs in top-10 |NES|", colour = NULL) + guides(colour = guide_legend(ncol = 1))
save_fig(fA4, "figA4_nonmonotone_contribution_to_nes", 9, 4.5)

write.table(do.call(rbind, nums), file.path(out, "appendix_numbers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("Apéndice: figuras en ", fig_dir, "; números en appendix_numbers.tsv")
