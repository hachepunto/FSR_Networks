# Evalúa las redes del cáncer completo (13) con NaRnEA y la firma de DESeq2.
#
# 1. Sustituto (cohortes grandes X reducidas a small_n): concordancia con la
#    verdad (red de X con todas sus muestras) de
#      pooled  — otros subtipos + small_n de X
#      own     — solo los small_n de X (mismos tumores)
#      swap    — red completa de OTRO subtipo del mismo cáncer (media sobre
#                subtipos; línea base sin tumores de X)
# 2. Todas las cohortes reales: NaRnEA con la red del cáncer completo frente a
#    su red propia, top-20 de ambas (production = las de n < min_n).
# 3. Especificidad: top-20 compartido entre subtipos del mismo cáncer cuando
#    todos se interrogan con la red pooled frente a cada uno con la suya. Si
#    la red pooled homogeneiza los subtipos, la especificidad recae solo en
#    la firma.
#
# Uso: Rscript scripts/network_benchmark/13b_evaluate_pooled.R
# Salida: <out_dir>/real/evaluation/pooled.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
pp <- nb$pooled
pool_dir <- file.path(dirs$out, "real", "pooled")
at <- analysis_table(config)
min_n <- network_min_n(config)

ref_reg <- function(id) topk_regulon(readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds"))), nb$top_k_main)
narnea <- function(reg, ges) run_narnea(to_narnea(reg, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed)
n_eval <- function(res) if (is.null(res)) 0L else nrow(res)
top_tfs <- function(res, k = 20) if (is.null(res)) character(0) else res$tf[order(-abs(res$nes))][seq_len(min(k, nrow(res)))]

# ── 1. Sustituto ──
sur <- list()
for (id in unlist(pp$cohorts)) {
  cancer <- at$cancer[at$analysis_id == id]
  ges <- read_deseq2_stat(config, id)
  truth <- narnea(ref_reg(id), ges)
  others <- setdiff(at$analysis_id[at$cancer == cancer], id)
  swap <- sapply(others, function(o) mr_concordance(truth, narnea(ref_reg(o), ges)))
  for (r in seq_len(pp$replicates)) {
    stem <- file.path(pool_dir, sprintf("%s_rep%d", id, r))
    if (!file.exists(paste0(stem, "_pooled.rds"))) next
    for (kind in c("pooled", "own")) {
      reg <- readRDS(paste0(stem, "_", kind, ".rds"))
      res <- narnea(reg, ges)
      cc <- mr_concordance(truth, res)
      sur[[length(sur) + 1]] <- data.frame(
        analysis_id = id, cancer = cancer, n_full = at$n_red[at$analysis_id == id], replicate = r,
        network = kind, n_network = attr(reg, "params")$n, n_regulators = length(reg), n_evaluable = n_eval(res),
        rho = cc["rho"], top20 = cc["top20"], top50 = cc["top50"],
        swap_top20_mean = mean(swap["top20", ], na.rm = TRUE), swap_rho_mean = mean(swap["rho", ], na.rm = TRUE),
        n_swap = length(others))
    }
  }
}
sur <- do.call(rbind, sur); rownames(sur) <- NULL

# ── 2 y 3. Todas las cohortes: TMRs con la red propia vs con la red del cáncer completo ──
# (la red del cáncer completo incluye a la propia cohorte; su peso en ella es frac_of_pool)
small <- at$analysis_id[at$n_red < min_n]
cancers <- unique(at$cancer)
allc <- list(); spec <- list(); tops <- list()
for (cancer in cancers) {
  f <- file.path(pool_dir, sprintf("%s_pooled_full.rds", cancer))
  if (!file.exists(f)) { message(cancer, ": falta la red del cáncer completo; se omite"); next }
  preg <- readRDS(f)
  ids <- at$analysis_id[at$cancer == cancer]
  res_pool <- res_own <- list()
  for (id in ids) {
    ges <- read_deseq2_stat(config, id)
    res_pool[[id]] <- narnea(preg, ges)
    res_own[[id]] <- narnea(ref_reg(id), ges)
    cc <- mr_concordance(res_own[[id]], res_pool[[id]])
    n_red <- at$n_red[at$analysis_id == id]
    allc[[length(allc) + 1]] <- data.frame(
      analysis_id = id, cancer = cancer, n_red = n_red, n_pooled = attr(preg, "params")$n,
      frac_of_pool = n_red / attr(preg, "params")$n, tier = analysis_tier(n_red, config),
      n_evaluable_own = n_eval(res_own[[id]]), n_evaluable_pooled = n_eval(res_pool[[id]]),
      rho_own_vs_pooled = cc["rho"], top20_own_vs_pooled = cc["top20"], top50_own_vs_pooled = cc["top50"])
    tops[[id]] <- list(own = top_tfs(res_own[[id]]), pooled = top_tfs(res_pool[[id]]))
  }
  for (pair in combn(ids, 2, simplify = FALSE)) {
    spec[[length(spec) + 1]] <- data.frame(
      cancer = cancer, x = pair[1], y = pair[2],
      any_small = any(pair %in% small),
      top20_own_networks = mr_concordance(res_own[[pair[1]]], res_own[[pair[2]]])["top20"],
      top20_pooled_network = mr_concordance(res_pool[[pair[1]]], res_pool[[pair[2]]])["top20"])
  }
}
allc <- do.call(rbind, allc); spec <- do.call(rbind, spec)
rownames(allc) <- rownames(spec) <- NULL
prod <- allc[allc$analysis_id %in% small, ]

saveRDS(list(surrogate = sur, production = prod, all_cohorts = allc, specificity = spec, top20 = tops),
        file.path(dirs$real_eval, "pooled.rds"))

message("\nSustituto — mediana del top-20 con la verdad (cohortes x réplicas):")
print(aggregate(cbind(top20, rho, n_evaluable) ~ network, sur, median), digits = 3)
message("Red de otro subtipo (media por cohorte), mediana: ", round(median(tapply(sur$swap_top20_mean, sur$analysis_id, mean)), 3))
per_coh <- aggregate(top20 ~ analysis_id + network, sur, mean)
w <- reshape(per_coh, idvar = "analysis_id", timevar = "network", direction = "wide")
w$swap <- tapply(sur$swap_top20_mean, sur$analysis_id, mean)[w$analysis_id]
print(w, digits = 2)
message("Pareado pooled - own por cohorte: Wilcoxon p = ",
        signif(wilcox.test(w$top20.pooled, w$top20.own, paired = TRUE, exact = FALSE)$p.value, 3))
message("\nTodas las cohortes, TMRs red propia vs red del cáncer completo:"); print(allc, digits = 3)
message("\nEspecificidad (top-20 entre subtipos del mismo cáncer):"); print(spec, digits = 2)
