# Evaluación de un analysis_id real: estructura de los regulones,
# concordancia con CollecTRI, NaRnEA con la firma del pipeline (stat de
# DESeq2) y concordancia de reguladores maestros entre esquemas.
#
# Regulones evaluados:
#   current      corto como está hoy en el pipeline (p fijo sobre log2TPM)
#   A_k<k>       esquema A: corto laxo sobre rangos (|r| >= r_min), top-k
#   A_r<r>       esquema A con cada umbral de r_min_grid, top-k_main
#   A_k<main>_awlik  igual que A_k<main> pero aw = likelihood cruda
#   ARACNe3 (solo si network_method == ARACNe3):
#     aracne3_full   red consolidada completa
#     aracne3_c<m>   aristas con soporte >= m subredes
#     B_k<main>      esquema B: ARACNe3 recortado a top-k (soporte, luego MI)
#   CollecTRI      regulón de referencia (curado)
#
# La comparación central es A_k<main> vs B_k<main>, contrastada con la
# variabilidad interna de ARACNe3 (B_k<main> vs aracne3_c<m>): dos
# filtrados razonables de la MISMA red.
#
# Uso: Rscript scripts/network_benchmark/04_evaluate_real.R <analysis_id>
# Salida: <out_dir>/real/evaluation/<analysis_id>.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/04_evaluate_real.R <analysis_id>")
analysis_id <- args[1]

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
k_main <- nb$top_k_main
min_count <- nb$aracne3_count_filter
n_subnets <- config$network_inference$aracne3$subnets

at <- analysis_table(config)
row <- at[at$analysis_id == analysis_id, ]
if (nrow(row) != 1) stop("analysis_id desconocido en config.yaml: ", analysis_id)

universe <- rownames(read_tumor_expr(config, analysis_id, row$cancer))
ges <- read_deseq2_stat(config, analysis_id)
ct <- read_collectri(dirs)

# ── Regulones ──
current <- readRDS(file.path(dirs$corto, paste0(analysis_id, "_current.rds")))
lenients <- lapply(setNames(unlist(nb$corto$r_min_grid), lenient_name(unlist(nb$corto$r_min_grid))),
                   function(r) readRDS(file.path(dirs$corto, paste0(analysis_id, "_", lenient_name(r), ".rds"))))
lenient <- lenients[[lenient_name(nb$corto$r_min)]]

regs <- list(current = current)
for (k in unlist(nb$top_k)) regs[[paste0("A_k", k)]] <- topk_regulon(lenient, k)
for (r in unlist(nb$corto$r_min_grid)) regs[[sprintf("A_r%.1f", r)]] <- topk_regulon(lenients[[lenient_name(r)]], k_main)

has_aracne3 <- row$network_method == "ARACNe3"
if (has_aracne3) {
  net <- read_aracne3_net(dirs, analysis_id)
  regs$aracne3_full <- aracne3_to_regulon(net, n_subnets)
  regs[[paste0("aracne3_c", min_count)]] <- aracne3_to_regulon(net, n_subnets, min_count)
  regs[[paste0("B_k", k_main)]] <- topk_regulon(regs$aracne3_full, k_main)
}

# ── Estructura ──
regulon_structure_tbl <- do.call(rbind, lapply(names(regs), function(nm) {
  s <- regulon_structure(regs[[nm]], universe, ct, nb$narnea$min_targets, names(ges))
  if (is.null(s)) return(NULL)
  cbind(analysis_id = analysis_id, regulon = nm, s)
}))

# ── NaRnEA ──
narnea <- lapply(regs, function(r) run_narnea(to_narnea(r, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed))
awlik <- paste0("A_k", k_main, "_awlik")
narnea[[awlik]] <- run_narnea(to_narnea(regs[[paste0("A_k", k_main)]], "likelihood"), ges,
                              nb$narnea$min_targets, nb$narnea$seed)
narnea$CollecTRI <- run_narnea(collectri_regulon(ct), ges, nb$narnea$min_targets, nb$narnea$seed)

# ── Concordancia entre esquemas ──
A <- paste0("A_k", k_main); B <- paste0("B_k", k_main); C <- paste0("aracne3_c", min_count)
pairs <- c(list(c(A, "current"), c(A, awlik), c(A, "CollecTRI"), c("current", "CollecTRI")),
           lapply(setdiff(unlist(nb$top_k), k_main), function(k) c(A, paste0("A_k", k))))
Ar <- sprintf("A_r%.1f", unlist(nb$corto$r_min_grid))
pairs <- c(pairs, lapply(Ar, function(a) c(a, "CollecTRI")), lapply(Ar, function(a) c(a, "current")))
if (has_aracne3) pairs <- c(pairs, list(c(A, B), c(B, C), c(B, "aracne3_full"), c(A, C),
                                        c(B, "CollecTRI"), c("aracne3_full", "CollecTRI")),
                            lapply(Ar, function(a) c(a, B)))
concordance <- do.call(rbind, lapply(pairs, function(p) {
  cc <- mr_concordance(narnea[[p[1]]], narnea[[p[2]]])
  data.frame(analysis_id = analysis_id, x = p[1], y = p[2], t(cc), check.names = FALSE)
}))

# ── No monotonía en ARACNe3: aristas con MI alta y |Spearman| baja ──
# Misma métrica que 07_evaluate_synthetic.R calcula sobre las redes ARACNe3
# sintéticas (control positivo con 20 % de aristas no monótonas verdaderas).
nonmono <- NULL
if (has_aracne3) {
  nonmono <- do.call(rbind, lapply(c(1, min_count, 50), function(m) {
    d <- net[net$count.values >= m, ]
    hi_mi <- d$mi.values >= quantile(d$mi.values, 0.75)
    data.frame(analysis_id = analysis_id, regulon = sprintf("aracne3_c%d", m), min_count = m, n_edges = nrow(d),
               frac_abs_scc_lt_0.2 = mean(abs(d$scc.values) < 0.2),
               spearman_mi_vs_abs_scc = cor(d$mi.values, abs(d$scc.values), method = "spearman"),
               n_top_quartile_mi = sum(hi_mi),
               n_top_quartile_mi_abs_scc_lt_0.2 = sum(hi_mi & abs(d$scc.values) < 0.2))
  }))
  eB <- regulon_edges(regs[[B]])
  nonmono <- rbind(nonmono, data.frame(analysis_id = analysis_id, regulon = B, min_count = NA, n_edges = nrow(eB),
                                       frac_abs_scc_lt_0.2 = mean(abs(eB$mode) < 0.2),
                                       spearman_mi_vs_abs_scc = NA, n_top_quartile_mi = NA,
                                       n_top_quartile_mi_abs_scc_lt_0.2 = NA))
}

saveRDS(list(analysis_id = analysis_id, n = row$n_red, network_method = row$network_method,
             structure = regulon_structure_tbl, narnea = narnea, concordance = concordance, nonmonotonic = nonmono,
             corto_params = c(list(current = attr(current, "params")), lapply(lenients, attr, "params"))),
        file.path(dirs$real_eval, paste0(analysis_id, ".rds")))

message(analysis_id, ": evaluación guardada")
print(concordance[, c("x", "y", "n_shared", "rho", "top20", "top50")], row.names = FALSE)
