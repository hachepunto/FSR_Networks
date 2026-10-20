# Evaluación de un dataset sintético contra la red verdadera.
#
# Regulones: current, A_k<k> (corto laxo, top-k), A_r<r> (cada umbral de
# r_min_grid, top-k_main), A_tri_k<main> (DPI por triángulos, paso 14, si existe), aracne3_full,
# aracne3_c<m>, B_k<main> (ARACNe3 top-k) y truth (la red verdadera, techo
# de referencia para la recuperación de reguladores maestros).
#
# Aristas (por regulón): precisión, recall total y por tipo de relación
# (lineal / monótona no lineal / no monótona), signo correcto en las
# monótonas recuperadas, y falsos positivos "indirectos" (el TF comparte
# módulo con un regulador verdadero del target: lo que el DPI debe podar),
# en número absoluto y por cada 100 aristas; la fracción sobre los FP se
# conserva pero su denominador depende del método. Además, precisión promedio (AP) del ranking completo de
# cada método (likelihood y desempate), independiente de k.
#
# Reguladores maestros: NaRnEA con la firma DESeq2 del dataset, con cada
# mínimo de targets de synthetic.narnea_min_targets. Métrica principal:
# precisión en el top-10 por |NES| (con 10 MR verdaderos = fracción de MR
# recuperados) y AUROC sobre todos los TFs (no evaluable = no detectado).
# Secundarias: AUROC solo entre TFs evaluables, signo de NES y cuántos MR
# fueron evaluables.
#
# Uso: Rscript scripts/network_benchmark/07_evaluate_synthetic.R <dataset_dir>
# Salida: <dataset_dir>/evaluation.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/07_evaluate_synthetic.R <dataset_dir>")
ds <- normalizePath(args[1])
name <- basename(ds)

config <- bench_config()
nb <- config$network_benchmark
k_main <- nb$top_k_main
min_count <- nb$aracne3_count_filter
n_subnets <- config$network_inference$aracne3$subnets

truth <- readRDS(file.path(ds, "truth.rds"))
te <- truth$edges
te$key <- paste(te$tf, te$target)
ges <- readRDS(file.path(ds, "signature.rds"))
mr <- truth$master_regulators
module <- truth$tf_module

# ── Regulones ──
lenients <- lapply(setNames(unlist(nb$corto$r_min_grid), lenient_name(unlist(nb$corto$r_min_grid))),
                   function(r) readRDS(file.path(ds, paste0("corto_", lenient_name(r), ".rds"))))
lenient <- lenients[[lenient_name(nb$corto$r_min)]]
net <- read.delim(file.path(ds, "aracne3", paste0("consolidated-net_", name, ".tsv")))
regs <- list(current = readRDS(file.path(ds, "corto_current.rds")))
for (k in unlist(nb$top_k)) regs[[paste0("A_k", k)]] <- topk_regulon(lenient, k)
for (r in unlist(nb$corto$r_min_grid)) regs[[sprintf("A_r%.1f", r)]] <- topk_regulon(lenients[[lenient_name(r)]], k_main)
# FSR con DPI por triángulos (paso 14), si existe
tri_file <- file.path(ds, paste0("corto_", lenient_name(nb$corto$r_min), "_triangle.rds"))
if (file.exists(tri_file)) regs[[paste0("A_tri_k", k_main)]] <- topk_regulon(readRDS(tri_file), k_main)
regs$aracne3_full <- aracne3_to_regulon(net, n_subnets)
regs[[paste0("aracne3_c", min_count)]] <- aracne3_to_regulon(net, n_subnets, min_count)
regs[[paste0("B_k", k_main)]] <- topk_regulon(regs$aracne3_full, k_main)
# B recortada al MISMO número de aristas que A (por soporte, luego MI), para
# que las métricas por arista (precisión, FP indirectos por 100) comparen
# a igual tamaño: A deja menos de k targets por regulador cuando su DPI poda.
n_edges_A <- sum(vapply(regs[[paste0("A_k", k_main)]], function(x) length(x$tfmode), integer(1)))
B_edges <- regulon_edges(regs$aracne3_full)
if (nrow(B_edges) > 0) {
  B_tb <- unlist(lapply(regs$aracne3_full, tiebreak_of), use.names = FALSE)
  keep <- order(B_edges$likelihood, B_tb, decreasing = TRUE)[seq_len(min(n_edges_A, nrow(B_edges)))]
  B_keep <- B_edges[keep, ]
  regs$B_matched <- lapply(split(B_keep, B_keep$tf), function(d) list(tfmode = setNames(d$mode, d$target), likelihood = setNames(d$likelihood, d$target)))
} else {
  regs$B_matched <- list()   # red ARACNe3 vacía: B emparejada también vacía
}

# Red verdadera en formato interno: modo = signo (0 en las no monótonas,
# cuyo efecto no tiene dirección), likelihood = 1.
truth_reg <- lapply(split(te, te$tf), function(d) {
  list(tfmode = setNames(ifelse(d$type == "non_monotone", 0, d$sign), d$target),
       likelihood = setNames(rep(1, nrow(d)), d$target))
})

# ── Métricas de aristas ──
# Reguladores verdaderos de cada target, para detectar FP indirectos.
true_regs_by_target <- split(te$tf, te$target)

edge_metrics <- function(reg) {
  e <- regulon_edges(reg)
  key <- paste(e$tf, e$target)
  tp <- key %in% te$key
  hit <- te$key %in% key
  fp <- e[!tp, ]
  indirect <- if (nrow(fp)) mapply(function(tf, tg) {
    tr <- true_regs_by_target[[tg]]
    !is.null(tr) && any(module[tr] == module[tf])
  }, fp$tf, fp$target, USE.NAMES = FALSE) else logical(0)
  mono <- tp & te$type[match(key, te$key)] != "non_monotone"
  data.frame(
    n_edges = nrow(e), precision = if (nrow(e)) mean(tp) else NA, recall = mean(hit),
    recall_linear = mean(hit[te$type == "linear"]),
    recall_monotone_nonlinear = mean(hit[te$type == "monotone_nonlinear"]),
    recall_non_monotone = mean(hit[te$type == "non_monotone"]),
    sign_accuracy_monotone = if (any(mono)) mean(sign(e$mode[mono]) == te$sign[match(key[mono], te$key)]) else NA,
    n_fp = nrow(fp), n_fp_indirect = sum(indirect),
    fp_indirect_per_100_edges = if (nrow(e)) 100 * sum(indirect) / nrow(e) else NA,
    frac_fp_indirect = if (nrow(fp)) mean(indirect) else NA)
}

# AP del ranking completo de un método: aristas ordenadas por likelihood
# y desempate; las verdaderas no reportadas nunca se alcanzan.
average_precision <- function(reg) {
  e <- regulon_edges(reg)
  if (nrow(e) == 0) return(0)
  tb <- unlist(lapply(reg, tiebreak_of), use.names = FALSE)
  o <- order(e$likelihood, tb, decreasing = TRUE)
  tp <- (paste(e$tf, e$target) %in% te$key)[o]
  sum((cumsum(tp) / seq_along(tp))[tp]) / nrow(te)
}

edges <- do.call(rbind, lapply(names(regs), function(nm) cbind(regulon = nm, edge_metrics(regs[[nm]]))))

# Control positivo de la métrica de no monotonía de 04_evaluate_real.R: la
# misma estadística sobre la red ARACNe3 sintética, donde SÍ hay aristas no
# monótonas verdaderas, más cuántas de ellas recupera y con qué |Spearman|.
nonmono <- do.call(rbind, lapply(c(1, min_count, 50), function(m) {
  d <- net[net$count.values >= m, ]
  if (!nrow(d)) return(NULL)
  key <- paste(d$regulator.values, d$target.values)
  tp_nm <- key %in% te$key[te$type == "non_monotone"]
  hi_mi <- d$mi.values >= quantile(d$mi.values, 0.75)
  data.frame(regulon = sprintf("aracne3_c%d", m), min_count = m, n_edges = nrow(d),
             frac_abs_scc_lt_0.2 = mean(abs(d$scc.values) < 0.2),
             spearman_mi_vs_abs_scc = cor(d$mi.values, abs(d$scc.values), method = "spearman"),
             n_top_quartile_mi = sum(hi_mi), n_top_quartile_mi_abs_scc_lt_0.2 = sum(hi_mi & abs(d$scc.values) < 0.2),
             n_true_non_monotone_recovered = sum(tp_nm),
             frac_true_non_monotone_abs_scc_lt_0.2 = if (any(tp_nm)) mean(abs(d$scc.values[tp_nm]) < 0.2) else NA,
             frac_edges_true_non_monotone = mean(tp_nm))
}))
eB <- regulon_edges(regs[[paste0("B_k", k_main)]])
if (nrow(eB)) nonmono <- rbind(nonmono, data.frame(regulon = paste0("B_k", k_main), min_count = NA, n_edges = nrow(eB),
  frac_abs_scc_lt_0.2 = mean(abs(eB$mode) < 0.2), spearman_mi_vs_abs_scc = NA, n_top_quartile_mi = NA,
  n_top_quartile_mi_abs_scc_lt_0.2 = NA,
  n_true_non_monotone_recovered = sum(paste(eB$tf, eB$target) %in% te$key[te$type == "non_monotone"]),
  frac_true_non_monotone_abs_scc_lt_0.2 = NA,
  frac_edges_true_non_monotone = mean(paste(eB$tf, eB$target) %in% te$key[te$type == "non_monotone"])))
ap <- data.frame(method = c(paste0("corto_", names(lenients)), "aracne3_full"),
                 average_precision = c(vapply(lenients, average_precision, numeric(1)), average_precision(regs$aracne3_full)))

# ── Contribución de las aristas no monótonas al NES, con la red VERDADERA ──
# Tres versiones del regulón verdadero: (i) completo con am = signo (0 en las
# no monótonas), (ii) sin las aristas no monótonas, (iii) no monótonas con
# |am| = 1 (un "signo correcto" que ninguna correlación puede dar; cota
# superior de lo que aportarían si el scoring las aprovechara).
truth_variants <- list(
  truth_full = truth_reg,
  truth_no_nonmono = lapply(split(te[te$type != "non_monotone", ], te$tf[te$type != "non_monotone"]), function(d)
    list(tfmode = setNames(d$sign, d$target), likelihood = setNames(rep(1, nrow(d)), d$target))),
  truth_nonmono_signed = lapply(split(te, te$tf), function(d)
    list(tfmode = setNames(d$sign, d$target), likelihood = setNames(rep(1, nrow(d)), d$target))))
truth_res <- lapply(truth_variants, function(r) run_narnea(to_narnea(r, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed))
nonmono_nes <- do.call(rbind, lapply(names(truth_res)[-1], function(v) {
  cc <- mr_concordance(truth_res$truth_full, truth_res[[v]])
  pos <- truth_res[[v]]$tf %in% mr$tf
  data.frame(variant = v, rho_vs_truth_full = cc["rho"], top20_vs_truth_full = cc["top20"],
             precision_top10 = sum(truth_res[[v]]$tf[order(-abs(truth_res[[v]]$nes))][1:10] %in% mr$tf) / nrow(mr))
}))

# ── Reguladores maestros ──
auroc <- function(score, positive) {
  r <- rank(score)
  np <- sum(positive); nn <- sum(!positive)
  (sum(r[positive]) - np * (np + 1) / 2) / (np * nn)
}

regs$truth <- truth_reg
mr_grid <- expand.grid(regulon = names(regs), min_targets = unlist(nb$synthetic$narnea_min_targets),
                       stringsAsFactors = FALSE)
mr_metrics <- do.call(rbind, lapply(seq_len(nrow(mr_grid)), function(i) {
  nm <- mr_grid$regulon[i]; mt <- mr_grid$min_targets[i]
  res <- run_narnea(to_narnea(regs[[nm]], "rank"), ges, mt, nb$narnea$seed)
  if (is.null(res)) return(data.frame(regulon = nm, min_targets = mt, n_regulators = 0, mr_evaluated = 0,
                                      auroc = NA, auroc_all = 0.5, precision_top10 = 0, mr_sign_correct = NA))
  pos <- res$tf %in% mr$tf
  # fracción de los n_MR verdaderos dentro del top-10 por |NES| (denominador
  # fijo: un MR no evaluable cuenta como no recuperado)
  top10 <- res$tf[order(-abs(res$nes))][1:min(10, nrow(res))]
  sgn <- res$nes[match(mr$tf, res$tf)]
  # AUROC sobre todos los TFs: un TF no evaluable (regulón < min_targets)
  # recibe |NES| = 0, así que un MR no evaluable cuenta como no detectado.
  all_score <- setNames(rep(0, length(module)), names(module))
  all_score[res$tf] <- abs(res$nes)
  data.frame(regulon = nm, min_targets = mt, n_regulators = nrow(res), mr_evaluated = sum(pos),
             auroc = if (sum(pos) > 0 && sum(!pos) > 0) auroc(abs(res$nes), pos) else NA,
             auroc_all = auroc(all_score, names(all_score) %in% mr$tf),
             precision_top10 = sum(top10 %in% mr$tf) / nrow(mr),
             mr_sign_correct = mean(sign(sgn) == mr$sign, na.rm = TRUE))
}))

meta <- data.frame(dataset = name, regime = truth$regime, n = truth$n_network, replicate = truth$replicate)
saveRDS(list(meta = meta, edges = cbind(meta, edges), average_precision = cbind(meta, ap),
             master_regulators = cbind(meta, mr_metrics),
             nonmonotonic = if (is.null(nonmono)) NULL else cbind(meta, nonmono),
             nonmono_nes_truth = cbind(meta, nonmono_nes)),
        file.path(ds, "evaluation.rds"))
message(name, ": evaluación guardada")
print(mr_metrics, row.names = FALSE)
