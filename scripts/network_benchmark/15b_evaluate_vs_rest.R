# Compara, con la red FSR propia de cada cohorte, los TMRs de la firma
# tumor vs control (run_deseq2.R) con los de la firma subtipo vs resto (15).
#
# Qué se mira (fijado antes de correr):
#   1. Controles positivos de biología conocida: rango (por |NES|) y signo de
#      TFs que definen subtipos. Si el contraste aporta, suben con el signo
#      esperado respecto a tumor vs control.
#   2. Especificidad: top-20 compartido entre subtipos del mismo cáncer; si el
#      contraste aísla lo propio del subtipo, baja respecto a tumor vs control.
#   3. Cuánto cambia el top-20 de cada cohorte entre las dos firmas.
#
# Uso: Rscript scripts/network_benchmark/15b_evaluate_vs_rest.R
# Salida: <out_dir>/real/evaluation/vs_rest.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
vr_dir <- file.path(dirs$out, "real", "vs_rest")

controls <- list(
  BRCA_basal = c(ESR1 = -1, FOXA1 = -1, GATA3 = -1, XBP1 = -1, SPDEF = -1),
  BRCA_luma  = c(ESR1 = 1, FOXA1 = 1, GATA3 = 1, XBP1 = 1, SPDEF = 1),
  COAD_cms4  = c(TWIST1 = 1, PRRX1 = 1, ZEB1 = 1, ZEB2 = 1, SNAI2 = 1))

ids <- sub("_stat\\.rds$", "", list.files(vr_dir, "_stat\\.rds$"))
at <- analysis_table(config)
own_reg <- function(id) topk_regulon(readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds"))), nb$top_k_main)
nar <- function(reg, ges) run_narnea(to_narnea(reg, "rank"), ges, nb$narnea$min_targets, nb$narnea$seed)

res <- list()
for (id in ids) {
  reg <- own_reg(id)
  res[[id]] <- list(vs_normal = nar(reg, read_deseq2_stat(config, id)),
                    vs_rest = nar(reg, readRDS(file.path(vr_dir, paste0(id, "_stat.rds")))))
}

# 3. Cambio por cohorte
per <- do.call(rbind, lapply(ids, function(id) {
  cc <- mr_concordance(res[[id]]$vs_normal, res[[id]]$vs_rest)
  data.frame(analysis_id = id, cancer = at$cancer[at$analysis_id == id], n = at$n_red[at$analysis_id == id],
             top20_normal_vs_rest = cc["top20"], rho_normal_vs_rest = cc["rho"])
}))

# 2. Especificidad entre subtipos del mismo cáncer
pairs <- do.call(rbind, lapply(split(ids, at$cancer[match(ids, at$analysis_id)]), function(v) {
  if (length(v) < 2) return(NULL)
  do.call(rbind, lapply(combn(v, 2, simplify = FALSE), function(p) data.frame(
    cancer = at$cancer[at$analysis_id == p[1]], x = p[1], y = p[2],
    top20_vs_normal = mr_concordance(res[[p[1]]]$vs_normal, res[[p[2]]]$vs_normal)["top20"],
    top20_vs_rest = mr_concordance(res[[p[1]]]$vs_rest, res[[p[2]]]$vs_rest)["top20"])))
}))

# 1. Controles positivos
rank_of <- function(r, tf) { o <- r$tf[order(-abs(r$nes))]; i <- match(tf, o); i }
ctrl <- do.call(rbind, lapply(intersect(names(controls), ids), function(id) {
  do.call(rbind, lapply(names(controls[[id]]), function(tf) {
    do.call(rbind, lapply(c("vs_normal", "vs_rest"), function(sig) {
      r <- res[[id]][[sig]]
      data.frame(analysis_id = id, tf = tf, signature = sig, expected_sign = controls[[id]][[tf]],
                 n_evaluable = nrow(r), rank = rank_of(r, tf),
                 nes = if (tf %in% r$tf) r$nes[r$tf == tf] else NA)
    }))
  }))
}))
ctrl$sign_ok <- sign(ctrl$nes) == ctrl$expected_sign

saveRDS(list(per_cohort = per, pairs = pairs, controls = ctrl, results = res), file.path(dirs$real_eval, "vs_rest.rds"))

message("Cambio del top-20 por cohorte (tumor vs control frente a subtipo vs resto):"); print(per, digits = 2, row.names = FALSE)
message("Mediana: ", median(per$top20_normal_vs_rest))
message("\nEspecificidad entre subtipos (mediana top-20): vs control ", median(pairs$top20_vs_normal),
        "; vs resto ", median(pairs$top20_vs_rest), " (", nrow(pairs), " pares)")
print(aggregate(cbind(top20_vs_normal, top20_vs_rest) ~ cancer, pairs, median), digits = 2)
message("\nControles positivos (rango por |NES| y signo):"); print(ctrl, digits = 3, row.names = FALSE)
