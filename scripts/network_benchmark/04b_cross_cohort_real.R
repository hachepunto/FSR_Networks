# Dos referencias para interpretar la concordancia A vs B (04_evaluate_real.R):
#
#   null   — línea base nula: NaRnEA con la firma de una cohorte y el regulón
#            de FSR de una cohorte de OTRO cáncer. Mide cuánta
#            concordancia se obtiene con una red que no corresponde a la
#            cohorte (la firma domina parte del NES).
#   within — efecto biológico: NaRnEA con el regulón A de cada cohorte,
#            comparando subtipos del MISMO cáncer (firma y red propias de
#            cada uno). Es la magnitud de las diferencias que el pipeline
#            quiere medir; el efecto del método (A vs B) se compara contra ella.
#   within_net — descomposición: firma propia, red A del otro subtipo del
#            mismo cáncer. Es el efecto de cambiar SOLO la red entre
#            subtipos: el yardstick emparejado con A vs B (que también
#            cambia solo la red, de A a B).
#   within_sig — firma del otro subtipo, red propia: efecto de cambiar
#            solo la firma.
#
# Uso: Rscript scripts/network_benchmark/04b_cross_cohort_real.R
# Requiere las evaluaciones de 04 para todos los analysis_id.
# Salida: <out_dir>/real/evaluation/cross_cohort.rds

suppressMessages(library(NaRnEA))
source("scripts/network_benchmark/lib_benchmark.R")
apply_narnea_overflow_patch()

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
A <- paste0("A_k", nb$top_k_main)
at <- analysis_table(config)

ev <- lapply(setNames(at$analysis_id, at$analysis_id), function(id) {
  f <- file.path(dirs$real_eval, paste0(id, ".rds"))
  if (!file.exists(f)) stop("Falta la evaluación de ", id, " (04_evaluate_real.R)")
  readRDS(f)
})
regA <- lapply(at$analysis_id, function(id) {
  lenient <- readRDS(file.path(dirs$corto, paste0(id, "_", lenient_name(nb$corto$r_min), ".rds")))
  to_narnea(topk_regulon(lenient, nb$top_k_main), "rank")
})
names(regA) <- at$analysis_id
ges <- lapply(setNames(at$analysis_id, at$analysis_id), function(id) read_deseq2_stat(config, id))
own <- lapply(ev, function(x) x$narnea[[A]])

# ── within: subtipos del mismo cáncer, cada uno con su red y su firma ──
within <- do.call(rbind, lapply(split(at$analysis_id, at$cancer), function(ids) {
  if (length(ids) < 2) return(NULL)
  pairs <- combn(ids, 2)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    cc <- mr_concordance(own[[pairs[1, j]]], own[[pairs[2, j]]])
    data.frame(cancer = sub("_.*", "", pairs[1, j]), x = pairs[1, j], y = pairs[2, j], t(cc), check.names = FALSE)
  }))
}))

# ── null: firma propia, red A de una cohorte de otro cáncer ──
# Para cada cohorte se eligen n_null cohortes de otros cánceres con
# semilla fija; la comparación es contra el resultado propio (firma + red propias).
# Las cohortes nulas se eligen todas ANTES de llamar a NaRnEA, que fija su
# propia semilla internamente y alteraría el estado del RNG entre cohortes.
set.seed(nb$narnea$seed)
n_null <- nb$null_networks_per_cohort
null_choice <- lapply(setNames(at$analysis_id, at$analysis_id), function(id)
  sample(at$analysis_id[at$cancer != at$cancer[at$analysis_id == id]], n_null))
null <- do.call(rbind, lapply(at$analysis_id, function(id) {
  others <- null_choice[[id]]
  do.call(rbind, lapply(others, function(o) {
    res <- run_narnea(regA[[o]], ges[[id]], nb$narnea$min_targets, nb$narnea$seed)
    cc <- mr_concordance(own[[id]], res)
    data.frame(analysis_id = id, n = at$n_red[at$analysis_id == id], network_from = o, t(cc), check.names = FALSE)
  }))
}))

# ── descomposición firma / red entre subtipos del mismo cáncer ──
decomp <- do.call(rbind, lapply(split(at$analysis_id, at$cancer), function(ids) {
  if (length(ids) < 2) return(NULL)
  pairs <- combn(ids, 2)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]
    do.call(rbind, lapply(list(c(a, b), c(b, a)), function(p) {
      net_swapped <- run_narnea(regA[[p[2]]], ges[[p[1]]], nb$narnea$min_targets, nb$narnea$seed)  # firma propia, red ajena
      sig_swapped <- run_narnea(regA[[p[1]]], ges[[p[2]]], nb$narnea$min_targets, nb$narnea$seed)  # red propia, firma ajena
      rbind(data.frame(comparison = "within_net", cancer = sub("_.*", "", p[1]), x = p[1], y = p[2], t(mr_concordance(own[[p[1]]], net_swapped)), check.names = FALSE),
            data.frame(comparison = "within_sig", cancer = sub("_.*", "", p[1]), x = p[1], y = p[2], t(mr_concordance(own[[p[1]]], sig_swapped)), check.names = FALSE))
    }))
  }))
}))

saveRDS(list(within = within, null = null, decomp = decomp), file.path(dirs$real_eval, "cross_cohort.rds"))
for (cmp in c("within_net", "within_sig")) {
  d <- decomp[decomp$comparison == cmp, ]
  message(cmp, " (", nrow(d), " pares dirigidos): mediana rho ", round(median(d$rho), 3), " top20 ", round(median(d$top20, na.rm = TRUE), 2))
}
message("within (", nrow(within), " pares): mediana rho ", round(median(within$rho), 3),
        " top20 ", round(median(within$top20, na.rm = TRUE), 2))
message("null (", nrow(null), " pares): mediana rho ", round(median(null$rho), 3),
        " top20 ", round(median(null$top20, na.rm = TRUE), 2))
