# Datos sintéticos con red de regulación conocida, para medir lo que los
# datos reales no permiten (no hay red verdadera): precisión/recall de
# aristas por tipo de relación y recuperación de reguladores maestros.
#
# Modelo generativo (parámetros en network_benchmark.synthetic):
#   - TFs latentes z ~ Normal multivariada, en módulos con correlación
#     tf_module_cor dentro del módulo: crea dependencias indirectas
#     TF -> target que el DPI debe podar.
#   - Cada target tiene entre regulators_per_target[1] y [2] reguladores,
#     elegidos con pesos Gamma(1) (grado de salida heterogéneo).
#   - Cada arista es lineal, monótona no lineal (sigmoide con saturación)
#     o no monótona (en U, Spearman ~ 0), con signo +1 (p = 0.7) o -1 y
#     peso U(0.5, 1). Cada f() se estandariza a varianza 1 bajo N(0, 1).
#   - target latente = suma de efectos + ruido N(0, noise_sd[regime]); escala log2 con
#     media por gen y cuentas NegBin(size = 10) con factores de librería.
#   - Grupo "tumor" (n_network, en el que se infiere la red) y grupo
#     "reference" (n_reference). En "tumor", n_master_regulators TFs (con
#     regulón verdadero >= narnea.min_targets) se desplazan +-mr_shift_sd:
#     son los reguladores maestros verdaderos.
#   - Firma: estadístico de Wald de DESeq2, tumor vs reference (igual que
#     run_deseq2.R). Matriz de red: log2(CPM + 1) del grupo tumor.
#
# Uso: Rscript scripts/network_benchmark/05_simulate_synthetic.R [regime n rep]
#   sin argumentos genera todos (regímenes x n_network x réplicas).
# Salida por dataset: <out_dir>/synthetic/<regime>_n<N>_rep<R>/
#   truth.rds (aristas verdaderas, tipo, signo; reguladores maestros),
#   expr.rds, matrix.txt + tfs.txt (entrada de ARACNe3), signature.rds

suppressMessages(library(DESeq2))
source("scripts/network_benchmark/lib_benchmark.R")

config <- bench_config()
dirs <- bench_dirs(config)
sc <- config$network_benchmark$synthetic

# f() estandarizadas: media 0 y varianza 1 bajo z ~ N(0, 1)
standardize_f <- function(f) {
  zz <- qnorm(ppoints(20000))
  v <- f(zz)
  m <- mean(v); s <- sd(v)
  function(z) (f(z) - m) / s
}

make_edge_function <- function(type) {
  switch(type,
    linear = function(z) z,
    monotone_nonlinear = { c0 <- runif(1, -1, 1); standardize_f(function(z) plogis(3 * (z - c0))) },
    non_monotone = standardize_f(function(z) z^2))
}

simulate_dataset <- function(regime, n_network, rep_id) {
  set.seed(sc$seed + 1000 * n_network + rep_id)
  noise_sd <- sc$noise_sd[[regime]]
  n_tf <- sc$n_tfs; n_tg <- sc$n_targets; n_ref <- sc$n_reference
  tf_names <- sprintf("TF%03d", seq_len(n_tf))
  tg_names <- sprintf("G%04d", seq_len(n_tg))

  # Red verdadera
  out_w <- rgamma(n_tf, shape = 1)
  types <- names(sc$edge_types); type_p <- unlist(sc$edge_types)
  truth <- do.call(rbind, lapply(seq_len(n_tg), function(g) {
    m <- sample(sc$regulators_per_target[[1]]:sc$regulators_per_target[[2]], 1)
    regs <- sample(n_tf, m, prob = out_w)
    data.frame(tf = tf_names[regs], target = tg_names[g],
               type = sample(types, m, replace = TRUE, prob = type_p),
               sign = sample(c(1, -1), m, replace = TRUE, prob = c(0.7, 0.3)),
               weight = runif(m, 0.5, 1), stringsAsFactors = FALSE)
  }))
  fns <- lapply(truth$type, make_edge_function)

  # TFs latentes con módulos correlacionados
  n_tot <- n_network + n_ref
  module <- rep(seq_len(sc$n_tf_modules), length.out = n_tf)
  z <- matrix(0, n_tot, n_tf, dimnames = list(NULL, tf_names))
  for (k in seq_len(sc$n_tf_modules)) {
    idx <- which(module == k)
    sigma <- matrix(sc$tf_module_cor, length(idx), length(idx)); diag(sigma) <- 1
    z[, idx] <- matrix(rnorm(n_tot * length(idx)), n_tot) %*% chol(sigma)
  }

  # Reguladores maestros: desplazados en el grupo tumor (primeras n_network
  # filas). Se eligen entre TFs con regulón verdadero evaluable por NaRnEA.
  out_degree <- table(factor(truth$tf, levels = tf_names))
  eligible <- names(out_degree)[out_degree >= config$network_benchmark$narnea$min_targets]
  mr <- sample(eligible, sc$n_master_regulators)
  mr_sign <- sample(c(1, -1), length(mr), replace = TRUE)
  is_tumor <- seq_len(n_tot) <= n_network
  for (i in seq_along(mr)) z[is_tumor, mr[i]] <- z[is_tumor, mr[i]] + mr_sign[i] * sc$mr_shift_sd

  # Targets
  y <- matrix(rnorm(n_tot * n_tg, sd = noise_sd), n_tot, n_tg, dimnames = list(NULL, tg_names))
  for (e in seq_len(nrow(truth))) {
    y[, truth$target[e]] <- y[, truth$target[e]] +
      truth$sign[e] * truth$weight[e] * fns[[e]](z[, truth$tf[e]])
  }

  # Escala log2 -> cuentas NegBin
  latent <- cbind(z, y)
  latent <- scale(latent)
  mu_log2 <- runif(ncol(latent), 5, 11)
  log2_mean <- sweep(0.8 * latent, 2, mu_log2, "+")
  lib <- exp(rnorm(n_tot, 0, 0.2))
  counts <- matrix(rnbinom(length(log2_mean), size = 10, mu = 2^log2_mean * lib), n_tot)
  counts <- t(counts)
  dimnames(counts) <- list(colnames(latent), c(sprintf("T%03d", seq_len(n_network)), sprintf("R%03d", seq_len(n_ref))))

  # Firma: DESeq2 tumor vs reference
  cond <- factor(ifelse(is_tumor, "tumor", "reference"), levels = c("reference", "tumor"))
  dds <- DESeqDataSetFromMatrix(counts, data.frame(condition = cond), ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  stat <- results(dds, contrast = c("condition", "tumor", "reference"))$stat
  names(stat) <- rownames(counts)
  stat <- stat[!is.na(stat)]

  # Matriz de red: log2(CPM + 1) del grupo tumor
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  expr <- log2(cpm[, is_tumor] + 1)

  out <- file.path(dirs$synthetic, sprintf("%s_n%d_rep%d", regime, n_network, rep_id))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(edges = truth, master_regulators = data.frame(tf = mr, sign = mr_sign),
               tf_module = setNames(module, tf_names), regime = regime, noise_sd = noise_sd, n_network = n_network,
               replicate = rep_id), file.path(out, "truth.rds"))
  saveRDS(expr, file.path(out, "expr.rds"))
  saveRDS(stat, file.path(out, "signature.rds"))
  write.table(data.frame(gene = rownames(expr), expr, check.names = FALSE),
              file.path(out, "matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(tf_names, file.path(out, "tfs.txt"))
  message(sprintf("%s n=%d rep=%d: %d aristas verdaderas (%s), MRs: %s", regime, n_network, rep_id, nrow(truth),
                  paste(names(table(truth$type)), table(truth$type), sep = "=", collapse = " "),
                  paste(mr, collapse = ",")))
}

args <- commandArgs(trailingOnly = TRUE)
grid <- if (length(args) == 3) {
  data.frame(regime = args[1], n = as.integer(args[2]), rep = as.integer(args[3]))
} else {
  expand.grid(regime = names(sc$noise_sd), n = unlist(sc$n_network), rep = seq_len(sc$replicates),
              stringsAsFactors = FALSE)
}
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  if (file.exists(file.path(dirs$synthetic, sprintf("%s_n%d_rep%d", g$regime, g$n, g$rep), "signature.rds"))) next
  simulate_dataset(g$regime, g$n, g$rep)
}
