# Utilidades compartidas del módulo network_benchmark (justificación del
# método de red). Se cargan con source() desde la raíz del repositorio.
#
# Config: read_config() de scripts/lib_common.R (config/config.yaml o MRA_CONFIG).

source("scripts/lib_common.R")

bench_config <- function() read_config()   # alias: el módulo lee el mismo config que el pipeline

# create = FALSE solo resuelve rutas (para leer sin efectos colaterales).
bench_dirs <- function(config, create = TRUE) {
  out <- config$network_benchmark$out_dir
  if (is.null(out)) out <- file.path(config$paths$work_root, "network_benchmark")
  dirs <- list(
    out        = out,
    reference  = file.path(out, "reference"),
    corto      = file.path(out, "real", "corto"),
    aracne3    = file.path(out, "real", "aracne3"),
    real_eval  = file.path(out, "real", "evaluation"),
    synthetic  = file.path(out, "synthetic"),
    tables     = file.path(config$paths$results_dir, "final", "tables"),
    figures    = file.path(config$paths$results_dir, "final", "figures")
  )
  if (create) for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# Un renglón por analysis_id, en el mismo orden que preprocess_expression.R.
# network_method es el de la REGLA ANTERIOR (ARACNe3 si n_red >=
# previous_rule.aracne3_min_n, corto si no): es lo que el benchmark compara,
# independientemente de lo que el pipeline use hoy (corto para todos).
analysis_table <- function(config) {
  min_n <- config$network_benchmark$previous_rule$aracne3_min_n
  rows <- list()
  for (project in config$projects) {
    cancer <- cancer_from_project(project)
    for (label in config$subtypes[[cancer]]$labels) {
      rows[[length(rows) + 1]] <- data.frame(
        analysis_id = paste0(cancer, "_", label$standardized_label), cancer = cancer,
        network_method = if (label$n_red >= min_n) "ARACNe3" else "corto",
        n_red = label$n_red, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

read_tf_list <- function(config) {
  f <- file.path(config$paths$work_root, "networks", "input", config$network_inference$tfs$file)
  tfs <- trimws(readLines(f))
  tfs[tfs != ""]
}

read_tumor_expr <- function(config, analysis_id, cancer) {
  readRDS(file.path(config$paths$work_root, "preprocessed", cancer,
                    paste0(analysis_id, "_tumor_expr.rds")))
}

# Firma de NaRnEA del pipeline: estadístico de Wald de DESeq2 (run_deseq2.R).
read_deseq2_stat <- function(config, analysis_id) {
  s <- readRDS(file.path(config$paths$work_root, "deseq2", analysis_id, "stat_vector.rds"))
  s[!is.na(s)]
}

# r2p(), rank_rows() y topk_regulon() viven en scripts/lib_common.R (las usa
# también run_corto.R del pipeline). Aquí solo lo específico del benchmark.

# corto aborta si ninguna arista pasa el umbral (pasa con el p fijo actual
# y n chica). Eso es un resultado, no un error: se devuelve un regulón
# vacío, que la evaluación reporta con 0 reguladores.
run_corto_or_empty <- function(mat, tfs, nbootstraps, p, nthreads) {
  tryCatch(
    corto::corto(inmat = mat, centroids = tfs, nbootstraps = nbootstraps, p = p,
                 nthreads = nthreads, verbose = FALSE),
    error = function(e) {
      if (!grepl("No edges passed", conditionMessage(e))) stop(e)
      message("corto: ninguna arista pasa p = ", signif(p, 3), "; regulón vacío")
      structure(list(), class = "regulon")
    })
}

# Nombre de la red laxa de corto para un umbral |r| dado (p.ej. "lenient_r0.4").
# Formatea a un decimal: r_min_grid debe usar valores con un decimal.
lenient_name <- function(r) sprintf("lenient_r%.1f", r)

# ── Regulones ──────────────────────────────────────────────
# Formato interno común (estilo corto/VIPER): lista por regulador con
# tfmode (modo con signo) y likelihood (confianza), nombrados por target.
# topk_regulon() y tiebreak_of() viven en scripts/lib_common.R.

# Red consolidada de ARACNe3 -> formato interno. tfmode = Spearman,
# likelihood = fracción de subredes con la arista, tiebreak = MI.
aracne3_to_regulon <- function(net, n_subnets, min_count = 1) {
  net <- net[net$count.values >= min_count, ]
  lapply(split(net, net$regulator.values), function(d) {
    list(tfmode = setNames(d$scc.values, d$target.values),
         likelihood = setNames(d$count.values / n_subnets, d$target.values),
         tiebreak = setNames(d$mi.values, d$target.values))
  })
}

read_aracne3_net <- function(dirs, analysis_id) {
  read.delim(file.path(dirs$aracne3, analysis_id, paste0("consolidated-net_", analysis_id, ".tsv")))
}

regulon_edges <- function(reg) {
  if (length(reg) == 0) return(data.frame(tf = character(0), target = character(0), mode = numeric(0),
                                           likelihood = numeric(0), stringsAsFactors = FALSE))
  data.frame(
    tf = rep(names(reg), vapply(reg, function(x) length(x$tfmode), integer(1))),
    target = unlist(lapply(reg, function(x) names(x$tfmode)), use.names = FALSE),
    mode = unlist(lapply(reg, function(x) unname(x$tfmode)), use.names = FALSE),
    likelihood = unlist(lapply(reg, function(x) unname(x$likelihood)), use.names = FALSE),
    stringsAsFactors = FALSE)
}

# aw por rango escalado a (0, 1] (receta de crc_mra / viñeta COAD de
# NaRnEA): el target con más soporte recibe 1, el último 1/n.
rank_aw <- function(likelihood, tiebreak) {
  n <- length(likelihood)
  o <- order(likelihood, tiebreak, decreasing = TRUE)
  aw <- numeric(n)
  aw[o] <- seq(n, 1) / n
  aw
}

to_narnea <- function(reg, aw = c("rank", "likelihood")) {
  aw <- match.arg(aw)
  lapply(reg, function(x) {
    w <- if (aw == "rank") rank_aw(x$likelihood, tiebreak_of(x)) else pmin(unname(x$likelihood), 1)
    list(aw = setNames(w, names(x$tfmode)), am = setNames(unname(x$tfmode), names(x$tfmode)))
  })
}

# ── NaRnEA ─────────────────────────────────────────────────

# NaRnEA::nes_covariance() calcula R * R con R entero (rangos): con más de
# 46 341 genes se desborda a NA. Se fuerza a double; nada más cambia. Con
# las firmas de este pipeline (~18 000 genes) no se dispara y el resultado
# es idéntico con y sin parche; se aplica por seguridad ante universos
# mayores (p.ej. incluir genes no codificantes).
apply_narnea_overflow_patch <- function() {
  fixed <- function(R, S, AW.AM.prod, AW.AM.abs.prod, E.r, E.rs, D.v, U.v, g) {
    storage.mode(R) <- "double"
    cov.prod.mat <- as.matrix(colSums(AW.AM.prod * AW.AM.abs.prod))
    E.r2s <- as.matrix((1 / g) * colSums(R * R * S))
    E.du <- cov.prod.mat %*% t(E.r2s)
    E.d.u <- cov.prod.mat %*% E.rs * E.r
    (E.du - E.d.u) / sqrt(D.v * U.v)
  }
  assignInNamespace("nes_covariance", fixed, ns = "NaRnEA")
}

# Una firma (vector nombrado) contra un regulón en formato NaRnEA. Targets
# fuera de la firma se descartan antes de aplicar min_targets.
run_narnea <- function(narnea_reg, ges, min_targets, seed) {
  narnea_reg <- lapply(narnea_reg, function(x) {
    keep <- names(x$aw) %in% names(ges)
    list(aw = x$aw[keep], am = x$am[keep])
  })
  narnea_reg <- narnea_reg[vapply(narnea_reg, function(x) length(x$aw), integer(1)) >= min_targets]
  if (length(narnea_reg) == 0) return(NULL)
  res <- NaRnEA::matrixNaRnEA(
    ges.mat = matrix(ges, ncol = 1, dimnames = list(names(ges), "signature")),
    regulon.list = narnea_reg, seed.val = seed, min.targets = min_targets)
  # Nombres desde rownames: con un solo regulador, [, "signature"] da un
  # escalar sin nombre.
  nes <- setNames(res$NES[, "signature"], rownames(res$NES))
  stopifnot("NaRnEA devolvió NES con NA" = !anyNA(nes))
  data.frame(tf = names(nes), nes = unname(nes), pes = unname(res$PES[rownames(res$NES), "signature"]),
             p = 2 * pnorm(-abs(unname(nes))), stringsAsFactors = FALSE)
}

# Concordancia entre dos resultados de NaRnEA sobre los reguladores comunes:
# Spearman de NES y fracción del top-k por |NES| compartida. Si hay menos
# de 2k reguladores comunes el top-k no discrimina y se devuelve NA.
mr_concordance <- function(x, y, tops = c(20, 50)) {
  if (is.null(x) || is.null(y)) return(c(n_shared = NA, rho = NA, setNames(rep(NA, length(tops)), paste0("top", tops))))
  m <- merge(x[, c("tf", "nes")], y[, c("tf", "nes")], by = "tf")
  top_overlap <- function(k) {
    if (nrow(m) < 2 * k) return(NA_real_)
    length(intersect(m$tf[order(-abs(m$nes.x))][1:k], m$tf[order(-abs(m$nes.y))][1:k])) / k
  }
  c(n_shared = nrow(m), rho = cor(m$nes.x, m$nes.y, method = "spearman"),
    setNames(vapply(tops, top_overlap, numeric(1)), paste0("top", tops)))
}

# ── CollecTRI (referencia externa) ─────────────────────────

read_collectri <- function(dirs) {
  ct <- read.delim(file.path(dirs$reference, "collectri.tsv"), stringsAsFactors = FALSE)
  ct <- ct[!grepl("^COMPLEX", ct$source) & !grepl("_", ct$source_genesymbol), ]
  ct$sign <- ifelse(ct$consensus_stimulation == "True" & ct$consensus_inhibition != "True", 1,
             ifelse(ct$consensus_inhibition == "True" & ct$consensus_stimulation != "True", -1, NA))
  ct <- ct[!duplicated(paste(ct$source_genesymbol, ct$target_genesymbol)), ]
  ct$key <- paste(ct$source_genesymbol, ct$target_genesymbol)
  ct
}

collectri_regulon <- function(ct) {
  lapply(split(ct, ct$source_genesymbol), function(d) {
    am <- ifelse(is.na(d$sign), 0, d$sign)
    list(aw = setNames(rep(1, nrow(d)), d$target_genesymbol), am = setNames(am, d$target_genesymbol))
  })
}

# Estructura de un regulón + concordancia con CollecTRI. El enriquecimiento
# compara la fracción de aristas en CollecTRI contra la esperada si cada
# regulador eligiera al azar el mismo número de targets del universo.
# signature_genes (opcional): universo de la firma; regulators_ge_min se
# cuenta tras intersectar los targets con él, como hace run_narnea().
regulon_structure <- function(reg, universe, ct, min_targets = 30, signature_genes = NULL) {
  e <- regulon_edges(reg)
  if (nrow(e) == 0) return(NULL)
  sizes <- table(e$tf)
  sizes_sig <- if (is.null(signature_genes)) sizes else table(e$tf[e$target %in% signature_genes])
  in_ct <- paste(e$tf, e$target) %in% ct$key
  ct_u <- ct[ct$target_genesymbol %in% universe, ]
  n_ct_tf <- as.numeric(table(factor(ct_u$source_genesymbol, levels = names(sizes))))
  expected <- sum(as.numeric(sizes) * n_ct_tf / length(universe)) / sum(sizes)
  s <- ct$sign[match(paste(e$tf, e$target)[in_ct], ct$key)]
  ok <- !is.na(s)
  data.frame(
    n_regulators = length(sizes), n_edges = nrow(e),
    median_targets = median(as.numeric(sizes)),
    regulators_ge_min = sum(sizes_sig >= min_targets),
    frac_negative = mean(e$mode < 0),
    median_likelihood = median(e$likelihood),
    collectri_edges = sum(in_ct),
    collectri_enrichment = mean(in_ct) / expected,
    collectri_sign_agreement = if (any(ok)) mean(sign(e$mode[in_ct][ok]) == s[ok]) else NA)
}
