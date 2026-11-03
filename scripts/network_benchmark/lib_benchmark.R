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

# Matriz tumoral del cáncer completo (preprocess_expression.R, <CANCER>_pooled).
read_pooled_expr <- function(config, cancer) {
  root <- config$network_benchmark$pooled$expr_root
  if (is.null(root)) root <- file.path(config$paths$work_root, "preprocessed")
  readRDS(file.path(root, cancer, paste0(cancer, "_pooled_tumor_expr.rds")))
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
                 nthreads = nthreads, verbose = FALSE, dpi = "target"),
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

# rank_aw(), to_narnea(), apply_narnea_overflow_patch() y run_narnea() viven
# en scripts/lib_common.R (las usa también scripts/run_narnea.R del pipeline).

# Concordancia entre dos resultados de NaRnEA sobre los reguladores comunes:
# Spearman de NES y fracción del top-k por |NES| compartida. Si hay menos
# de 2k reguladores comunes el top-k no discrimina y se devuelve NA.
# tops: barrido de k (propuesta de Daniel); el texto cita k = 20 y la figura de
# barrido muestra que las conclusiones no dependen de ese valor.
K_SWEEP <- c(5, 10, 15, 20, 30, 40, 50)
mr_concordance <- function(x, y, tops = K_SWEEP) {
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
