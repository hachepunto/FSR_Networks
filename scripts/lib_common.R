# Utilidades compartidas entre scripts del pipeline MRA multicáncer.
# Se cargan con source("scripts/lib_common.R") desde la raíz del repositorio.

# Una sola muestra tumoral o control por paciente: se ordena por
# sample_submitter_id (incluye vial/portion) y se conserva la primera
# por orden alfabético — regla simple, determinista y documentada.
# Usada tanto para clasificar el subtipo (classify_subtypes.R) como
# para preparar la expresión (preprocess_expression.R); debe ser la
# MISMA regla en ambos pasos o la muestra usada para etiquetar el
# subtipo y la usada para su expresión podrían divergir.
dedup_one_per_patient <- function(se) {
  ord <- order(se$sample_submitter_id)
  se <- se[, ord]
  se[, !duplicated(se$patient)]
}

# Lee config/config.yaml, o el archivo indicado en la variable de entorno
# MRA_CONFIG (p.ej. una copia con otras rutas de trabajo). Se lee con
# readLines + yaml.load porque read_yaml trunca el archivo en el primer
# carácter no ASCII cuando el locale no es UTF-8.
read_config <- function() {
  path <- Sys.getenv("MRA_CONFIG", "config/config.yaml")
  yaml::yaml.load(readLines(path, warn = FALSE, encoding = "UTF-8"))
}

# Deriva el nombre corto del tumor ("brca") y el código de cáncer
# ("BRCA") a partir de un project id de GDC ("TCGA-BRCA"). Centralizado
# aquí para que un cambio en la convención de nombres de proyecto no pueda
# desincronizarse entre los scripts que después necesitan coincidir
# exactamente (p.ej. al indexar config$subtypes[[cancer]]).
tumor_from_project <- function(project) tolower(sub("^TCGA-", "", project))
cancer_from_project <- function(project) toupper(tumor_from_project(project))

# Método de inferencia de red. Desde octubre de 2026 es corto para TODA
# cohorte, sin importar n (FSR; justificación en
# docs/supplementary_network_method.md). La regla anterior (ARACNe3 si
# n >= 95, corto si no) se conserva en ARACNE3_MIN_N solo como referencia
# histórica y para el módulo de benchmark; ningún script la usa para
# decidir. classify_subtypes.R y preprocess_expression.R siguen
# escribiendo la columna network_method por compatibilidad.
ARACNE3_MIN_N <- 95

network_method_for_n <- function(n) rep("corto", length(n))

# ── Cohortes chicas ──
# n_red < network_inference.min_n: red propia, resultados exploratorios
# (suplementario §3.6). Decisión de 2026-10-18.
network_min_n <- function(config) {
  m <- config$network_inference$min_n
  if (is.null(m)) 0L else as.integer(m)
}

analysis_tier <- function(n_red, config) ifelse(n_red >= network_min_n(config), "main", "exploratory")

# Matrices del cáncer completo (<CANCER>_pooled): todos los cánceres si
# network_inference.cancer_wide_network (red de sensibilidad, suplementario §3.8)
# o network_benchmark.pooled.write_matrices (paso 13 del benchmark).
pooled_cancers <- function(config) {
  if (isTRUE(config$network_inference$cancer_wide_network) ||
      isTRUE(config$network_benchmark$pooled$write_matrices)) names(config$subtypes) else character(0)
}

# Filtros de genes de preprocess_expression.R, aplicados a una matriz tumoral
# (y a su control, para que ambos conserven los mismos genes): sin NA, con
# varianza en el tumor, y log2TPM > 1 en >= 20 % de las muestras tumorales.
filter_tumor_genes <- function(tumor_expr, normal_expr) {
  complete <- rowSums(is.na(tumor_expr)) == 0 & rowSums(is.na(normal_expr)) == 0
  tumor_expr <- tumor_expr[complete, , drop = FALSE]
  zero_var <- matrixStats::rowSds(tumor_expr) == 0
  tumor_expr <- tumor_expr[!zero_var, , drop = FALSE]
  expressed <- rowSums(tumor_expr > 1) >= ceiling(0.2 * ncol(tumor_expr))
  tumor_expr <- tumor_expr[expressed, , drop = FALSE]
  list(tumor = tumor_expr, normal = normal_expr[rownames(tumor_expr), , drop = FALSE],
       n_na = sum(!complete), n_zero_var = sum(zero_var, na.rm = TRUE), n_low_expr = sum(!expressed))
}

# ── FSR: red con corto sobre rangos, umbral de efecto y top-k ──
# Compartidas por run_corto.R (pipeline) y scripts/network_benchmark/.

# |r| -> p bilateral para n muestras; inversa exacta de corto:::p2r, de
# modo que corto(p = r2p(r, n)) equivale a exigir |r| >= r sea cual sea n.
r2p <- function(r, n) {
  t <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t), df = n - 2)
}

# Rangos por gen (fila): la Pearson que calcula corto pasa a ser Spearman.
rank_rows <- function(m) {
  r <- t(apply(m, 1, rank))
  dimnames(r) <- dimnames(m)
  r
}

# Recorta cada regulón a sus k targets de mayor likelihood (fracción de
# bootstraps en que la arista ganó el DPI de corto). Empates por el campo
# opcional `tiebreak` (el módulo de benchmark lo llena con la MI en
# regulones de ARACNe3) o, si no existe, por |tfmode|. Iguala el tamaño de
# los regulones entre cohortes de n distinta. Única definición: la usa
# run_corto.R y scripts/network_benchmark/.
tiebreak_of <- function(x) if (is.null(x$tiebreak)) abs(unname(x$tfmode)) else unname(x$tiebreak)

topk_regulon <- function(reg, k) {
  out <- lapply(reg, function(x) {
    o <- order(x$likelihood, tiebreak_of(x), decreasing = TRUE)
    o <- o[seq_len(min(k, length(o)))]
    out <- list(tfmode = x$tfmode[o], likelihood = x$likelihood[o])
    if (!is.null(x$tiebreak)) out$tiebreak <- x$tiebreak[o]
    out
  })
  class(out) <- class(reg)
  out
}

# ── Regulón de corto -> NaRnEA ──
# Compartidas por run_narnea.R (pipeline) y scripts/network_benchmark/: el
# benchmark que justifica FSR corre NaRnEA exactamente así.

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

# Escribe una fila de manifiesto en formato TSV, agregando el
# encabezado solo si el archivo todavía no existe. Usado por
# download_data.R y preprocess_expression.R para no reimplementar el
# mismo patrón de apertura/append por separado.
append_manifest_row <- function(row, manifest_file) {
  write.table(
    row, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = !file.exists(manifest_file), append = file.exists(manifest_file)
  )
}

# Colapsa filas con el mismo rowname (p.ej. varios Ensembl ID mapeando al
# mismo símbolo HGNC) aplicando `fun` por columna dentro de cada grupo de
# filas duplicadas. No se usa directamente — ver los dos wrappers de
# abajo, uno por criterio de agregación.
collapse_duplicate_rownames <- function(mat, fun) {
  df <- as.data.frame(mat)
  df$.sym <- rownames(mat)
  collapsed <- do.call(rbind,
    lapply(split(df[, colnames(df) != ".sym", drop = FALSE], df$.sym),
           function(chunk) apply(chunk, 2, fun, na.rm = TRUE))
  )
  collapsed[, colnames(mat), drop = FALSE]
}

# Mediana por columna — para TPM/log2TPM, ya normalizado. Usada tanto por
# preprocess_expression.R (build_symbol_tpm, genes protein_coding) como por
# compute_thca_brs() más abajo — antes duplicada en ambos lugares.
collapse_duplicate_rownames_median <- function(mat) collapse_duplicate_rownames(mat, median)

# Suma por columna — el criterio correcto para CUENTAS crudas (dos
# Ensembl ID con el mismo símbolo son cuentas independientes que deben
# sumarse, no promediarse). Usada por run_deseq2.R al construir la
# matriz de cuentas para DESeq2.
collapse_duplicate_rownames_sum <- function(mat) collapse_duplicate_rownames(mat, sum)

# Reconstruye a qué cáncer/proyecto/label_spec corresponde un analysis_id
# ("BRCA_basal") sin asumir su estructura de texto (varios standardized_label
# ya tienen guión bajo, p.ej. "prox_inflam") — recorre config$subtypes y
# arma el mismo analysis_id que preprocess_expression.R y run_mra.R hasta
# encontrar el que coincide. Usada por scripts de línea de comandos por
# analysis_id (p.ej. run_deseq2.R) que reciben solo el analysis_id y
# necesitan derivar el resto (cáncer, proyecto, network_method) desde
# config.yaml, igual que ya lo hacen los scripts que recorren todos.
find_analysis_spec <- function(analysis_id, config) {
  for (cancer in names(config$subtypes)) {
    spec <- config$subtypes[[cancer]]
    for (original_value in names(spec$labels)) {
      label_spec <- spec$labels[[original_value]]
      if (identical(paste0(cancer, "_", label_spec$standardized_label), analysis_id)) {
        project <- config$projects[
          vapply(config$projects, function(p) cancer_from_project(p) == cancer, logical(1))
        ]
        return(list(
          cancer = cancer,
          project = project[1],
          original_value = original_value,
          label_spec = label_spec
        ))
      }
    }
  }
  NULL
}

# Reemplaza los rownames de una matriz por el gene_id de Ensembl sin
# versión (p.ej. "ENSG00000141510.16" -> "ENSG00000141510"), y descarta
# duplicados quedándose con la primera fila. Usada por classify_coad_cmscaller
# y classify_stad_gcclassifier en classify_subtypes.R.
assign_ensembl_rownames <- function(mat, gene_ids) {
  rownames(mat) <- sub("\\..*$", "", gene_ids)
  mat[!duplicated(rownames(mat)), , drop = FALSE]
}

# ============================================================
# THCA — BRAF-RAS Score (BRS), recalculado
#
# Agrawal et al. 2014, Cell (TCGA-THCA), Extended Experimental
# Procedures sección 14.1: BRS(t) = |v(t)-c(B)|^2 - |v(t)-c(R)|^2,
# distancia euclidiana normalizada (z-score) al centroide de 71
# genes en tumores BRAF-mutantes (B) vs RAS-mutantes (R).
#
# El paper solo publicó el score para 391/496 pacientes porque el
# exoma (necesario para DEFINIR B y R) solo estaba disponible para
# esos 391 — pero una vez fijados los centroides, clasificar una
# muestra nueva solo requiere su RNA-seq. Aquí no se reutilizan los
# centroides exactos del paper (no publicados como valores
# numéricos, solo como heatmap en la Figura S7A) — se recalculan
# usando como referencia los 391 pacientes cuyo BRAF_RAF_class ya
# publicó TCGAquery_subtype(), sobre nuestra propia expresión.
#
# La lista de 71 genes se transcribió manualmente de la Figura S7A
# (papers/NIHMS633017-supplement-7.pdf, panel A). 3 símbolos
# desactualizados resueltos por alias (org.Hs.eg.db) a su nombre
# HGNC actual; FLJ23867 no se encontró bajo ningún alias y se
# excluyó (70/71 genes finales, cobertura 100% en nuestra matriz).
# Ver results/final/tables/subtype_classification_methods.tsv para
# la referencia completa y la verificación de cobertura.
THCA_BRS_GENES_ORIGINAL <- c(
  "ANKRD46", "CYB561", "GNA14", "HGD", "KATNAL2", "KCNAB1", "KCNIP3", "LGI3", "MLEC", "NQO1",
  "SFTPC", "SLC4A4", "SORBS2", "ABTB2", "AHR", "ANKLE2", "ANXA1", "ANXA2P2", "ARNTL", "ASAP2",
  "BID", "CDC42EP1", "COL8A2", "CREB5", "CTSC", "CYP1B1", "DTX4", "DUSP5", "ETHE1", "FAM176A",
  "FAM20C", "FCHO1", "FLJ23867", "FN1", "FSTL3", "GABRB2", "GBP2", "ITGA3", "ITGB8", "KCNN4",
  "LAMB3", "LLGL1", "LY6E", "MDFIC", "MET", "PDE5A", "PDLIM4", "PLCD3", "PLEKHA6", "PNPLA5",
  "PPL", "PRICKLE1", "PROS1", "PTPRE", "PVRL4", "RASGEF1B", "RUNX1", "RUNX2", "SDC4", "SEL1L3",
  "SFTPB", "SLC35F2", "SOX4", "SPOCK2", "STK17B", "SYT12", "TACSTD2", "TBC1D2", "TGFBR1",
  "TM7SF4", "TMEM43"
)

THCA_BRS_ALIAS_TO_CURRENT <- c(FAM176A = "EVA1A", PVRL4 = "NECTIN4", TM7SF4 = "DCSTAMP")

THCA_BRS_GENES <- setdiff(
  unname(ifelse(
    THCA_BRS_GENES_ORIGINAL %in% names(THCA_BRS_ALIAS_TO_CURRENT),
    THCA_BRS_ALIAS_TO_CURRENT[THCA_BRS_GENES_ORIGINAL],
    THCA_BRS_GENES_ORIGINAL
  )),
  "FLJ23867"
)

# tumor_se: SummarizedExperiment ya deduplicado a Primary Tumor, una
# muestra por paciente (ver dedup_one_per_patient). subtype_table: la
# tabla de TCGAquery_subtype("thca").
#
# CORREGIDO (revisión externa, ver results/final/tables/subtype_classification_methods.tsv):
#   1. Grupos de referencia: se usa BRAFV600E_RAS (genotipo real por
#      mutación: BRAF_V600E / RAS / OTHER), NO BRAF_RAF_class. Esta
#      última es sign(BRAF_RAF_score) -- la salida del propio
#      clasificador -- así que ajustar los centroides con ella sería
#      circular. "OTHER" (drivers distintos de BRAF-V600E puro/RAS) se
#      excluye del ajuste, igual que en el paper.
#   2. Fórmula: BRS(t) = |v(t)-c(B)|_2 - |v(t)-c(R)|_2, la DIFERENCIA DE
#      DISTANCIAS EUCLIDIANAS (norma L2), no de distancias al cuadrado.
#      Confirmado contra el XML crudo del .docx original (Extended
#      Experimental Procedures 14.1): el "2" en |x-y|_2 está marcado
#      como subscript (notación de norma), no superscript.
#
# published_class (BRAF_RAF_class, "" tratada como NA) y reference_group
# (BRAFV600E_RAS para quienes definieron los centroides, NA para el
# resto) se devuelven por separado para que la validación pueda
# distinguir resustitución (pacientes que sí entraron al ajuste) de
# validación genuina (pacientes con published_class que NO entraron).
compute_thca_brs <- function(tumor_se, subtype_table) {
  rd <- as.data.frame(rowData(tumor_se))
  gene_idx <- rd$gene_name %in% THCA_BRS_GENES

  tpm <- assay(tumor_se[gene_idx, ], "tpm_unstrand")
  rownames(tpm) <- rd$gene_name[gene_idx]
  colnames(tpm) <- tumor_se$patient

  tpm_dedup <- collapse_duplicate_rownames_median(tpm)

  log2tpm <- log2(tpm_dedup + 1)

  genotype <- setNames(as.character(subtype_table$BRAFV600E_RAS), subtype_table$patient)
  genotype[trimws(genotype) == ""] <- NA_character_
  genotype <- genotype[colnames(log2tpm)]

  published_class <- setNames(as.character(subtype_table$BRAF_RAF_class), subtype_table$patient)
  published_class[trimws(published_class) == ""] <- NA_character_
  published_class <- published_class[colnames(log2tpm)]

  is_braf <- !is.na(genotype) & genotype == "BRAF_V600E"
  is_ras  <- !is.na(genotype) & genotype == "RAS"
  ref_samples <- is_braf | is_ras

  ref_mean <- rowMeans(log2tpm[, ref_samples, drop = FALSE])
  ref_sd   <- apply(log2tpm[, ref_samples, drop = FALSE], 1, sd)
  valid_genes <- ref_sd > 0 & !is.na(ref_sd)

  z <- sweep(log2tpm[valid_genes, , drop = FALSE], 1, ref_mean[valid_genes], "-")
  z <- sweep(z, 1, ref_sd[valid_genes], "/")

  centroid_braf <- rowMeans(z[, is_braf, drop = FALSE])
  centroid_ras  <- rowMeans(z[, is_ras, drop = FALSE])

  euclid_dist <- function(mat, centroid) sqrt(colSums((mat - centroid)^2))
  brs_score <- euclid_dist(z, centroid_braf) - euclid_dist(z, centroid_ras)

  data.frame(
    patient = colnames(z),
    reference_group = ifelse(ref_samples, genotype, NA_character_),
    published_class = published_class[colnames(z)],
    brs_score = brs_score,
    brs_class = ifelse(brs_score < 0, "Braf-like", "Ras-like"),
    stringsAsFactors = FALSE
  )
}
