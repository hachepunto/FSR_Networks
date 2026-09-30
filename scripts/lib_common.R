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

# Deriva el nombre corto del tumor ("brca") y el código de cáncer
# ("BRCA") a partir de un project id de GDC ("TCGA-BRCA"). Se usa en
# download_data.R, classify_subtypes.R, preprocess_expression.R,
# run_mra.R y audit_tcga_cohorts.R — centralizado aquí para que un cambio en la
# convención de nombres de proyecto no pueda desincronizarse entre
# scripts que después necesitan coincidir exactamente (p.ej. al
# indexar config$subtypes[[cancer]]).
tumor_from_project <- function(project) tolower(sub("^TCGA-", "", project))
cancer_from_project <- function(project) toupper(tumor_from_project(project))

# Umbral de N real (Sección 5 del plan) que decide el método de
# inferencia de red: por debajo de este valor, ARACNe3 pierde
# estabilidad en el subsampling/bootstrap y se usa corto en su lugar.
# Centralizado aquí porque classify_subtypes.R lo usa para calcular
# network_method y preprocess_expression.R lo vuelve a evaluar sobre
# el n_red final tras el filtrado de genes — ambos deben coincidir.
ARACNE3_MIN_N <- 95

network_method_for_n <- function(n) ifelse(n < ARACNE3_MIN_N, "corto", "ARACNe3")

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
