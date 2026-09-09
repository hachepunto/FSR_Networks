library(TCGAbiolinks)
library(SummarizedExperiment)

# ============================================================
# PATHS
# ============================================================
tcga_raw_dir <- "data/raw/tcga"
metabric_raw_dir <- "data/raw/metabric"
tcga_proc_dir <- "data/processed/tcga"
metabric_proc_dir <- "data/processed/metabric"

tcga_subtypes_file <- file.path(tcga_proc_dir, "tcga_brca_subtypes.rds")
tcga_rnaseq_file <- file.path(tcga_proc_dir, "tcga_brca_rnaseq_se.rds")

# Crear directorios de salida si no existen
dir.create(tcga_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metabric_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tcga_proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metabric_proc_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# DESCARGA DE DATOS
# TCGA-BRCA + METABRIC
#
# Descarga datos clínicos, subtipos y expresión para el
# análisis de Master Regulators.
# ============================================================

# ============================================================
# SALIDAS
#
# TCGA:
#   data/processed/tcga/tcga_brca_subtypes.rds
#   data/processed/tcga/tcga_brca_rnaseq_se.rds
#
# METABRIC:
#   data/processed/metabric/mtbrc_clinical.rds
#   data/processed/metabric/mtbrc_clinical_mrna.rds
#   data/processed/metabric/mtbrc_expr.rds
# ============================================================

# Aumentar el tiempo máximo para descargas grandes
options(timeout = 1800)

# ============================================================
# SECCIÓN 1: TCGA-BRCA
# ============================================================

if (file.exists(tcga_subtypes_file)) {
    message("El archivo de subtipos TCGA ya existe. Se omite su descarga")
} else {
    message("Falta el archivo de subtipos TCGA. Iniciando descarga y procesamiento")
    query_clin <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
)

GDCdownload(query_clin, directory = tcga_raw_dir)
clinical_data <- GDCprepare(query_clin, directory = tcga_raw_dir)

names(clinical_data)

patient_table <- clinical_data[["clinical_patient_brca"]]

# ── 1.2  Subtipos PAM50 ─────────────────────────────────────

brca_subtypes <- TCGAquery_subtype(tumor = "BRCA")
table(brca_subtypes$BRCA_Subtype_PAM50, useNA = "always")


# ── 1.3  Cruzar PAM50 con estado TNBC por IHC ──────────────

patient_table$patient_id_short <- substr(patient_table$bcr_patient_barcode, 1, 12)

tcga_brca_subtypes <- merge(
  patient_table[, c("patient_id_short", "er_status_by_ihc",
                    "pr_status_by_ihc",  "her2_status_by_ihc")],
  brca_subtypes[, c("patient", "BRCA_Subtype_PAM50")],
  by.x = "patient_id_short",
  by.y = "patient"
)

tcga_brca_subtypes$is_TNBC <-
  tcga_brca_subtypes$er_status_by_ihc   %in% "Negative" &
  tcga_brca_subtypes$pr_status_by_ihc   %in% "Negative" &
  tcga_brca_subtypes$her2_status_by_ihc %in% "Negative"

table(tcga_brca_subtypes$BRCA_Subtype_PAM50, tcga_brca_subtypes$is_TNBC, useNA = "always")

saveRDS(
    tcga_brca_subtypes,
    tcga_subtypes_file
    )
}

message("Subtipos TCGA guardados:", tcga_subtypes_file)

# ── 1.4  Descarga RNA-seq STAR Counts ───────────────────────

if (file.exists(tcga_rnaseq_file)) {
    message("El archivo RNA-seq TCGA ya existe. Se omite su descarga")
} else {
    message(format(Sys.time()), " — Iniciando descarga TCGA-BRCA RNA-seq")

    query_rnaseq <- GDCquery(
      project       = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )

    message(format(Sys.time()), " — Query OK: ", nrow(getResults(query_rnaseq)), " archivos")

    GDCdownload(query_rnaseq, directory = tcga_raw_dir, files.per.chunk = 100)

    message(format(Sys.time()), " — Descarga completa. Preparando SummarizedExperiment...")

    tcga_brca_rnaseq_se <- GDCprepare(query_rnaseq, directory = tcga_raw_dir)

    message(format(Sys.time()), " — SE creado: ",
            nrow(tcga_brca_rnaseq_se), " genes × ", ncol(tcga_brca_rnaseq_se), " muestras")

    saveRDS(
        tcga_brca_rnaseq_se,
        tcga_rnaseq_file
    )
}


message("TCGA RNA-seq guardado: ", tcga_rnaseq_file)

# ============================================================
# SECCIÓN 2: METABRIC
# ============================================================

# ── 2.1  Descarga desde cBioPortal ──────────────────────────

if (!file.exists(file.path(metabric_raw_dir, "brca_metabric.tar.gz"))) {

    download.file(
        url = "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz",
        destfile = file.path(metabric_raw_dir, "brca_metabric.tar.gz"),
        mode = "wb"
    )

}

if (!dir.exists(file.path(metabric_raw_dir, "brca_metabric"))) {

    untar(
        file.path(metabric_raw_dir, "brca_metabric.tar.gz"),
        exdir = metabric_raw_dir
    )

}

message("METABRIC extracted in data/raw/metabric/brca_metabric/")

# ── 2.2  Clínica ─────────────────────────────────────────────

mtbrc_patient <- read.delim(
    file.path(metabric_raw_dir, "brca_metabric", "data_clinical_patient.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
)

mtbrc_sample <- read.delim(
    file.path(metabric_raw_dir, "brca_metabric", "data_clinical_sample.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "#"
)

mtbrc_clinical <- merge(
  mtbrc_patient,
  mtbrc_sample[, c("PATIENT_ID", "SAMPLE_ID", "ER_STATUS", "HER2_STATUS", "PR_STATUS")],
  by    = "PATIENT_ID",
  all.x = TRUE
)

table(mtbrc_clinical$CLAUDIN_SUBTYPE, useNA = "always")
table(mtbrc_clinical$INTCLUST,        useNA = "always")

saveRDS(
    mtbrc_clinical,
    file.path(metabric_proc_dir, "mtbrc_clinical.rds")
)

message(
    "Clinical data saved: ",
    file.path(metabric_proc_dir, "mtbrc_clinical.rds"),
    " (", nrow(mtbrc_clinical), " patients)"
)

# ── 2.3  Clínica restringida a muestras con expresión ────────

mrna_ids_raw <- readLines(
	file.path(metabric_raw_dir, "brca_metabric", "case_lists", "cases_mRNA.txt")
)
mrna_ids_line <- grep("^case_list_ids:", mrna_ids_raw, value = TRUE)
mrna_ids      <- strsplit(sub("^case_list_ids:\\s*", "", mrna_ids_line), "\t")[[1]]

mtbrc_clinical_mrna <- mtbrc_clinical[mtbrc_clinical$SAMPLE_ID %in% mrna_ids, ]

table(mtbrc_clinical_mrna$CLAUDIN_SUBTYPE, useNA = "always")

saveRDS(
    mtbrc_clinical_mrna,
    file.path(metabric_proc_dir, "mtbrc_clinical_mrna.rds")
)

message("Clinical data with expression saved: ", nrow(mtbrc_clinical_mrna), " samples")

# ── 2.4  Matriz de expresión (microarray Illumina) ───────────

message("Leyendo data_mrna_illumina_microarray.txt (~658 MB, ~1 min)...")

mtbrc_expr_raw <- read.delim(
    file.path(metabric_raw_dir, "brca_metabric", "data_mrna_illumina_microarray.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Convertir a matriz numérica (columnas 3 en adelante son muestras)
expr_mat <- as.matrix(mtbrc_expr_raw[, -(1:2)])
rownames(expr_mat) <- mtbrc_expr_raw$Hugo_Symbol

invalid_genes <- is.na(rownames(expr_mat)) | trimws(rownames(expr_mat)) == ""
expr_mat <- expr_mat[!invalid_genes, , drop = FALSE]

message(sum(invalid_genes), " genes sin símbolo HGNC eliminados")

# Resolver duplicados: colapsar por mediana de columnas
dup_genes <- unique(rownames(expr_mat)[duplicated(rownames(expr_mat))])
message(length(dup_genes), " símbolos duplicados — colapsando por mediana")

df_expr <- as.data.frame(expr_mat)
df_expr$.symbol <- rownames(expr_mat)
mtbrc_expr <- do.call(rbind,
  lapply(
    split(df_expr[, colnames(df_expr) != ".symbol"], df_expr$.symbol),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)

message("Dimensiones finales: ", nrow(mtbrc_expr), " genes × ", ncol(mtbrc_expr), " muestras")

saveRDS(
    mtbrc_expr,
    file.path(metabric_proc_dir, "mtbrc_expr.rds")
)

message("Expression matrix saved.")
message("\n========================================")
message("Download completed successfully")
message("========================================")

message("Processed files:")

message("  data/processed/tcga/tcga_brca_subtypes.rds")
message("  data/processed/tcga/tcga_brca_rnaseq_se.rds")

message("  data/processed/metabric/mtbrc_clinical.rds")
message("  data/processed/metabric/mtbrc_clinical_mrna.rds")
message("  data/processed/metabric/mtbrc_expr.rds")
