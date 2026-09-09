library(SummarizedExperiment)
library(sva)

# ============================================================
# PATHS
# ============================================================

tcga_proc_dir <- "data/processed/tcga"
metabric_proc_dir <- "data/processed/metabric"
aracne_input_dir <- "networks/input"

tcga_rnaseq_file <- file.path(tcga_proc_dir, "tcga_brca_rnaseq_se.rds")
mtbrc_clinical_file <- file.path(metabric_proc_dir, "mtbrc_clinical_mrna.rds")
mtbrc_expression_file <- file.path(metabric_proc_dir, "mtbrc_expr.rds")

subtype <- snakemake@wildcards[["subtype"]]

tcga_expr_file <- file.path(tcga_proc_dir, paste0("tcga_", subtype, "_expr.rds"))
tcga_normal_expr_file <- file.path(tcga_proc_dir, paste0("tcga_", subtype, "_normal_expr.rds"))
mtbrc_expr_file <- file.path(metabric_proc_dir, paste0("mtbrc_", subtype, "_expr.rds"))
mtbrc_normal_expr_file <- file.path(metabric_proc_dir, paste0("mtbrc_", subtype, "_normal_expr.rds"))
tcga_aracne_file <- file.path(aracne_input_dir, paste0("tcga_", subtype, "_matrix.txt"))
mtbrc_aracne_file <- file.path(aracne_input_dir, paste0("mtbrc_", subtype, "_matrix.txt"))

# Crear directorios de salida si no existen
dir.create(aracne_input_dir, recursive = TRUE, showWarnings = FALSE)
# ============================================================
# PRE-PROCESAMIENTO PARA ARACNE-AP
# TCGA-BRCA + METABRIC
#
# Genera matrices de expresión por subtipo molecular
# para construcción de redes regulatorias con ARACNe-AP
#
# Subtipos:
#   Basal
#   Luminal A
#   Luminal B
#   HER2
# ============================================================

# ============================================================
# ENTRADAS
#
# TCGA:
#   data/processed/tcga/tcga_brca_rnaseq_se.rds
#
# METABRIC:
#   data/processed/metabric/mtbrc_clinical_mrna.rds
#   data/processed/metabric/mtbrc_expr.rds
# ============================================================

# ============================================================
# SALIDAS
#
# matrices ARACNe:
#   networks/input/tcga_<subtype>_matrix.txt
#   networks/input/mtbrc_<subtype>_matrix.txt
#
# matrices expresión:
#   data/processed/tcga/tcga_<subtype>_expr.rds
#   data/processed/metabric/mtbrc_<subtype>_expr.rds
#
# controles normales:
#   data/processed/tcga/tcga_<subtype>_normal_expr.rds
#   data/processed/metabric/mtbrc_<subtype>_normal_expr.rds
# ============================================================

# ============================================================
# SECCIÓN 1: TCGA-BRCA
# ============================================================

# ── 1.1  Cargar SE y metadatos ──────────────────────────────

se <- readRDS(
    tcga_rnaseq_file
)

cd <- as.data.frame(colData(se))
rd <- as.data.frame(rowData(se))

# ── 1.2  Mapa Ensembl ID → símbolo HGNC (solo protein_coding) ─

pc_idx <- rd$gene_type == "protein_coding"
message("Genes protein_coding: ", sum(pc_idx))

# Ensembl IDs duplicados al mismo símbolo: conservar el de mayor expresión media
se_pc <- se[pc_idx, ]
rd_pc <- rd[pc_idx, ]

tpm_all <- assay(se_pc, "tpm_unstrand")
rownames(tpm_all) <- rd_pc$gene_name

invalid_genes <- is.na(rownames(tpm_all)) | trimws(rownames(tpm_all)) == ""

tpm_all <- tpm_all[!invalid_genes, , drop = FALSE]

message("Genes sin símbolo HGNC eliminados de TCGA: ", sum(invalid_genes))

dup_syms <- rownames(tpm_all)[duplicated(rownames(tpm_all))]
message("Símbolos duplicados (protein_coding): ", length(unique(dup_syms)))

df_tpm <- as.data.frame(tpm_all)
df_tpm$.sym <- rownames(tpm_all)
tpm_dedup <- do.call(rbind,
  lapply(
    split(df_tpm[, colnames(df_tpm) != ".sym"], df_tpm$.sym),
    function(chunk) apply(chunk, 2, median, na.rm = TRUE)
  )
)
message("Genes únicos tras colapso: ", nrow(tpm_dedup))

# ── 1.3  Asignar barcodes TCGA como colnames ────────────────

# Usar los primeros 16 caracteres del barcode (incluye tipo de muestra)
colnames(tpm_dedup) <- substr(cd$barcode, 1, 16)

# ── 1.4  Filtrar muestras por subtipo ───────────────────────

tcga_subtype_map <- c(Basal = "Basal", LumA = "LumA", LumB = "LumB", HER2 = "Her2")

if (!subtype %in% names(tcga_subtype_map)) {
    stop("Subtype not supported: ", subtype)
}

tcga_subtype <- tcga_subtype_map[[subtype]]

subtype_bc <- unique(substr(cd$barcode[
    cd$paper_BRCA_Subtype_PAM50 %in% tcga_subtype &
    cd$shortLetterCode %in% "TP"
], 1, 16))

normal_bc <- unique(substr(cd$barcode[
    cd$shortLetterCode %in% "NT"
], 1, 16))

if (length(subtype_bc) == 0) {
    stop("No TCGA tumor samples found for subtype: ", subtype)
}

if (length(normal_bc) == 0) {
    stop("No TCGA normal samples found")
}

cat("Muestras ", subtype, " TCGA:", length(subtype_bc), "\n")
cat("Muestras Normal TCGA:", length(normal_bc), "\n")
# ── 1.5  log2(TPM + 1) ───────────────────────────────────────────────────────

tcga_subtype_expr <- log2(tpm_dedup[, subtype_bc, drop = FALSE] + 1)
tcga_normal <- log2(tpm_dedup[, normal_bc, drop = FALSE] + 1)

# ── 1.6  Eliminar genes con NA ───────────────────────────────────────────────

complete_genes <- rowSums(is.na(tcga_subtype_expr)) == 0 & rowSums(is.na(tcga_normal)) == 0

tcga_subtype_expr <- tcga_subtype_expr[complete_genes, , drop = FALSE]
tcga_normal <- tcga_normal[complete_genes, , drop = FALSE]

message("Genes con NA eliminados de TCGA: ", sum(!complete_genes))

# ── 1.7  Filtro de expresión: ≥20 % del subtipo con log2TPM > 1 ──────────────

expressed <- rowSums(tcga_subtype_expr > 1) >= ceiling(0.2 * ncol(tcga_subtype_expr))

tcga_subtype_expr <- tcga_subtype_expr[expressed, , drop = FALSE]
tcga_normal <- tcga_normal[rownames(tcga_subtype_expr), , drop = FALSE]

if (nrow(tcga_subtype_expr) == 0) {
    stop("No genes passed the TCGA expression filter for subtype: ", subtype)
}

message("Genes tras filtro expresión TCGA: ", nrow(tcga_subtype_expr))

# ── 1.8  Guardar ────────────────────────────────────────────

saveRDS(tcga_subtype_expr, tcga_expr_file)
saveRDS(tcga_normal, tcga_normal_expr_file)

message(tcga_expr_file, " — ", nrow(tcga_subtype_expr), " genes × ", ncol(tcga_subtype_expr), " muestras")
message(tcga_normal_expr_file, " — ", nrow(tcga_normal), " genes × ", ncol(tcga_normal), " muestras")

# ── 1.9  Escribir matriz ARACNe-AP ──────────────────────────
# Formato: primera fila = "gene" \t sample1 \t sample2 ...

write_aracne_matrix <- function(mat, path) {
    df <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    message(path, " escrito — ", nrow(mat), " genes × ", ncol(mat), " muestras")
}

write_aracne_matrix(tcga_subtype_expr, tcga_aracne_file)

# ============================================================
# SECCIÓN 2: METABRIC
# ============================================================

# ── 2.1  Cargar datos ───────────────────────────────────────

cl <- readRDS(mtbrc_clinical_file)
expr <- readRDS(mtbrc_expression_file)

mtbrc_subtype <- subtype

if (!mtbrc_subtype %in% unique(cl$CLAUDIN_SUBTYPE)) {
    stop("No METABRIC samples found for subtype: ", subtype)
}

subtype_ids <- unique(cl$SAMPLE_ID[
    cl$CLAUDIN_SUBTYPE == mtbrc_subtype
])

normal_ids <- unique(cl$SAMPLE_ID[
    cl$CLAUDIN_SUBTYPE == "Normal"
])

if (length(subtype_ids) == 0) {
    stop("No METABRIC tumor samples found for subtype: ", subtype)
}

if (length(normal_ids) == 0) {
    stop("No METABRIC normal samples found")
}

cat("Muestras ", subtype, " METABRIC:", length(subtype_ids), "\n")
cat("Muestras Normal METABRIC:", length(normal_ids), "\n")

# ── 2.2  Subconjunto de muestras ────────────────────────────

mtbrc_subtype_expr <- expr[, colnames(expr) %in% subtype_ids, drop = FALSE]
mtbrc_normal <- expr[, colnames(expr) %in% normal_ids, drop = FALSE]

# ── 2.3  Eliminar genes con NA ──────────────────────────────

complete_genes <- rowSums(is.na(mtbrc_subtype_expr)) == 0 & rowSums(is.na(mtbrc_normal)) == 0

mtbrc_subtype_expr <- mtbrc_subtype_expr[complete_genes, , drop = FALSE]
mtbrc_normal <- mtbrc_normal[complete_genes, , drop = FALSE]

message("Genes con NA eliminados de METABRIC: ", sum(!complete_genes))

# ── 2.4  Filtro de expresión: señal log2 > 6 en ≥20 % del subtipo ───────────

expressed <- rowSums(mtbrc_subtype_expr > 6) >= ceiling(0.2 * ncol(mtbrc_subtype_expr))

mtbrc_subtype_expr <- mtbrc_subtype_expr[expressed, , drop = FALSE]
mtbrc_normal <- mtbrc_normal[rownames(mtbrc_subtype_expr), , drop = FALSE]

if (nrow(mtbrc_subtype_expr) == 0) {
    stop("No genes passed the METABRIC expression filter for subtype: ", subtype)
}

message("Genes tras filtro expresión METABRIC: ", nrow(mtbrc_subtype_expr))

# ── 2.5  ComBat: corrección de efecto de cohorte ─────────────
# METABRIC tiene 5 cohortes. ComBat preserva la señal biológica
# del subtipo versus Normal mientras elimina el efecto de batch.

combined <- cbind(mtbrc_subtype_expr, mtbrc_normal)

group_labels <- c(
    rep(subtype, ncol(mtbrc_subtype_expr)),
    rep("Normal", ncol(mtbrc_normal))
)

cohort_map <- setNames(cl$COHORT, cl$SAMPLE_ID)
batch_vec <- unname(cohort_map[colnames(combined)])

if (any(is.na(batch_vec))) {
    stop("Missing COHORT values found for METABRIC samples")
}

if (length(unique(batch_vec)) < 2) {
    stop("ComBat requires at least two METABRIC cohorts")
}

if (any(table(batch_vec) < 2)) {
    stop("Each METABRIC cohort must contain at least two samples for ComBat")
}

cat("Distribución de cohortes en ", subtype, " + Normal:\n")
print(table(batch_vec, group_labels))

mod <- model.matrix(~ group_labels)

message("Aplicando ComBat (batch = COHORT, mod = ~grupo)...")

combined_corrected <- ComBat(
    dat = combined,
    batch = batch_vec,
    mod = mod
)

mtbrc_subtype_expr <- combined_corrected[, group_labels == subtype, drop = FALSE]
mtbrc_normal <- combined_corrected[, group_labels == "Normal", drop = FALSE]

message("ComBat aplicado")

# ── 2.6  Guardar ────────────────────────────────────────────

saveRDS(mtbrc_subtype_expr, mtbrc_expr_file)
saveRDS(mtbrc_normal, mtbrc_normal_expr_file)

message(mtbrc_expr_file, " — ", nrow(mtbrc_subtype_expr), " genes × ", ncol(mtbrc_subtype_expr), " muestras")
message(mtbrc_normal_expr_file, " — ", nrow(mtbrc_normal), " genes × ", ncol(mtbrc_normal), " muestras")

write_aracne_matrix(mtbrc_subtype_expr, mtbrc_aracne_file)
