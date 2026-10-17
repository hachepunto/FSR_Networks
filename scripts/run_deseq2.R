suppressMessages({
  library(yaml)
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

source("scripts/lib_common.R")

# ============================================================
# EXPRESIÓN DIFERENCIAL (DESeq2) — tumor vs. control, por analysis_id
#
# Basado en un script de referencia de expresión diferencial con DESeq2
# (INMEGEN 2024); el flujo estadístico (diseño ~condition, filtrado por
# cuenta mínima, results()+contraste) no cambia. La transformación para
# las figuras de control de calidad sí: aquí se usa vst() en vez de
# rlog() (ver Sección 5) porque varios analysis_id tienen cientos de
# muestras, y la propia documentación de DESeq2 recomienda vst() sobre
# rlog() por velocidad pasadas ~50-100 muestras. Lo demás que sí cambia,
# para encajar con el resto de
# este pipeline:
#
#   - Cuentas crudas: se leen directamente del SummarizedExperiment ya
#     descargado (assay "unstranded"), no de un featureCounts.tsv +
#     GTF aparte — ya tenemos gene_name/gene_type en su rowData, igual
#     que preprocess_expression.R.
#   - Muestras: exactamente las mismas que ya fijó preprocess_expression.R
#     (mismas columnas de <analysis_id>_tumor_expr.rds y
#     control_expr.rds), no un archivo de metadatos aparte — así la
#     expresión diferencial y la matriz que ve la red parten de la
#     misma selección de muestras.
#   - Un solo analysis_id por corrida (recibido como argumento), para
#     que Snakemake lo corra en paralelo a la inferencia de redes: este
#     paso no depende de la red y la red no depende de él.
#   - Se quitan los pasos de la referencia que ya no aplican aquí:
#     import de GTF (rtracklayer), metadatos con más de 2 condiciones
#     (dds_all separado del subconjunto filtrado — aquí solo hay
#     tumor/control), y varios paquetes cargados mas no usados en el
#     script original (BUSpaRse, PCAtools, ComplexHeatmap — el propio
#     script llamaba pheatmap(), de otro paquete, no de ComplexHeatmap
#     — AnnotationDbi, enrichR, tidyr, tibble, purrr, magrittr, ggrepel,
#     dplyr).
#
# La columna que le interesa al resto del pipeline es `stat` — el
# estadístico de Wald de DESeq2 (combina log2FC y su significancia) —
# que reemplaza al t-test de msVIPER como firma de expresión de entrada
# para NaRnEA (instrucciones aparte, todavía no implementadas).
# ============================================================

# ============================================================
# ARGUMENTOS
#
# 1. analysis_id (p.ej. BRCA_basal) — el resto (cáncer, proyecto,
#    rutas) se deriva de config.yaml, igual que en el resto del
#    pipeline.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Uso: Rscript scripts/run_deseq2.R <analysis_id>")
}

analysis_id <- args[1]

config <- read_config()  # scripts/lib_common.R; respeta MRA_CONFIG

data_root   <- config$paths$data_root
work_root   <- config$paths$work_root
results_dir <- config$paths$results_dir

deseq2_cfg   <- config$differential_expression$deseq2
min_count    <- deseq2_cfg$min_count
min_samples  <- deseq2_cfg$min_samples
log2fc_th    <- deseq2_cfg$log2fc_threshold
padj_th      <- deseq2_cfg$padj_threshold

preproc_dir <- file.path(work_root, "preprocessed")
deseq2_dir  <- file.path(work_root, "deseq2")
tables_dir  <- file.path(results_dir, "final", "tables")
figures_dir <- file.path(results_dir, "final", "figures")

# ============================================================
# ENTRADAS
#
#   config/config.yaml
#   <data_root>/processed/<PROJECT>/<tumor>_rnaseq_se.rds   (cuentas crudas)
#   <work_root>/preprocessed/<CANCER>/<analysis_id>_tumor_expr.rds  (solo para las columnas)
#   <work_root>/preprocessed/<CANCER>/control_expr.rds              (solo para las columnas)
# ============================================================

# ============================================================
# SALIDAS
#
#   <work_root>/deseq2/<analysis_id>/dds.rds
#   <work_root>/deseq2/<analysis_id>/vst_matrix.rds
#   <work_root>/deseq2/<analysis_id>/stat_vector.rds   (gene -> stat; entrada de NaRnEA)
#   <work_root>/deseq2/<analysis_id>/session_info.log
#
#   results/final/tables/<analysis_id>_deseq2_results.tsv   (todos los genes probados)
#   results/final/tables/<analysis_id>_deseq2_deg.tsv       (subconjunto |log2FC|/padj, reporte)
#
#   results/final/figures/<analysis_id>_deseq2_pca.pdf
#   results/final/figures/<analysis_id>_deseq2_sample_correlation.pdf
#   results/final/figures/<analysis_id>_deseq2_volcano.pdf
# ============================================================

spec <- find_analysis_spec(analysis_id, config)

if (is.null(spec)) {
  stop("analysis_id desconocido en config.yaml: ", analysis_id)
}

cancer  <- spec$cancer
project <- spec$project
tumor   <- tumor_from_project(project)

deseq2_analysis_dir <- file.path(deseq2_dir, analysis_id)
dir.create(deseq2_analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

dds_file          <- file.path(deseq2_analysis_dir, "dds.rds")
vst_file          <- file.path(deseq2_analysis_dir, "vst_matrix.rds")
stat_file         <- file.path(deseq2_analysis_dir, "stat_vector.rds")
session_info_file <- file.path(deseq2_analysis_dir, "session_info.log")

results_table_file <- file.path(tables_dir, paste0(analysis_id, "_deseq2_results.tsv"))
deg_table_file      <- file.path(tables_dir, paste0(analysis_id, "_deseq2_deg.tsv"))

pca_figure_file         <- file.path(figures_dir, paste0(analysis_id, "_deseq2_pca.pdf"))
correlation_figure_file <- file.path(figures_dir, paste0(analysis_id, "_deseq2_sample_correlation.pdf"))
volcano_figure_file     <- file.path(figures_dir, paste0(analysis_id, "_deseq2_volcano.pdf"))

# ============================================================
# SECCIÓN 1: Muestras — las mismas que ya fijó preprocess_expression.R
# ============================================================

cancer_dir <- file.path(preproc_dir, cancer)
tumor_expr_file   <- file.path(cancer_dir, paste0(analysis_id, "_tumor_expr.rds"))
control_expr_file <- file.path(cancer_dir, "control_expr.rds")

if (!file.exists(tumor_expr_file) || !file.exists(control_expr_file)) {
  stop(
    analysis_id, ": faltan matrices preprocesadas (", tumor_expr_file, "). ",
    "Correr scripts/preprocess_expression.R primero."
  )
}

tumor_samples   <- colnames(readRDS(tumor_expr_file))
control_samples <- colnames(readRDS(control_expr_file))
target_samples  <- c(tumor_samples, control_samples)

message(analysis_id, ": ", length(tumor_samples), " muestras tumorales, ",
        length(control_samples), " de control (mismas que preprocess_expression.R)")

# ============================================================
# SECCIÓN 2: Cuentas crudas desde el SummarizedExperiment ya descargado
# ============================================================

se_file <- file.path(data_root, "processed", project, paste0(tumor, "_rnaseq_se.rds"))
se <- readRDS(se_file)

# El SE crudo puede tener varias columnas (alícuotas) con el mismo
# sample_submitter_id. Se aplica la misma deduplicación que
# preprocess_expression.R para quedarse con la misma columna que usó.
se <- cbind(
  dedup_one_per_patient(se[, se$sample_type == "Primary Tumor"]),
  dedup_one_per_patient(se[, se$sample_type == "Solid Tissue Normal"])
)

se_sel <- se[, se$sample_submitter_id %in% target_samples]

# No basta con comparar cuántas columnas coinciden: un sample_submitter_id
# duplicado en el SummarizedExperiment crudo (o uno de target_samples
# ausente ahí) podría dar el mismo ncol() por coincidencia y terminar
# asignando una muestra a la columna equivocada más abajo. Se exige una
# biyección real: ni duplicados ni conjuntos distintos.
if (anyDuplicated(se_sel$sample_submitter_id)) {
  stop(
    analysis_id, ": sample_submitter_id duplicado en el SummarizedExperiment crudo ",
    "para las muestras de este analysis_id"
  )
}
if (!setequal(se_sel$sample_submitter_id, target_samples)) {
  stop(
    analysis_id, ": las muestras del SummarizedExperiment crudo (", se_file, ") no ",
    "coinciden exactamente con las que ya fijó preprocess_expression.R — ¿se volvió a ",
    "descargar el proyecto después de correr ese script?"
  )
}

colnames(se_sel) <- se_sel$sample_submitter_id
se_sel <- se_sel[, target_samples]

rd <- as.data.frame(rowData(se_sel))
pc_idx <- rd$gene_type == "protein_coding"

counts_raw <- assay(se_sel[pc_idx, ], "unstranded")
rownames(counts_raw) <- rd$gene_name[pc_idx]

invalid_genes <- is.na(rownames(counts_raw)) | trimws(rownames(counts_raw)) == ""
counts_raw <- counts_raw[!invalid_genes, , drop = FALSE]

# Símbolos duplicados: SUMAR, no promediar (ver collapse_duplicate_rownames_sum
# en lib_common.R) — son cuentas, a diferencia del log2TPM que usa
# preprocess_expression.R para las redes.
counts_mat <- collapse_duplicate_rownames_sum(counts_raw)
storage.mode(counts_mat) <- "integer"

# ============================================================
# SECCIÓN 3: DESeq2
# ============================================================

condition <- factor(
  ifelse(colnames(counts_mat) %in% tumor_samples, "tumor", "control"),
  levels = c("control", "tumor")
)
col_data <- data.frame(condition = condition, row.names = colnames(counts_mat))

dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~condition)

keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]

message(analysis_id, ": ", nrow(dds), " genes tras el filtro (cuenta >= ", min_count,
        " en >= ", min_samples, " muestras)")

dds <- DESeq(dds)

# El orden de la comparación importa: tumor vs control, para que un
# estadístico positivo signifique "más alto en tumor" — misma
# convención que el t-test de msVIPER (rowTtest(x = tumor, y = normal))
# que este paso reemplaza.
res <- results(dds, contrast = c("condition", "tumor", "control"))
res_ordered <- res[order(res$pvalue), ]
res_df <- data.frame(gene = rownames(res_ordered), as.data.frame(res_ordered), row.names = NULL)

write.table(res_df, results_table_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(results_table_file, " guardado — ", nrow(res_df), " genes")

saveRDS(dds, dds_file)

# Firma de expresión para NaRnEA: gene -> stat, sin filtrar por
# significancia (NaRnEA, como msVIPER, necesita el ranking completo).
stat_vector <- setNames(res_df$stat, res_df$gene)
stat_vector <- stat_vector[!is.na(stat_vector)]
saveRDS(stat_vector, stat_file)
message(stat_file, " guardado — ", length(stat_vector), " genes con stat")

# ============================================================
# SECCIÓN 4: Subconjunto de genes diferencialmente expresados (reporte)
# ============================================================

res_subset <- subset(res, abs(res$log2FoldChange) >= log2fc_th)
deg <- subset(res_subset, !is.na(res_subset$padj) & res_subset$padj < padj_th)
deg_ordered <- deg[order(deg$pvalue), ]
deg_df <- data.frame(gene = rownames(deg_ordered), as.data.frame(deg_ordered), row.names = NULL)

write.table(deg_df, deg_table_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(deg_table_file, " guardado — ", nrow(deg_df),
        " genes (|log2FC| >= ", log2fc_th, ", padj < ", padj_th, ")")

# ============================================================
# SECCIÓN 5: Figuras de control de calidad
# ============================================================

# vst(), no rlog(): varios analysis_id de config.yaml superan por mucho
# las ~50-100 muestras a partir de las cuales la propia documentación de
# DESeq2 recomienda vst() en vez de rlog() por velocidad (p.ej.
# BRCA_luma, 562 tumor + 113 control) — con rlog() esas corridas hubieran
# sido drásticamente más lentas que la prueba de humo con STAD_ebv (n=63).
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)
saveRDS(vst_matrix, vst_file)

pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle(paste("PCA:", analysis_id)) +
  theme_bw(base_size = 12)

ggsave(pca_figure_file, pca_plot, width = 7, height = 6)
message(pca_figure_file, " guardado")

ordered_samples <- rownames(col_data)[order(col_data$condition)]
spearman_mat <- cor(vst_matrix[, ordered_samples, drop = FALSE], method = "spearman")
spearman_mat[upper.tri(spearman_mat)] <- NA

pdf(correlation_figure_file, width = 7, height = 6)
pheatmap(
  spearman_mat,
  main = paste("Correlación de Spearman:", analysis_id),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),
  fontsize = 6,
  na_col = "white",
  annotation_col = data.frame(
    condition = col_data[ordered_samples, "condition"],
    row.names = ordered_samples
  ),
  show_rownames = FALSE,
  show_colnames = FALSE,
  gaps_col = cumsum(table(col_data[ordered_samples, "condition"]))
)
dev.off()
message(correlation_figure_file, " guardado")

res_df$direction <- ifelse(
  !is.na(res_df$padj) & res_df$log2FoldChange > log2fc_th & res_df$padj < padj_th, "Activado",
  ifelse(!is.na(res_df$padj) & res_df$log2FoldChange < -log2fc_th & res_df$padj < padj_th,
         "Reprimido", "No significativo")
)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  geom_point(size = 1, na.rm = TRUE) +
  scale_color_manual(values = c(
    "Activado" = "#C0392B", "Reprimido" = "#2980B9", "No significativo" = "grey70"
  )) +
  geom_vline(xintercept = c(-log2fc_th, log2fc_th), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_th), linetype = "dashed", color = "grey50") +
  labs(
    title = paste("Volcano plot:", analysis_id),
    subtitle = paste0("tumor vs control — |log2FC| >= ", log2fc_th, ", padj < ", padj_th),
    x = "log2FoldChange", y = expression(-log[10]~"(padj)"), color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(volcano_figure_file, volcano_plot, width = 7, height = 6)
message(volcano_figure_file, " guardado")

# ============================================================
# SECCIÓN 6: Info de sesión
# ============================================================

sink(session_info_file)
message("Tabla de condiciones utilizada para el objeto dds de DESeq2")
print(col_data)
message("\nResumen del objeto results de DESeq2")
summary(res)
message("\nInformación de la sesión de R")
print(sessionInfo())
sink()

message("\n== ", analysis_id, ": expresión diferencial completa ==")
