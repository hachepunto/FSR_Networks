library(SummarizedExperiment)
library(yaml)
library(CMScaller)
library(GCclassifier)

source("scripts/lib_common.R")

# ============================================================
# CONFIGURACIÓN Y RUTAS
# ============================================================

config <- read_config()  # scripts/lib_common.R; respeta MRA_CONFIG y lee UTF-8

data_root <- config$paths$data_root
results_dir <- config$paths$results_dir

proc_dir <- file.path(data_root, "processed")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

samples_file  <- file.path(results_dir, "subtype_assignments_by_sample.tsv")
summary_file  <- file.path(results_dir, "subtype_n_real.tsv")
control_file  <- file.path(results_dir, "control_n_real.tsv")

# ============================================================
# ASIGNACIÓN DE SUBTIPO POR MUESTRA — N REAL
#
# Aplica, por proyecto, el método de clasificación decidido y
# documentado en results/final/tables/subtype_classification_methods.tsv:
#   - tcga_annotated: reutiliza la columna candidata ya calculada
#     por TCGAquery_subtype(), cruzada con pacientes que sí tienen
#     expresión (STAR - Counts) y con una sola muestra por paciente.
#   - recomputed_classifier: recalcula el subtipo aplicando un
#     clasificador publicado (CMScaller para COAD, GCclassifier
#     para STAD) sobre la matriz de expresión real descargada.
#
# n_red (Sección 5 del plan) = número de muestras tumorales finales
# por subtipo, tras deduplicar a una muestra por paciente y excluir
# valores NA / indeterminados / no aprobados.
# ============================================================

# ============================================================
# ENTRADAS
#
#   <data_root>/processed/<PROJECT>/<tumor>_rnaseq_se.rds
#   <data_root>/processed/<PROJECT>/<tumor>_subtype.rds (cuando aplica)
# ============================================================

# ============================================================
# SALIDAS
#
#   results/subtype_assignments_by_sample.tsv  (una fila por muestra tumoral final)
#   results/subtype_n_real.tsv                 (n_red por cáncer y subtipo)
#   results/control_n_real.tsv                 (n de Solid Tissue Normal por cáncer)
# ============================================================

# Las decisiones de método por cáncer (classification_method,
# subtype_column, excluded_values, gcclassifier_scheme) se leen
# directamente de config$subtypes — NO se duplican aquí. La
# justificación completa de cada decisión vive en
# results/final/tables/subtype_classification_methods.tsv y
# results/final/tables/subtype_label_mapping.tsv.

# ── Utilidades comunes ──────────────────────────────────────

load_rnaseq <- function(project, tumor) {
  readRDS(file.path(proc_dir, project, paste0(tumor, "_rnaseq_se.rds")))
}

load_subtype_table <- function(project, tumor) {
  readRDS(file.path(proc_dir, project, paste0(tumor, "_subtype.rds")))
}

# ── Clasificación por proyecto ──────────────────────────────

classify_tcga_annotated <- function(tumor_se, project, tumor, spec) {
  subtype_table <- load_subtype_table(project, tumor)

  values <- as.character(subtype_table[[spec$subtype_column]])
  values[trimws(values) == ""] <- NA_character_
  names(values) <- subtype_table$patient

  values[tumor_se$patient]
}

# excluded_values aplica por igual a los 3 métodos de clasificación (no solo
# a tcga_annotated) — config.yaml lo declara como un campo general por
# cáncer, así que se aplica aquí una sola vez, después del dispatch, en vez
# de que cada función classify_* tenga que acordarse de manejarlo por su
# cuenta.
apply_excluded_values <- function(assigned, spec) {
  excluded_values <- unlist(spec$excluded_values)
  if (length(excluded_values) > 0) {
    assigned[!is.na(assigned) & assigned %in% excluded_values] <- NA_character_
  }
  assigned
}

classify_coad_cmscaller <- function(tumor_se) {
  counts <- assign_ensembl_rownames(assay(tumor_se, "unstranded"), rowData(tumor_se)$gene_id)

  pred <- CMScaller(
    emat = counts,
    rowNames = "ensg",
    RNAseq = TRUE,
    doPlot = FALSE,
    verbose = TRUE,
    seed = 42
  )

  as.character(pred$prediction)
}

classify_thca_brs <- function(tumor_se, project, tumor) {
  subtype_table <- load_subtype_table(project, tumor)
  result <- compute_thca_brs(tumor_se, subtype_table)  # scripts/lib_common.R

  assigned <- setNames(result$brs_class, result$patient)
  unname(assigned[tumor_se$patient])
}

classify_stad_gcclassifier <- function(tumor_se, scheme = "TCGA") {
  tpm <- assign_ensembl_rownames(assay(tumor_se, "tpm_unstrand"), rowData(tumor_se)$gene_id)
  tpm <- tpm[rowSums(is.na(tpm)) == 0, , drop = FALSE]

  log2tpm <- log2(tpm + 1)

  set.seed(42)  # classifyGC() usa impute.knn() internamente para genes de firma faltantes
  pred <- classifyGC(
    Expr = log2tpm,
    method = scheme,
    idType = "ENSEMBL",
    verbose = TRUE
  )

  assigned <- setNames(as.character(pred$subtype), pred$sample)
  unname(assigned[colnames(tumor_se)])
}

# ============================================================
# EJECUCIÓN POR PROYECTO
# ============================================================

sample_rows <- list()
summary_rows <- list()
control_rows <- list()

for (project in config$projects) {

  tumor <- tumor_from_project(project)
  cancer <- cancer_from_project(project)
  spec <- config$subtypes[[cancer]]

  if (is.null(spec)) {
    message(project, ": sin especificación de clasificación en config.yaml. Se omite.")
    next
  }

  message("\n== ", project, " (método: ", spec$classification_method, ") ==")

  se <- load_rnaseq(project, tumor)

  tumor_se  <- dedup_one_per_patient(se[, se$sample_type == "Primary Tumor"])
  normal_se <- dedup_one_per_patient(se[, se$sample_type == "Solid Tissue Normal"])

  message(project, ": ", ncol(tumor_se), " tumores únicos por paciente, ",
          ncol(normal_se), " controles STN únicos por paciente")

  control_rows[[project]] <- data.frame(
    cancer = cancer,
    project = project,
    n_control_real = ncol(normal_se),
    stringsAsFactors = FALSE
  )

  if (spec$classification_method == "tcga_annotated") {
    assigned <- classify_tcga_annotated(tumor_se, project, tumor, spec)
  } else if (spec$classification_method == "cmscaller") {
    assigned <- classify_coad_cmscaller(tumor_se)
  } else if (spec$classification_method == "gcclassifier") {
    assigned <- classify_stad_gcclassifier(tumor_se, scheme = spec$gcclassifier_scheme)
  } else if (spec$classification_method == "thca_brs_recomputed") {
    assigned <- classify_thca_brs(tumor_se, project, tumor)
  } else {
    stop("Método de clasificación desconocido: ", spec$classification_method)
  }

  assigned <- apply_excluded_values(assigned, spec)

  n_excluded_na <- sum(is.na(assigned))

  sample_rows[[project]] <- data.frame(
    cancer = cancer,
    project = project,
    patient = tumor_se$patient,
    sample_submitter_id = tumor_se$sample_submitter_id,
    classification_method = spec$classification_method,
    subtype_label = assigned,
    stringsAsFactors = FALSE
  )

  subtype_counts <- table(assigned, useNA = "no")

  summary_rows[[project]] <- data.frame(
    cancer = cancer,
    project = project,
    classification_method = spec$classification_method,
    subtype_label = names(subtype_counts),
    n_red = as.integer(subtype_counts),
    network_method = network_method_for_n(as.integer(subtype_counts)),
    stringsAsFactors = FALSE
  )

  message(project, ": ", n_excluded_na, " muestras excluidas por subtipo NA/indeterminado/no aprobado")
}

samples_df <- do.call(rbind, sample_rows)
summary_df <- do.call(rbind, summary_rows)
control_df <- do.call(rbind, control_rows)

write.table(samples_df, samples_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(summary_df, summary_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(control_df, control_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message("\n==== N REAL POR SUBTIPO (cáncer, subtipo, n_red, método de red) ====")
print(summary_df, row.names = FALSE)

message("\n==== N REAL DE CONTROL (Solid Tissue Normal, único por paciente) ====")
print(control_df, row.names = FALSE)

message("\nGuardado: ", samples_file)
message("Guardado: ", summary_file)
message("Guardado: ", control_file)
