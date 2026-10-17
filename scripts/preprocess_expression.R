library(SummarizedExperiment)
library(yaml)
library(matrixStats)

source("scripts/lib_common.R")

# ============================================================
# CONFIGURACIÓN Y RUTAS
# ============================================================

config <- yaml::read_yaml("config/config.yaml")

data_root <- config$paths$data_root
work_root <- config$paths$work_root
results_dir <- config$paths$results_dir

proc_dir <- file.path(data_root, "processed")
preproc_dir <- file.path(work_root, "preprocessed")
network_input_dir <- file.path(work_root, "networks", "input")

assignments_file <- file.path(results_dir, "subtype_assignments_by_sample.tsv")
manifest_file <- file.path(results_dir, "preprocess_manifest.tsv")

# ============================================================
# PREPROCESAMIENTO POR analysis_id (cáncer × subtipo)
#
# Para cada proyecto (una vez):
#   - se construye la matriz TPM de genes protein_coding con símbolo
#     HGNC (mismo criterio que la versión BRCA-only original);
#   - se seleccionan los controles Solid Tissue Normal, una muestra
#     por paciente — se preparan una sola vez y se comparten entre
#     todos los subtipos de ese proyecto.
#
# Para cada analysis_id = <CANCER>_<standardized_label> (uno por
# subtipo con included=TRUE en config$subtypes):
#   - se seleccionan las muestras tumorales de ese subtipo (según
#     results/subtype_assignments_by_sample.tsv, ya congelado por
#     classify_subtypes.R — este script NO vuelve a clasificar nada);
#   - log2(TPM + 1);
#   - se eliminan genes con NA, sin varianza o sin expresión
#     suficiente (log2TPM > 1 en ≥20% de las muestras tumorales);
#   - tumor y control quedan con los mismos genes y el mismo orden;
#   - se guarda la matriz tumoral, la de control y la entrada ARACNe;
#   - se registra en el manifiesto cuántas muestras/genes se
#     incluyeron y excluyeron, y n_red final.
# ============================================================

# ============================================================
# ENTRADAS
#
#   config/config.yaml
#   results/subtype_assignments_by_sample.tsv
#   <data_root>/processed/<PROJECT>/<tumor>_rnaseq_se.rds
# ============================================================

# ============================================================
# SALIDAS
#
#   <work_root>/preprocessed/<CANCER>/control_expr.rds
#   <work_root>/preprocessed/<CANCER>/<analysis_id>_tumor_expr.rds
#   <work_root>/preprocessed/<CANCER>/<analysis_id>_normal_expr.rds
#   <work_root>/networks/input/<analysis_id>_matrix.txt
#   results/preprocess_manifest.tsv
# ============================================================

dir.create(preproc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_input_dir, recursive = TRUE, showWarnings = FALSE)

# A diferencia de download_data.R (que se puede reanudar y por eso
# acumula su manifiesto entre corridas), este script siempre reprocesa
# todos los proyectos configurados desde cero en cada ejecución — el
# manifiesto debe reiniciarse aquí o cada corrida duplicaría sus filas.
if (file.exists(manifest_file)) file.remove(manifest_file)

assignments <- read.delim(assignments_file, stringsAsFactors = FALSE, na.strings = "NA")

# ── Utilidades comunes ──────────────────────────────────────
# dedup_one_per_patient() se carga de scripts/lib_common.R — es la
# MISMA regla que usa classify_subtypes.R para elegir muestra por
# paciente; no debe reimplementarse por separado aquí.

# Mapa Ensembl ID -> símbolo HGNC, restringido a protein_coding,
# colapsando símbolos duplicados por mediana (igual que la versión
# BRCA-only original).
build_symbol_tpm <- function(se) {
  rd <- as.data.frame(rowData(se))
  pc_idx <- rd$gene_type == "protein_coding"

  se_pc <- se[pc_idx, ]
  rd_pc <- rd[pc_idx, ]

  tpm <- assay(se_pc, "tpm_unstrand")
  rownames(tpm) <- rd_pc$gene_name
  colnames(tpm) <- se_pc$sample_submitter_id

  invalid_genes <- is.na(rownames(tpm)) | trimws(rownames(tpm)) == ""
  tpm <- tpm[!invalid_genes, , drop = FALSE]

  dup_syms <- rownames(tpm)[duplicated(rownames(tpm))]

  tpm_dedup <- collapse_duplicate_rownames_median(tpm)

  list(
    tpm = tpm_dedup,
    n_genes_protein_coding = sum(pc_idx),
    n_genes_no_symbol = sum(invalid_genes),
    n_symbols_duplicated = length(unique(dup_syms))
  )
}

write_aracne_matrix <- function(mat, path) {
  df <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

failed_manifest_row <- function(cancer, project, analysis_id, original_value, reason) {
  data.frame(
    cancer = cancer, project = project, analysis_id = analysis_id,
    subtype_original_value = original_value,
    n_tumor_samples = 0L, n_control_samples = NA_integer_,
    n_genes_protein_coding = NA_integer_, n_genes_removed_no_symbol = NA_integer_,
    n_symbols_collapsed_duplicated = NA_integer_, n_genes_removed_na = NA_integer_,
    n_genes_removed_zero_variance = NA_integer_, n_genes_removed_low_expression = NA_integer_,
    n_genes_final = 0L, network_method = paste0("FAILED_", reason),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# EJECUCIÓN POR PROYECTO
# ============================================================

for (project in config$projects) {

  tumor <- tumor_from_project(project)
  cancer <- cancer_from_project(project)
  spec <- config$subtypes[[cancer]]

  if (is.null(spec)) {
    message(project, ": sin especificación de subtipos en config.yaml. Se omite.")
    next
  }

  message("\n== ", project, " ==")

  se <- readRDS(file.path(proc_dir, project, paste0(tumor, "_rnaseq_se.rds")))

  cancer_dir <- file.path(preproc_dir, cancer)
  dir.create(cancer_dir, recursive = TRUE, showWarnings = FALSE)

  tumor_se  <- dedup_one_per_patient(se[, se$sample_type == "Primary Tumor"])
  normal_se <- dedup_one_per_patient(se[, se$sample_type == "Solid Tissue Normal"])

  if (ncol(normal_se) == 0) {
    warning(project, ": sin muestras de control Solid Tissue Normal. Se omite el proyecto completo.",
            call. = FALSE)
    for (original_value in names(spec$labels)) {
      analysis_id <- paste0(cancer, "_", spec$labels[[original_value]]$standardized_label)
      append_manifest_row(
        failed_manifest_row(cancer, project, analysis_id, original_value, "no_control_samples"),
        manifest_file
      )
    }
    next
  }

  # Se construye la matriz de símbolos protein_coding UNA sola vez para
  # tumor+control combinados (misma rowData en ambos, ya que vienen del
  # mismo objeto se) y se separa después — evita repetir el colapso de
  # símbolos duplicados (el paso más lento) dos veces por proyecto.
  combined_se <- cbind(tumor_se, normal_se)
  combined_built <- build_symbol_tpm(combined_se)
  combined_log2 <- log2(combined_built$tpm + 1)

  normal_log2 <- combined_log2[, normal_se$sample_submitter_id, drop = FALSE]
  tumor_log2_all <- combined_log2[, tumor_se$sample_submitter_id, drop = FALSE]

  control_file <- file.path(cancer_dir, "control_expr.rds")
  saveRDS(normal_log2, control_file)
  message(project, ": control STN preparado — ", nrow(normal_log2), " genes × ",
          ncol(normal_log2), " muestras (", control_file, ")")

  project_assignments <- assignments[assignments$project == project & !is.na(assignments$subtype_label), ]

  for (original_value in names(spec$labels)) {

    label_spec <- spec$labels[[original_value]]
    standardized_label <- label_spec$standardized_label
    analysis_id <- paste0(cancer, "_", standardized_label)

    sample_ids <- project_assignments$sample_submitter_id[
      project_assignments$subtype_label == original_value
    ]
    sample_ids <- intersect(sample_ids, colnames(tumor_log2_all))

    # Se requieren al menos 2 muestras: con 1 sola, sd() no está definida y
    # matrixStats::rowSds() más abajo devolvería NA para todos los genes,
    # lo que insertaría una fila de NA en vez de filtrarla (indexar una
    # matriz con NA en un vector lógico agrega una fila NA, no la quita).
    if (length(sample_ids) < 2) {
      warning(analysis_id, ": ", length(sample_ids), " muestra(s) tumoral(es) asignada(s) — ",
              "insuficiente para calcular varianza. Se omite.", call. = FALSE)
      append_manifest_row(
        failed_manifest_row(cancer, project, analysis_id, original_value, "insufficient_tumor_samples"),
        manifest_file
      )
      next
    }

    tumor_expr <- tumor_log2_all[, sample_ids, drop = FALSE]
    normal_expr <- normal_log2

    # -- Genes con NA en tumor o control --
    complete_genes <- rowSums(is.na(tumor_expr)) == 0 & rowSums(is.na(normal_expr)) == 0
    n_na_removed <- sum(!complete_genes)
    tumor_expr <- tumor_expr[complete_genes, , drop = FALSE]
    normal_expr <- normal_expr[rownames(tumor_expr), , drop = FALSE]

    # -- Genes sin varianza en el tumor --
    zero_var <- matrixStats::rowSds(tumor_expr) == 0
    n_zero_var_removed <- sum(zero_var, na.rm = TRUE)
    tumor_expr <- tumor_expr[!zero_var, , drop = FALSE]
    normal_expr <- normal_expr[rownames(tumor_expr), , drop = FALSE]

    # -- Filtro de expresión: log2TPM > 1 en >= 20% de las muestras tumorales --
    expressed <- rowSums(tumor_expr > 1) >= ceiling(0.2 * ncol(tumor_expr))
    n_low_expr_removed <- sum(!expressed)
    tumor_expr <- tumor_expr[expressed, , drop = FALSE]
    normal_expr <- normal_expr[rownames(tumor_expr), , drop = FALSE]

    if (nrow(tumor_expr) == 0) {
      warning(analysis_id, ": ningún gen pasó los filtros. Se omite.", call. = FALSE)
      append_manifest_row(
        failed_manifest_row(cancer, project, analysis_id, original_value, "no_genes_passed_filters"),
        manifest_file
      )
      next
    }

    tumor_file  <- file.path(cancer_dir, paste0(analysis_id, "_tumor_expr.rds"))
    normal_file <- file.path(cancer_dir, paste0(analysis_id, "_normal_expr.rds"))
    aracne_file <- file.path(network_input_dir, paste0(analysis_id, "_matrix.txt"))

    saveRDS(tumor_expr, tumor_file)
    saveRDS(normal_expr, normal_file)
    write_aracne_matrix(tumor_expr, aracne_file)

    n_red_final <- ncol(tumor_expr)
    network_method <- network_method_for_n(n_red_final)

    if (n_red_final != label_spec$n_red) {
      warning(
        analysis_id, ": n_red final (", n_red_final,
        ") difiere del calculado en classify_subtypes.R (", label_spec$n_red, ") — ",
        "actualizar n_red en config.yaml. (El método de red ya no depende de n; es corto para todos.)",
        call. = FALSE
      )
    }

    append_manifest_row(data.frame(
      cancer = cancer,
      project = project,
      analysis_id = analysis_id,
      subtype_original_value = original_value,
      n_tumor_samples = n_red_final,
      n_control_samples = ncol(normal_expr),
      n_genes_protein_coding = combined_built$n_genes_protein_coding,
      n_genes_removed_no_symbol = combined_built$n_genes_no_symbol,
      n_symbols_collapsed_duplicated = combined_built$n_symbols_duplicated,
      n_genes_removed_na = n_na_removed,
      n_genes_removed_zero_variance = n_zero_var_removed,
      n_genes_removed_low_expression = n_low_expr_removed,
      n_genes_final = nrow(tumor_expr),
      network_method = network_method,
      stringsAsFactors = FALSE
    ), manifest_file)

    message(analysis_id, ": ", nrow(tumor_expr), " genes x ", n_red_final,
            " tumor / ", ncol(normal_expr), " control -> ", network_method)
  }
}

message("\nManifiesto: ", manifest_file)
