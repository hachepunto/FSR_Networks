library(TCGAbiolinks)
library(SummarizedExperiment)
library(yaml)

source("scripts/lib_common.R")

# ============================================================
# CONFIGURACIÓN Y RUTAS
# ============================================================

config <- yaml::read_yaml("config/config.yaml")

data_root <- config$paths$data_root
projects <- config$projects

tcga_raw_dir  <- file.path(data_root, "raw")
tcga_proc_dir <- file.path(data_root, "processed")
tcga_log_dir  <- file.path(data_root, "logs")

# Aumentar el tiempo máximo para descargas grandes
options(timeout = 1800)

# ============================================================
# DESCARGA DE DATOS TCGA — proyectos configurados en config$projects
#
# Para cada proyecto de config$projects se descarga:
#   - expresión Gene Expression Quantification, STAR - Counts,
#     restringida a Primary Tumor + Solid Tissue Normal;
#   - clínica (Clinical Supplement, BCR Biotab), para vincular
#     barcodes con metadatos de paciente;
#   - tabla de subtipos vía TCGAquery_subtype(), cuando el
#     proyecto tenga una disponible.
#
# Un proyecto se omite por completo solo cuando sus tres
# salidas procesadas ya existen. No hay reanudación parcial
# dentro de un mismo paso.
# ============================================================

# ============================================================
# SALIDAS (por proyecto, p.ej. TCGA-BRCA -> tumor "brca")
#
#   <data_root>/processed/<PROJECT>/<tumor>_rnaseq_se.rds
#   <data_root>/processed/<PROJECT>/<tumor>_clinical.rds
#   <data_root>/processed/<PROJECT>/<tumor>_subtype.rds
#   <data_root>/logs/<tumor>_manifest.tsv
# ============================================================

download_project <- function(project) {

  tumor <- tumor_from_project(project)

  proj_raw_dir  <- file.path(tcga_raw_dir, project)
  proj_proc_dir <- file.path(tcga_proc_dir, project)

  rnaseq_file   <- file.path(proj_proc_dir, paste0(tumor, "_rnaseq_se.rds"))
  clinical_file <- file.path(proj_proc_dir, paste0(tumor, "_clinical.rds"))
  subtype_file  <- file.path(proj_proc_dir, paste0(tumor, "_subtype.rds"))
  manifest_file <- file.path(tcga_log_dir, paste0(tumor, "_manifest.tsv"))

  if (file.exists(rnaseq_file) && file.exists(clinical_file) && file.exists(subtype_file)) {
    message(project, ": las tres salidas procesadas ya existen. Se omite por completo.")
    return(invisible(NULL))
  }

  dir.create(proj_raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(proj_proc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tcga_log_dir, recursive = TRUE, showWarnings = FALSE)

  log_manifest_row <- function(step, data_category, data_type = NA_character_,
                                workflow_type = NA_character_, sample_type = NA_character_,
                                n_files = NA_integer_, n_rows = NA_integer_,
                                n_cols = NA_integer_, n_tables = NA_integer_,
                                output_file = NA_character_) {
    row <- data.frame(
      date = as.character(Sys.time()),
      project = project,
      step = step,
      data_category = data_category,
      data_type = data_type,
      workflow_type = workflow_type,
      sample_type = sample_type,
      n_files = n_files,
      n_rows = n_rows,
      n_cols = n_cols,
      n_tables = n_tables,
      tcgabiolinks_version = as.character(packageVersion("TCGAbiolinks")),
      output_file = output_file,
      stringsAsFactors = FALSE
    )

    # Se escribe de inmediato (no al final de download_project) para que un
    # error en un paso posterior no descarte el registro de un paso ya
    # completado y guardado en disco.
    append_manifest_row(row, manifest_file)
  }

  # ── 1. Expresión RNA-seq: Primary Tumor + Solid Tissue Normal ──

  if (file.exists(rnaseq_file)) {
    message(project, ": expresión ya existe. Se omite su descarga")
  } else {
    message(format(Sys.time()), " — ", project, ": iniciando descarga RNA-seq")

    rnaseq_ok <- tryCatch({
      query_rnaseq <- GDCquery(
        project       = project,
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        access        = "open",
        sample.type   = c("Primary Tumor", "Solid Tissue Normal")
      )

      n_files <- nrow(getResults(query_rnaseq))
      message(project, ": ", n_files, " archivos STAR - Counts encontrados (tumor + control)")

      GDCdownload(query_rnaseq, directory = proj_raw_dir, files.per.chunk = 100)

      rnaseq_se <- GDCprepare(query_rnaseq, directory = proj_raw_dir)

      message(
        project, ": SummarizedExperiment creado — ",
        nrow(rnaseq_se), " genes × ", ncol(rnaseq_se), " muestras — ",
        "tabla(sample_type): ",
        paste(capture.output(print(table(rnaseq_se$sample_type, useNA = "always"))), collapse = " | ")
      )

      saveRDS(rnaseq_se, rnaseq_file)

      log_manifest_row(
        step = "rnaseq",
        data_category = "Transcriptome Profiling",
        data_type = "Gene Expression Quantification",
        workflow_type = "STAR - Counts",
        sample_type = "Primary Tumor;Solid Tissue Normal",
        n_files = n_files,
        n_rows = nrow(rnaseq_se),
        n_cols = ncol(rnaseq_se),
        output_file = rnaseq_file
      )
      TRUE
    }, error = function(error) {
      message(project, ": ERROR al descargar RNA-seq — ", conditionMessage(error))
      log_manifest_row(
        step = "rnaseq",
        data_category = paste0("ERROR: ", conditionMessage(error))
      )
      FALSE
    })

    if (!rnaseq_ok) return(invisible(NULL))
  }

  # ── 2. Clínica (BCR Biotab), para vincular barcodes ──

  if (file.exists(clinical_file)) {
    message(project, ": clínica ya existe. Se omite su descarga")
  } else {
    message(format(Sys.time()), " — ", project, ": iniciando descarga clínica")

    tryCatch({
      query_clin <- GDCquery(
        project       = project,
        data.category = "Clinical",
        data.type     = "Clinical Supplement",
        data.format   = "BCR Biotab"
      )

      GDCdownload(query_clin, directory = proj_raw_dir)
      clinical_data <- GDCprepare(query_clin, directory = proj_raw_dir)

      saveRDS(clinical_data, clinical_file)

      message(project, ": clínica guardada — ", length(clinical_data), " tablas BCR Biotab")

      log_manifest_row(
        step = "clinical",
        data_category = "Clinical",
        data_type = "Clinical Supplement",
        workflow_type = "BCR Biotab",
        n_tables = length(clinical_data),
        output_file = clinical_file
      )
    }, error = function(error) {
      message(project, ": ERROR al descargar clínica — ", conditionMessage(error))
      log_manifest_row(
        step = "clinical",
        data_category = paste0("ERROR: ", conditionMessage(error))
      )
    })
  }

  # ── 3. Tabla de subtipos (TCGAquery_subtype), si existe ──

  if (file.exists(subtype_file)) {
    message(project, ": tabla de subtipos ya existe. Se omite su descarga")
  } else {
    message(format(Sys.time()), " — ", project, ": consultando tabla de subtipos")

    subtype_table <- tryCatch(
      TCGAquery_subtype(tumor = tumor),
      error = function(error) {
        message(project, ": sin tabla de subtipos disponible — ", conditionMessage(error))
        NULL
      }
    )

    if (!is.null(subtype_table)) {
      saveRDS(subtype_table, subtype_file)
      message(project, ": tabla de subtipos guardada — ", nrow(subtype_table), " registros")
    }

    log_manifest_row(
      step = "subtype",
      data_category = "TCGAquery_subtype",
      n_rows = if (!is.null(subtype_table)) nrow(subtype_table) else NA_integer_,
      output_file = if (!is.null(subtype_table)) subtype_file else NA_character_
    )
  }

  invisible(NULL)
}

for (project in projects) {
  tryCatch(
    download_project(project),
    error = function(error) {
      message(project, ": ERROR inesperado, se omite el resto de este proyecto — ", conditionMessage(error))
    }
  )
}

message("\n========================================")
message("Descarga TCGA completada para los ", length(projects), " proyectos configurados")
message("========================================")
