# Firma subtipo vs resto: DESeq2 de los tumores de un subtipo contra los de
# los demás subtipos del mismo cáncer (los controles no entran). Mismas
# muestras, cuentas, filtros y estadístico (Wald, columna stat) que
# run_deseq2.R, que hace tumor vs control.
#
# Uso: Rscript scripts/network_benchmark/15_subtype_vs_rest.R <analysis_id>
# Salida: <out_dir>/real/vs_rest/<analysis_id>_stat.rds (gene -> stat)

suppressMessages({
  library(SummarizedExperiment)
  library(DESeq2)
})
source("scripts/network_benchmark/lib_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/15_subtype_vs_rest.R <analysis_id>")
analysis_id <- args[1]

config <- bench_config()
dirs <- bench_dirs(config)
out_dir <- file.path(dirs$out, "real", "vs_rest"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(analysis_id, "_stat.rds"))
if (file.exists(out_file)) { message(basename(out_file), ": ya existe, se omite"); quit(save = "no") }

spec <- find_analysis_spec(analysis_id, config)
at <- analysis_table(config)
cancer <- spec$cancer
others <- setdiff(at$analysis_id[at$cancer == cancer], analysis_id)
if (length(others) == 0) stop(analysis_id, ": su cáncer no tiene otros subtipos")

test_samples <- colnames(read_tumor_expr(config, analysis_id, cancer))
rest_samples <- unique(unlist(lapply(others, function(o) colnames(read_tumor_expr(config, o, cancer)))))
stopifnot("Una muestra en dos subtipos" = !length(intersect(test_samples, rest_samples)))

# Cuentas crudas, misma deduplicación y colapso de símbolos que run_deseq2.R
tumor <- tumor_from_project(spec$project)
se <- readRDS(file.path(config$paths$data_root, "processed", spec$project, paste0(tumor, "_rnaseq_se.rds")))
se <- dedup_one_per_patient(se[, se$sample_type == "Primary Tumor"])
se <- se[, se$sample_submitter_id %in% c(test_samples, rest_samples)]
stopifnot(!anyDuplicated(se$sample_submitter_id), setequal(se$sample_submitter_id, c(test_samples, rest_samples)))
colnames(se) <- se$sample_submitter_id
rd <- as.data.frame(rowData(se)); pc <- rd$gene_type == "protein_coding"
counts <- assay(se[pc, ], "unstranded"); rownames(counts) <- rd$gene_name[pc]
counts <- counts[!(is.na(rownames(counts)) | trimws(rownames(counts)) == ""), , drop = FALSE]
counts <- collapse_duplicate_rownames_sum(counts); storage.mode(counts) <- "integer"

group <- factor(ifelse(colnames(counts) %in% test_samples, "subtype", "rest"), levels = c("rest", "subtype"))
dds <- DESeqDataSetFromMatrix(counts, data.frame(group = group, row.names = colnames(counts)), ~group)
dcfg <- config$differential_expression$deseq2
dds <- dds[rowSums(counts(dds) >= dcfg$min_count) >= dcfg$min_samples, ]
t0 <- Sys.time()
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "subtype", "rest"))
stat <- setNames(res$stat, rownames(res)); stat <- stat[!is.na(stat)]
attr(stat, "params") <- list(n_subtype = length(test_samples), n_rest = length(rest_samples),
                             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(stat, out_file)
message(sprintf("%s vs resto (%d vs %d tumores): %d genes (%.1f min)", analysis_id, length(test_samples),
                length(rest_samples), length(stat), attr(stat, "params")$minutes))
