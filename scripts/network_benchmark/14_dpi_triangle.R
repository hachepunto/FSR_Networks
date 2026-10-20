# FSR con DPI por triángulos (regla de ARACNe, corto dpi = "triangle") en vez
# de la regla por blanco. Mismos parámetros que FSR (rangos, |r| >= r_min,
# 100 bootstraps); el recorte a top-k lo hacen 07 y 14b.
#
# Uso: Rscript scripts/network_benchmark/14_dpi_triangle.R synthetic <dataset_dir>
#      Rscript scripts/network_benchmark/14_dpi_triangle.R real <analysis_id>
# Salida: <dataset_dir>/corto_lenient_r<r>_triangle.rds
#         <out_dir>/real/corto_triangle/<analysis_id>_lenient_r<r>_triangle.rds

source("scripts/network_benchmark/lib_benchmark.R")
config <- bench_config()
nb <- config$network_benchmark
lib <- nb$dpi_triangle$corto_lib
if (!is.null(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(corto))
if (!"dpi" %in% names(formals(corto))) stop("Este corto (", packageVersion("corto"), ") no tiene el argumento dpi; ver network_benchmark.dpi_triangle.corto_lib")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || !(args[1] %in% c("synthetic", "real")))
  stop("Uso: 14_dpi_triangle.R synthetic <dataset_dir> | real <analysis_id>")
dirs <- bench_dirs(config)
stem <- paste0(lenient_name(nb$corto$r_min), "_triangle.rds")

if (args[1] == "synthetic") {
  ds <- args[2]
  out_file <- file.path(ds, paste0("corto_", stem))
  expr <- readRDS(file.path(ds, "expr.rds"))
  tfs <- intersect(readLines(file.path(ds, "tfs.txt")), rownames(expr))
} else {
  id <- args[2]
  out_dir <- file.path(dirs$out, "real", "corto_triangle"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(id, "_", stem))
  at <- analysis_table(config)
  expr <- read_tumor_expr(config, id, at$cancer[at$analysis_id == id])
  tfs <- intersect(read_tf_list(config), rownames(expr))
}
if (file.exists(out_file)) { message(basename(out_file), ": ya existe, se omite"); quit(save = "no") }

t0 <- Sys.time()
reg <- tryCatch(
  corto(rank_rows(expr), centroids = tfs, nbootstraps = nb$corto$nbootstraps,
        p = r2p(nb$corto$r_min, ncol(expr)), nthreads = nb$threads, dpi = "triangle"),
  error = function(e) { if (!grepl("No edges passed", conditionMessage(e))) stop(e); structure(list(), class = "regulon") })
attr(reg, "params") <- list(minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                            corto_version = as.character(packageVersion("corto")))
saveRDS(reg, out_file)
message(sprintf("%s: %d reguladores, %d aristas (%.1f min)", basename(out_file), length(reg),
                sum(vapply(reg, function(x) length(x$tfmode), integer(1))), attr(reg, "params")$minutes))
