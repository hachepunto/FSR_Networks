# corto sobre la matriz tumoral de un analysis_id, con dos configuraciones:
#
#   current — la regla anterior del pipeline: log2(TPM+1), p fijo
#             (network_benchmark.previous_rule.fixed_p; el |r| efectivo
#             depende de n).
#   lenient_r<r> — la propuesta (FSR) antes del recorte top-k: rangos
#             por gen (Pearson sobre rangos = Spearman) y |r| >= r, mismo
#             umbral de efecto para cualquier n; una red por cada r de
#             corto.r_min_grid. El top-k se aplica después, en la evaluación.
#
# Uso: Rscript scripts/network_benchmark/02_run_corto_real.R <analysis_id>
# Salida: <out_dir>/real/corto/<analysis_id>_{current,lenient_r<r>}.rds

suppressMessages(library(corto))
source("scripts/network_benchmark/lib_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/02_run_corto_real.R <analysis_id>")
analysis_id <- args[1]

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark

at <- analysis_table(config)
row <- at[at$analysis_id == analysis_id, ]
if (nrow(row) != 1) stop("analysis_id desconocido en config.yaml: ", analysis_id)

expr <- read_tumor_expr(config, analysis_id, row$cancer)
tfs <- intersect(read_tf_list(config), rownames(expr))
n <- ncol(expr)

ranks <- rank_rows(expr)
runs <- list(current = list(mat = expr, p = nb$previous_rule$fixed_p))
for (r in unlist(nb$corto$r_min_grid)) runs[[lenient_name(r)]] <- list(mat = ranks, p = r2p(r, n))

for (name in names(runs)) {
  out_file <- file.path(dirs$corto, paste0(analysis_id, "_", name, ".rds"))
  if (file.exists(out_file)) {
    message(analysis_id, " ", name, ": ya existe, se omite")
    next
  }
  t0 <- Sys.time()
  reg <- run_corto_or_empty(runs[[name]]$mat, tfs, nb$corto$nbootstraps, runs[[name]]$p, nb$threads)
  attr(reg, "params") <- list(analysis_id = analysis_id, config = name, n = n, p = runs[[name]]$p,
                              r_threshold = corto:::p2r(runs[[name]]$p, n),
                              n_genes = nrow(expr), n_tfs_in_matrix = length(tfs),
                              nbootstraps = nb$corto$nbootstraps,
                              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                              corto_version = as.character(packageVersion("corto")))
  saveRDS(reg, out_file)
  message(sprintf("%s %s: n=%d p=%.3g |r|>=%.3f -> %d reguladores (%.1f min)", analysis_id, name, n,
                  runs[[name]]$p, corto:::p2r(runs[[name]]$p, n), length(reg),
                  attr(reg, "params")$minutes))
}
