# Submuestreo de cohortes grandes reales: infiere la red del esquema A sobre
# n muestras elegidas al azar (sin reemplazo) de la matriz tumoral, con los
# mismos parámetros que el pipeline. La firma NO se recalcula: el efecto que
# se aísla es el de n sobre la RED. Una corrida = una (cohorte, n, réplica);
# idempotente.
#
# Uso: Rscript scripts/network_benchmark/12_subsample_real.R <analysis_id> <n> <replicate>
# Salida: <out_dir>/real/subsample/<analysis_id>_n<n>_rep<r>.rds (regulón top-k)
#         + _samples.txt (las muestras usadas)

suppressMessages(library(corto))
source("scripts/network_benchmark/lib_benchmark.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Uso: Rscript scripts/network_benchmark/12_subsample_real.R <analysis_id> <n> <replicate>")
analysis_id <- args[1]; n_sub <- as.integer(args[2]); rep_id <- as.integer(args[3])

config <- bench_config()
dirs <- bench_dirs(config)
nb <- config$network_benchmark
sub_dir <- file.path(dirs$out, "real", "subsample"); dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(sub_dir, sprintf("%s_n%d_rep%d.rds", analysis_id, n_sub, rep_id))
if (file.exists(out_file)) { message(basename(out_file), ": ya existe, se omite"); quit(save = "no") }

at <- analysis_table(config)
row <- at[at$analysis_id == analysis_id, ]
expr <- read_tumor_expr(config, analysis_id, row$cancer)
if (n_sub >= ncol(expr)) stop(analysis_id, ": n = ", n_sub, " no es menor que la cohorte (", ncol(expr), ")")
tfs <- intersect(read_tf_list(config), rownames(expr))

# Semilla determinista por (cohorte, n, réplica), independiente del orden de
# corrida. El índice de la cohorte entra multiplicado por 1e6 para que dos
# cohortes no compartan semilla (una suma de códigos de caracteres colisionaba).
set.seed(nb$subsample$seed + 1e6 * match(analysis_id, unlist(nb$subsample$cohorts)) + 1000 * n_sub + rep_id)
samples <- sort(sample(colnames(expr), n_sub))
sub <- expr[, samples, drop = FALSE]

# Mismo procedimiento que run_corto.R (esquema A)
t0 <- Sys.time()
reg_full <- run_corto_or_empty(rank_rows(sub), tfs, nb$corto$nbootstraps, r2p(nb$corto$r_min, n_sub), nb$threads)
reg <- topk_regulon(reg_full, nb$top_k_main)
attr(reg, "params") <- list(analysis_id = analysis_id, n = n_sub, replicate = rep_id,
                            minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(reg, out_file)
writeLines(samples, sub(".rds$", "_samples.txt", out_file))
message(sprintf("%s n=%d rep=%d: %d reguladores (%.1f min)", analysis_id, n_sub, rep_id, length(reg),
                attr(reg, "params")$minutes))
