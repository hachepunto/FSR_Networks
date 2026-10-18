# Descarga CollecTRI (interacciones TF -> target curadas, vía la API de
# OmniPath) como referencia externa del módulo. Se registra la fecha de
# descarga y el número de interacciones, porque la API sirve la versión
# vigente, no una fija.
#
# Salida: <out_dir>/reference/collectri.tsv y collectri_provenance.tsv

source("scripts/network_benchmark/lib_benchmark.R")

config <- bench_config()
dirs <- bench_dirs(config)

out_file <- file.path(dirs$reference, "collectri.tsv")

archived <- file.path(dirs$tables, "network_benchmark_collectri_edges.tsv")
if (file.exists(out_file)) {
  message("CollecTRI ya existe: ", out_file, ". Se omite la descarga.")
} else if (file.exists(archived)) {
  # La API sirve la versión vigente; para reproducir la corrida publicada se
  # reutiliza la copia archivada en results/final/tables.
  file.copy(archived, out_file); file.copy(sub("_edges", "_provenance", archived), file.path(dirs$reference, "collectri_provenance.tsv"))
  message("CollecTRI copiado del archivo versionado: ", archived)
} else {
  download.file(config$network_benchmark$collectri_url, out_file, mode = "wb", quiet = TRUE)
  ct <- read.delim(out_file)
  stopifnot(c("source_genesymbol", "target_genesymbol", "consensus_stimulation",
              "consensus_inhibition") %in% colnames(ct))
  write.table(data.frame(source = config$network_benchmark$collectri_url,
                         downloaded = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                         n_interactions = nrow(ct)),
              file.path(dirs$reference, "collectri_provenance.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("CollecTRI: ", nrow(ct), " interacciones -> ", out_file)
}
