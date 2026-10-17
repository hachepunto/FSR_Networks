# Convierte appendix_numbers.tsv (10_appendix_statistics.R) y
# network_benchmark_key_numbers.tsv (08_summarize.R) en macros LaTeX
# (\newcommand) para el apéndice técnico, de modo que ningún número se
# transcriba a mano. Los ids se convierten a nombres de macro válidos
# (solo letras): dígitos a palabras, el resto se elimina.
#
# Uso: Rscript scripts/network_benchmark/11_appendix_numbers_tex.R <appendix_dir>
# Salida: <appendix_dir>/numbers.tex

source("scripts/network_benchmark/lib_benchmark.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Uso: Rscript scripts/network_benchmark/11_appendix_numbers_tex.R <appendix_dir>")
out <- args[1]
config <- bench_config(); dirs <- bench_dirs(config, create = FALSE)

key <- read.delim(file.path(dirs$tables, "network_benchmark_key_numbers.tsv"), stringsAsFactors = FALSE)
app <- read.delim(file.path(out, "appendix_numbers.tsv"), stringsAsFactors = FALSE)
all <- rbind(data.frame(id = paste0("K_", key$id), value = key$value), data.frame(id = app$id, value = app$value))

digits_to_words <- c("0" = "Zero", "1" = "One", "2" = "Two", "3" = "Three", "4" = "Four",
                     "5" = "Five", "6" = "Six", "7" = "Seven", "8" = "Eight", "9" = "Nine")
macro_name <- function(id) {
  s <- gsub("\\.", "dot", id)
  for (d in names(digits_to_words)) s <- gsub(d, digits_to_words[[d]], s, fixed = TRUE)
  gsub("[^A-Za-z]", "", s)
}
fmt <- function(v) {
  if (is.na(v)) return("NA")
  if (v == round(v) && abs(v) < 1e7) return(format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
  format(signif(v, 3), scientific = FALSE, trim = TRUE)
}
names_tex <- macro_name(all$id)
stopifnot("nombres de macro duplicados" = !any(duplicated(names_tex)))
lines <- c("% Generado por scripts/network_benchmark/11_appendix_numbers_tex.R. No editar a mano.",
           sprintf("\\newcommand{\\%s}{%s}", names_tex, vapply(all$value, fmt, character(1))))
writeLines(lines, file.path(out, "numbers.tex"))
message(length(names_tex), " macros escritas en ", file.path(out, "numbers.tex"))
