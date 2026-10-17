# Rellena docs/supplementary_network_method.template.md con los valores de
# results/final/tables/network_benchmark_key_numbers.tsv y escribe
# docs/supplementary_network_method.md. Falla si la plantilla cita un id
# que no está en la tabla: ningún número del texto se transcribe a mano.
#
# Formato: {{id}} usa el valor tal cual (3 cifras significativas);
# {{id|pct}} lo escribe como porcentaje; {{id|int}} redondea a entero.
#
# Uso: Rscript scripts/network_benchmark/09_render_text.R

source("scripts/network_benchmark/lib_benchmark.R")

config <- bench_config()
dirs <- bench_dirs(config)

key <- read.delim(file.path(dirs$tables, "network_benchmark_key_numbers.tsv"), stringsAsFactors = FALSE)
values <- setNames(key$value, key$id)

template <- readLines("docs/supplementary_network_method.template.md", encoding = "UTF-8")
text <- paste(template, collapse = "\n")
text <- gsub("(?s)<!--.*?-->\n*", "", text, perl = TRUE)   # notas internas de la plantilla

tokens <- unique(regmatches(text, gregexpr("\\{\\{[^}]+\\}\\}", text))[[1]])
missing <- character(0)
for (tok in tokens) {
  spec <- strsplit(gsub("[{}]", "", tok), "|", fixed = TRUE)[[1]]
  id <- spec[1]; fmt <- if (length(spec) > 1) spec[2] else ""
  if (!id %in% names(values) || is.na(values[[id]])) { missing <- c(missing, id); next }
  v <- values[[id]]
  out <- switch(fmt,
    pct = paste0(formatC(100 * v, format = "f", digits = if (abs(100 * v) < 1) 2 else if (abs(100 * v) < 10) 1 else 0), " %"),
    int = format(round(v), big.mark = ","),
    format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
  text <- gsub(tok, out, text, fixed = TRUE)
}
if (length(missing)) stop("ids de la plantilla sin valor en key_numbers: ", paste(unique(missing), collapse = ", "))

writeLines(text, "docs/supplementary_network_method.md", useBytes = TRUE)
message("docs/supplementary_network_method.md escrito (", length(tokens), " valores insertados)")
