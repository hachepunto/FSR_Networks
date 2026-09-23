# ============================================================
# PATHS
# ============================================================

library(yaml)

config <- yaml::read_yaml("config/config.yaml")

data_raw_dir <- file.path(config$paths$data_root, "raw", "tfs")
aracne_input_dir <- file.path(config$paths$work_root, "networks", "input")

lambert_raw_file <- file.path(data_raw_dir, "lambert2018_DatabaseExtract_v1.01.csv")
tf_file <- file.path(aracne_input_dir, "TFs_Lambert2018.txt")
lambert_url <- "https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv"

# Crear unicamente los directorios donde se escribiran archivos
dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aracne_input_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PREPARACION DE FACTORES DE TRANSCRIPCION
# Lambert et al. 2018
#
# Descarga la base de datos de Lambert si no existe y genera
# la lista de factores de transcripcion utilizada por ARACNe-AP.
# ============================================================

# ============================================================
# ENTRADAS
#
# Archivo descargado desde Lambert et al. 2018:
#   data/raw/lambert2018_DatabaseExtract_v1.01.csv
#
# URL:
#   https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
# ============================================================

# ============================================================
# SALIDAS
#
# Lista de factores de transcripcion:
#   networks/input/TFs_Lambert2018.txt
# ============================================================

# ============================================================
# SECCION 1: Descargar la base de datos Lambert
# ============================================================

if (file.exists(lambert_raw_file)) {
    message("La base de datos Lambert ya existe. Se omite la descarga.")
} else {
    message("Descargando la base de datos Lambert...")

    download.file(
        url = lambert_url,
        destfile = lambert_raw_file,
        mode = "wb"
    )

    message(lambert_raw_file, " descargado")
}

# ============================================================
# SECCION 2: Generar la lista de factores de transcripcion
# ============================================================

lambert_db <- read.csv(
    lambert_raw_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_columns <- c("HGNC symbol", "Is TF?")

if (!all(required_columns %in% colnames(lambert_db))) {
    stop(
        "La base de datos Lambert no contiene las columnas requeridas: ",
        paste(required_columns, collapse = ", ")
    )
}

lambert_symbols <- lambert_db[["HGNC symbol"]][
    lambert_db[["Is TF?"]] == "Yes"
]

lambert_symbols <- unique(trimws(lambert_symbols))
lambert_symbols <- lambert_symbols[
    !is.na(lambert_symbols) &
    lambert_symbols != ""
]

if (length(lambert_symbols) == 0) {
    stop("No se encontraron factores de transcripcion en la base Lambert")
}

writeLines(lambert_symbols, tf_file)

message(
    tf_file,
    " escrito - ",
    length(lambert_symbols),
    " factores de transcripcion"
)
