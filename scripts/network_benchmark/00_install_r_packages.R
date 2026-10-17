# Instala en el entorno activo (config/env_network_benchmark.yaml) los dos
# paquetes que no están en conda, a versiones fijas. Idempotente.
#
#   corto 1.3.1   — CRAN (misma versión que usa el pipeline, tcga_env).
#   NaRnEA 0.2.0  — GitHub califano-lab/NaRnEA, commit fijado abajo.
#                   install_github() falla: el campo License del
#                   DESCRIPTION es un placeholder de usethis sin editar que
#                   rompe el parser DCF. Se baja el tarball, se reemplaza
#                   ese campo por "GPL-3" y se instala desde la fuente.

CORTO_VERSION <- "1.3.1"
NARNEA_COMMIT <- "36b8fdbc5312bc709a9af90abbe8c6f3f44b1575"

if (!requireNamespace("corto", quietly = TRUE) ||
    as.character(packageVersion("corto")) != CORTO_VERSION) {
  remotes::install_version("corto", version = CORTO_VERSION,
                           repos = "https://cloud.r-project.org", upgrade = "never")
}

# Si ya está instalado, se exige el mismo commit (lo registra la instalación
# de abajo en el DESCRIPTION como RemoteSha).
narnea_ok <- requireNamespace("NaRnEA", quietly = TRUE) &&
  identical(packageDescription("NaRnEA")$RemoteSha, NARNEA_COMMIT)
if (!narnea_ok) {
  tmpdir <- tempfile()
  tarball <- tempfile(fileext = ".tar.gz")
  download.file(sprintf("https://github.com/califano-lab/NaRnEA/archive/%s.tar.gz", NARNEA_COMMIT),
                destfile = tarball, mode = "wb")
  untar(tarball, exdir = tmpdir)
  pkg_dir <- list.dirs(tmpdir, recursive = FALSE)[1]

  desc_path <- file.path(pkg_dir, "DESCRIPTION")
  desc <- readLines(desc_path)
  i <- grep("^License:", desc)
  stopifnot(length(i) == 1)
  if (grepl("^\\s+license\\s*$", desc[i + 1])) desc <- desc[-(i + 1)]
  desc[i] <- "License: GPL-3"
  desc <- c(desc, paste0("RemoteSha: ", NARNEA_COMMIT))
  writeLines(desc, desc_path)

  install.packages(pkg_dir, repos = NULL, type = "source")
}

stopifnot(as.character(packageVersion("corto")) == CORTO_VERSION,
          identical(packageDescription("NaRnEA")$RemoteSha, NARNEA_COMMIT))
message("corto ", packageVersion("corto"), " y NaRnEA ", packageVersion("NaRnEA"), " listos")
