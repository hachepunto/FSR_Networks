# FSR_Networks

Código y manuscrito de *Homogeneous network inference for subtype-stratified master regulator analysis in cancer* (Tovar, Barrios, Pérez-Calixto, Hernández-Lemus). El artículo compara un solo método de red para todas las cohortes, fixed-size Spearman regulons (FSR, con `corto`), contra información mutua (ARACNe3) en análisis de reguladores maestros transcripcionales con NaRnEA sobre cohortes de TCGA estratificadas por cáncer y subtipo.

## Contenido

- `scripts/` — preparación de datos (pasos 1–5 abajo), compilación y corrida de ARACNe3 (`setup_aracne3.sh`, `run_aracne3.sh`) y `lib_common.R`, la biblioteca compartida.
- `scripts/network_benchmark/` — pasos 00–15b del benchmark (`run_all.sh` encadena 01–09 y 12/12b; `lib_benchmark.R` con las funciones comunes).
- `config/config.yaml` — única fuente de configuración: rutas, cohortes TCGA, clasificación de subtipos ya congelada, parámetros de redes, DESeq2 y benchmark. `config/env_network_benchmark.yaml` — entorno conda del benchmark.
- `results/final/` — tablas y figuras finales del benchmark (`network_benchmark_*`), entre ellas `network_benchmark_key_numbers.tsv`.
- `docs/supplementary_network_method*.md` — nota suplementaria (plantilla y versión renderizada).
- `manuscript/` — `FSR_TMR.tex` (texto principal), `numbers.tex` y `figures/`. Conectado a Overleaf por la sincronización con GitHub; el documento principal es `manuscript/FSR_TMR.tex`.

Fuera de git (en `.gitignore`): `data/` (descargas de TCGA en `data/tcga`, intermedios en `data/work`) y `software/ARACNe3`. Todos los scripts asumen que el directorio de trabajo es la raíz del repositorio.

## Ejecución

```bash
# 1–5. Datos: TCGA → subtipos → matrices → firmas de DESeq2
Rscript scripts/download_data.R          # RNA-seq STAR-Counts, clínica y subtipos de los proyectos de config (≈30 GB)
Rscript scripts/prepare_tfs.R            # lista de TFs de Lambert et al. 2018
Rscript scripts/classify_subtypes.R      # subtipo por muestra y n por cáncer × subtipo
Rscript scripts/preprocess_expression.R  # matrices log2(TPM+1) por analysis_id y por cáncer completo
for id in <analysis_id>...; do Rscript scripts/run_deseq2.R $id; done   # firma tumor vs. normal (34 analysis_id)

# ARACNe3 (califano-lab/ARACNe3 en el commit de config.yaml)
bash scripts/setup_aracne3.sh

# Benchmark
Rscript scripts/network_benchmark/00_install_r_packages.R   # corto y NaRnEA a versiones fijas
bash scripts/network_benchmark/run_all.sh [corto_jobs] [aracne3_jobs] [aracne3_threads]

# Números y figuras del manuscrito (generados; no editar a mano)
Rscript scripts/network_benchmark/10_appendix_statistics.R manuscript
Rscript scripts/network_benchmark/11_appendix_numbers_tex.R manuscript
```

Compilar el manuscrito: `make -C manuscript` (latexmk) o `tectonic manuscript/FSR_TMR.tex`.

Los pasos 1–5 corren en un entorno R con TCGAbiolinks 2.38.0, SummarizedExperiment, DESeq2, CMScaller y GCclassifier (GitHub, commits fijados), org.Hs.eg.db, matrixStats y yaml; ese entorno todavía no tiene receta versionada. El benchmark usa `config/env_network_benchmark.yaml`.

## Origen

Extraído de `mra_multicancer_pipeline` (octubre de 2026) con su historial de git. El pipeline de MRA no depende de este repositorio.
