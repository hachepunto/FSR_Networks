# FSR_Networks

Code, results and manuscript for **"Homogeneous network inference for subtype-stratified master regulator analysis in cancer"** (Hugo Tovar, Camila Barrios, Daniel Pérez-Calixto, Enrique Hernández-Lemus; manuscript in preparation).

This repository contains everything needed to regenerate the manuscript from public data. It downloads TCGA, classifies tumours into molecular subtypes, infers the regulatory networks, runs the benchmark, and writes every number and figure that the manuscript cites.

## Status

The results, tables, figures and manuscript in `main` are frozen at the tag `congelado-2026-11-05-pam50-tcga`; this is the version submitted for publication. They use the subtype rules of `config/config.yaml` in `main`. For BRCA this is the PAM50 column of `TCGAquery_subtype()` (Basal 190, Her2 82, LumA 562, LumB 209). On 2026-11-04 the MRA pipeline (commit `d054286`) replaced that column with PAM50 recomputed with genefu and ER-balanced centring (192 / 91 / 543 / 239; 44 tumours change subtype). That change is not incorporated here. The branch `brca-pam50-genefu` carries the synchronised classification code for a future re-run of the benchmark; none of its results have been recomputed.

## The question

Transcriptional master regulator (TMR) analysis with NaRnEA needs one gene regulatory network per cohort. When TCGA cancers are stratified by molecular subtype, cohort sizes range from 27 to 562 tumours. Mutual-information networks inferred with ARACNe3 are the usual choice, but they cannot be inferred for the smallest cohorts. Filling the gap with a correlation network would mix two inference methods inside the very comparison between subtypes that the analysis is meant to make.

The manuscript proposes inferring every network with a single procedure, **fixed-size Spearman regulons (FSR)**. It then quantifies what this choice concedes and what it provides with respect to ARACNe3, labelling each claim as a standard result, a heuristic, or an empirical result. FSR is presented as an adequate and uniform choice for multi-cohort analysis, not as a statistically superior one.

## Main results

- **Agreement with ARACNe3.** Where both methods are feasible (n ≥ 95, 15 cohorts), FSR and ARACNe3-top100 share a median 70 % of the top-20 TMRs. This is the same agreement as between two filterings of the same ARACNe3 network, and far above a null in which the network comes from another cancer. The rank correlation of the full NES vector shows a measurable method effect in the tail of the ranking, not at its head.
- **Synthetic data with known regulators.** From n = 100 on, FSR recovers the planted master regulators within 0.06 of the true-network ceiling. ARACNe3 produces no edges at n = 30 and too few to evaluate at n = 60.
- **Minimum cohort size.** Subsampling large real cohorts shows that networks from fewer than 60 tumours do not reliably reproduce their own cohort's TMRs, and that 60–100 is a grey zone. A cancer-wide network does not rescue small subtypes.
- **What FSR concedes.**
  - Rank correlation is more powerful than the MI estimators we simulated for n ≲ 200, but it cannot detect non-monotone dependence. In a synthetic worst case this costs about 0.05–0.10 of the master regulators recovered.
  - The per-target data-processing inequality (DPI) is a heuristic, with a measurable disagreement rate against the triangle DPI used by ARACNe.

## Methods in brief

**FSR.** FSR is inferred with [`corto`](https://cran.r-project.org/package=corto) 1.4.0, using one set of parameters for every cohort, whatever its size:

1. Expression is converted to per-gene ranks, so Pearson correlation equals Spearman correlation.
2. An effect-size threshold of |ρ| ≥ 0.2 is applied identically for every n. The p-value passed to `corto` is derived from `r_min` and n, which makes the effective correlation threshold independent of cohort size.
3. Per-target DPI pruning (`dpi = "target"`): for each target, only the regulator with the largest |ρ| is kept in each bootstrap.
4. 100 bootstraps; the support of an edge is the fraction of bootstraps that retain it.
5. Each regulator is truncated to its 100 best-supported targets (`top_k`).

**Comparator.** [ARACNe3](https://github.com/califano-lab/ARACNe3), pinned to commit `f0fa4fb`, run with 100 subnetworks, seed 42 and α = 0.05. This commit is the merge of califano-lab PR #3, which fixes data races in the MI computation and makes the null model reproducible. ARACNe3 has no dedicated publication; its algorithm is described in Griffin et al. (2023, *Entropy* 25:542), Section 2.2. The benchmark evaluates three versions of the ARACNe3 network:
- the full consolidated network;
- ARACNe3-top100, truncated like FSR;
- a support-filtered network (≥ 10 of 100 subnetworks).

The benchmark also evaluates the previous rule of the pipeline: `corto` on log2(TPM + 1) at a fixed p < 10⁻⁸, whose effective correlation threshold depends on n.

**Data and signatures.**
- **Cohorts:** ten TCGA projects (BRCA, LUAD, LUSC, THCA, PRAD, COAD, KIRC, STAD, UCEC, LIHC), split into 34 cancer × subtype cohorts. Subtype assignment rules are frozen in `config/config.yaml`.
- **Expression:** STAR counts, one tumour sample per patient, log2(TPM + 1).
- **Signatures:** the DESeq2 Wald statistic of tumour against solid-tissue normal.
- **Regulators:** the human transcription factors of Lambert et al. (2018).
- **TMR inference:** NaRnEA 0.2.0 with a minimum of 30 targets per regulator.

**Benchmark design.**

| Analysis | What it measures | Steps |
|---|---|---|
| Real cohorts | Regulon structure, enrichment for CollecTRI, TMR agreement between schemes (top-k overlap, ρ of NES) | 01–04 |
| Biological reference and null | Agreement between subtypes of the same cancer (network and signature components) and with networks from another cancer | 04b |
| Synthetic data | 200 TFs in 20 correlated modules, 5,000 targets, linear / monotone non-linear / non-monotone edges, n = 30–400, two signal strengths, 5 replicates; edge precision and recall, and recovery of 10 planted master regulators | 05–07 |
| Subsampling | Smallest n at which a network reproduces the TMRs of the full cohort, in 5 large real cohorts | 12, 12b |
| Cancer-wide network | Whether pooling all subtypes of a cancer helps a small subtype | 13, 13b |
| Triangle DPI | FSR with ARACNe's triangle rule (`corto` `dpi = "triangle"`) instead of the per-target rule | 14, 14b |
| Robustness to k | Every comparison repeated for k = 5–50, with bootstrap confidence intervals clustered by cohort or by cancer | 08 |
| Technical appendix | Estimator efficiency, empirical validity of the DPI on real triplets, indirect edges against partial correlation (graphical lasso), contribution of non-monotone edges to NES | 10 |

## Repository layout

```
config/
  config.yaml                    single source of configuration: paths, TCGA projects, frozen subtype rules,
                                 network parameters, DESeq2 and benchmark settings (comments in Spanish)
  env_network_benchmark.yaml     conda environment for the benchmark (env name: mra_bench)
scripts/
  download_data.R ... run_deseq2.R   data preparation (steps 1–5 below)
  setup_aracne3.sh, run_aracne3.sh   build ARACNe3 at the pinned commit; run one network
  lib_common.R                       functions shared with the MRA pipeline (rank transform, r -> p, top-k)
  network_benchmark/                 benchmark steps 00–15b, run_all.sh and lib_benchmark.R
results/final/
  tables/network_benchmark_*.tsv     final tables; network_benchmark_key_numbers.tsv holds every number cited
  figures/network_benchmark_*.pdf    supplementary figures S1–S9
docs/
  supplementary_network_method.template.md   Supplementary Note with {{key}} placeholders
  supplementary_network_method.md            rendered Supplementary Note (generated)
manuscript/
  FSR_TMR.tex                    main text and technical appendix
  numbers.tex                    LaTeX macros (\K...) with every number in the text (generated)
  figures/                       appendix figures A1–A6 (generated)
  appendix_*.tsv                 appendix tables (generated)
```

Not under version control: `data/` (TCGA downloads in `data/tcga`, intermediates in `data/work`, about 30 GB) and `software/ARACNe3`. All scripts assume the repository root as working directory.

## No number is typed by hand

Every number in the Supplementary Note and in the manuscript is generated:
1. `08_summarize.R` writes `results/final/tables/network_benchmark_key_numbers.tsv`.
2. `09_render_text.R` fills the Supplementary Note template from that table, and fails if the template cites a key that does not exist.
3. `10_appendix_statistics.R` and `11_appendix_numbers_tex.R` write the appendix tables and figures and turn all numbers into macros in `manuscript/numbers.tex`.

For this reason, `numbers.tex`, `manuscript/figures/` and `docs/supplementary_network_method.md` must not be edited by hand: they are overwritten every time the results are regenerated. To change a number, change the analysis. To change the wording, edit the template or `FSR_TMR.tex`.

`08_summarize.R` also contains guards that stop the run if a qualitative statement in the text no longer holds for the data. For example, it checks that the ordering methods > subtypes > null fails only at the smallest k.

## Reproducing the results

**1. Data (steps 1–5).** These steps need an R environment with TCGAbiolinks 2.38.0, SummarizedExperiment, DESeq2, ggplot2, pheatmap, matrixStats, yaml, CMScaller and GCclassifier (GitHub, commits pinned in the scripts). This environment does not yet have a versioned recipe.

```bash
Rscript scripts/download_data.R          # STAR counts, clinical data and subtype annotations (~30 GB)
Rscript scripts/prepare_tfs.R            # TF list of Lambert et al. 2018
Rscript scripts/classify_subtypes.R      # subtype per sample and n per cancer x subtype
Rscript scripts/preprocess_expression.R  # log2(TPM+1) matrices per cohort and per cancer
for id in <analysis_id>...; do Rscript scripts/run_deseq2.R $id; done   # tumour vs normal signature, 34 cohorts
```

**2. ARACNe3.** Clones califano-lab/ARACNe3 into `software/ARACNe3`, checks out the pinned commit and builds it. It needs CMake ≥ 3.5 and a C++ compiler with OpenMP. `run_aracne3.sh` refuses to run a binary built from any other commit.

```bash
bash scripts/setup_aracne3.sh
```

**3. Benchmark.** The arguments of `run_all.sh` are:
- `corto_jobs`: parallel `corto` runs, each with `network_benchmark.threads` threads (default 8);
- `aracne3_jobs`: parallel ARACNe3 runs;
- `aracne3_threads`: threads per ARACNe3 run.

```bash
conda env create -f config/env_network_benchmark.yaml && conda activate mra_bench
Rscript scripts/network_benchmark/00_install_r_packages.R   # corto 1.4.0 (CRAN) and NaRnEA 0.2.0 (GitHub, pinned commit)
bash scripts/network_benchmark/run_all.sh 6 2 18
```

`run_all.sh` runs steps 01–09, 12/12b, 13/13b and 14/14b, and stops if `corto` is older than 1.4.0.

Network steps (02, 03, 06, 12, 13, 14) skip outputs that already exist, so a run can be resumed. If parameters change, delete the affected outputs first, because the scripts do not compare saved parameters with the config.

ARACNe3 dominates the computing cost. In our run, the 15 real ARACNe3 networks took 5.5 h of wall-clock time at 2 × 18 threads (about 100 core-hours), against 1.8 h for FSR on the same cohorts.

Steps 15/15b (subtype-vs-rest signatures) are a side analysis that is not cited in the manuscript, and `run_all.sh` does not run them.

A different configuration file can be used through `MRA_CONFIG=/path/to/config.yaml`. All scripts honour it.

**4. Manuscript.**

```bash
Rscript scripts/network_benchmark/10_appendix_statistics.R manuscript
Rscript scripts/network_benchmark/11_appendix_numbers_tex.R manuscript
make -C manuscript            # latexmk; or: tectonic manuscript/FSR_TMR.tex
```

## Editing the manuscript in Overleaf

The repository is linked to an Overleaf project through Overleaf's GitHub synchronisation. The main document is `manuscript/FSR_TMR.tex`, compiled with pdfLaTeX. Synchronisation is manual:
- *Menu → GitHub → Pull* before editing;
- *Push* at the end of a session.

Overleaf drops the executable bit of shell scripts on push. This is harmless, because the README always invokes the scripts through `bash` or `Rscript`.

## Software versions

| Software | Version |
|---|---|
| R | 4.5 |
| corto | 1.4.0 (CRAN) |
| NaRnEA | 0.2.0 (califano-lab, commit `36b8fdb`) |
| DESeq2 | 1.50 |
| ARACNe3 | califano-lab commit `f0fa4fb` |
| TCGAbiolinks | 2.38.0 |
| CollecTRI | OmniPath API (download date recorded in `network_benchmark_collectri_provenance.tsv`) |

## Known limitations of reproducibility

- **Data preparation (steps 1–5).**
  - The R environment of steps 1–5 is not yet pinned by a lock file.
  - The GDC data release is not pinned. Downloads of a later release may change cohort sizes. `preprocess_expression.R` warns when they change but does not stop.
- **CollecTRI.** CollecTRI is served by the OmniPath API without versioning. The file used is archived as `results/final/tables/network_benchmark_collectri_edges.tsv`.

## Origin

This repository was extracted, with its git history, from the multi-cancer MRA pipeline [`mra_multicancer_pipeline`](https://github.com/camilasauria/mra_multicancer_pipeline) in October 2026. The MRA pipeline does not depend on this repository.
