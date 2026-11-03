# Supplementary Note: choice of the gene regulatory network inference method

## Summary

Master regulator (MR) analysis needs a tumour-specific regulatory network (regulon) for each of the 34 cancer × subtype cohorts in this study, whose sizes range from 27 to 562 primary tumours. The previous design inferred networks with ARACNe3 (mutual information, MI) when n ≥ 95 and with corto (correlation) at a fixed p-value when n < 95. We compared that design with a single procedure for every cohort:

- **fixed-size Spearman regulons (FSR)**: corto on rank-transformed expression, |r| ≥ 0.2 before the data-processing inequality, and the 100 best-supported targets per regulator;
- **ARACNe3-top100**: ARACNe3 harmonised to the same 100 targets per regulator, for the 15 cohorts with n ≥ 95;
- **corto at fixed p** (p < 10⁻⁸ on log2 TPM), the rule previously used for n < 95, here applied to every cohort for comparison.

Evaluation used the pipeline's own DESeq2 signatures and NaRnEA on all real cohorts, and synthetic data with a known network. We adopt FSR. In brief:

- **Mixing methods would confound between-subtype comparisons.** Changing the inference method within a cohort (FSR vs ARACNe3-top100, same signature) turns over 30 % of the top-20 MRs. Changing only the network between two subtypes of the same cancer (same signature, FSR network of the other subtype) turns over 45 %, and two subtypes analysed each with their own signature and network differ by 55 %. The method effect is smaller than the network component of the biological effect, and about half of the total. A single method removes it entirely.
- **Where both methods can be run, FSR and ARACNe3-top100 agree as closely as two filterings of the same ARACNe3 network.** Median top-20 overlap 70 % (FSR vs ARACNe3-top100) against 70 % (ARACNe3-top100 vs the support-filtered ARACNe3 network); both far above a null in which the network comes from another cancer (median 35 %, maximum 60 %).
- **corto at a fixed p-value does not give comparable regulons across cohort sizes**, which was the rule's weakness for n < 95: regulon size tracks n, and in cohorts with n < 95 a median of only 518 regulators reach the 30 targets NaRnEA needs, against 1,270 under FSR.
- **In synthetic data with known MRs, FSR recovers them at every cohort size** (within 0.06 of the true-network upper bound from n = 100), whereas ARACNe3 and corto at fixed p recover none at n ≤ 60. Synthetic data are, however, uninformative about FSR versus ARACNe3-top100 at the sizes where both apply in real data (Section 3.4).
- **Non-monotonic dependencies**: ARACNe3 detects them in synthetic data at large n. In real tumours we cannot tell whether they are present; a positive control shows that our statistics would not distinguish the real networks from synthetic ones containing 20 % non-monotone edges. Their contribution to NaRnEA scores is bounded by design, because the association mode is a Spearman correlation.
- **Subtype-specific networks are the primary design; a cancer-wide network is the sensitivity analysis.** For cohorts with n ≥ 60 the two share 65 % of the top-20 TMRs, close to the agreement between methods, and the differences between subtypes are carried mostly by their signatures (Section 3.8). TMRs found with both networks are reported as the robust core.
- **Computing cost**: 14.6 core-hours for FSR against 99.5 for ARACNe3 on the same 15 cohorts.

## 1. Background and question

The pipeline infers one network per cohort from the tumour samples of that cohort, then scores MRs with NaRnEA, using the DESeq2 Wald statistic of tumour versus matched normal tissue as the signature. Results are compared across subtypes and cancers. Three questions:

1. Does a fixed p-value threshold yield regulons that are comparable across cohort sizes?
2. Does mixing two inference methods introduce method-driven differences between subtypes, and is anything lost by using a correlation-based method everywhere?
3. MI can detect non-monotonic dependencies that a correlation cannot. Does this matter for the networks and MRs we infer?

## 2. Methods

### 2.1 Network schemes

All schemes use the same tumour expression matrix as the pipeline (log2(TPM + 1), protein-coding genes, expression filters from `preprocess_expression.R`) and the same regulator list (Lambert et al. 2018).

- **corto at fixed p.** corto 1.3.1 on log2(TPM + 1), 100 bootstraps, p < 10⁻⁸. The effective correlation threshold therefore depends on n: |r| ≥ 0.859 at n = 27, |r| ≥ 0.549 at n = 94, |r| ≥ 0.239 at n = 562.
- **Fixed-size Spearman regulons (FSR).** corto 1.3.1 on per-gene ranks, so that Pearson correlation equals Spearman correlation; 100 bootstraps; a threshold |r| ≥ r_min applied before the data-processing inequality (DPI), with the p-value derived for each n so that the effect-size threshold is constant; then, per regulator, the k = 100 targets with the highest bootstrap likelihood (ties broken by |r|). corto applies the DPI in a simplified per-target form rather than ARACNe's triangle form: for each target and bootstrap, among all regulators that pass the threshold only the one with the largest |r| survives ("winner takes all"), whether or not the competing regulators are co-expressed with each other, and the likelihood is the fraction of bootstraps in which an edge wins. It is therefore a measure of competition between regulators for a target, not of edge stability alone; a target with two independent true regulators keeps only the stronger one, which bounds recall. r_min = 0.2 was chosen after an initial run with 0.4, by a grid over {0, 0.1, 0.2, 0.3, 0.4} on the same data (Section 3.5). k was varied in {50, 100, 150}.
- **ARACNe3** (commit f0fa4fb of github.com/califano-lab/ARACNe3, which includes the fix of the data races in the MI computation): 100 subnetworks, seed 42, α = 0.05. Three versions were evaluated: the full consolidated network; edges supported by ≥ 10 of 100 subnetworks; and **ARACNe3-top100**, the 100 best-supported targets per regulator (subnetwork count, then MI).

In NaRnEA regulons, the association mode (`am`) is the correlation (Spearman for ARACNe3, rank-Pearson for corto), and the association weight (`aw`) is the within-regulator rank of edge support scaled to (0, 1]. As a sensitivity analysis, `aw` was also set to the raw bootstrap likelihood.

### 2.2 Real data

All 34 cohorts of the pipeline (10 TCGA projects; subtype definitions in `subtype_classification_methods.tsv`). Signatures are the DESeq2 Wald statistics from `run_deseq2.R`. MRs were scored with NaRnEA 0.2.0 (`matrixNaRnEA`, minimum 30 targets). Without a ground-truth network, we used four kinds of evidence:

- **External reference.** CollecTRI (curated TF → target interactions, via OmniPath; the retrieved file is archived with the results). *CollecTRI enrichment* is the fraction of a regulon's edges found in CollecTRI, divided by the fraction expected if every regulator chose the same number of targets at random among expressed genes (1 = random). It depends on regulon size, so comparisons are made at equal k and over the same cohorts.
- **Agreement of MR results between schemes.** Spearman correlation of NES over shared regulators and overlap of the top-20 and top-50 regulators by |NES|. Two references calibrate these numbers: a *within-method* baseline (ARACNe3-top100 vs the support-filtered ARACNe3 network, two filterings of the same network) and a *null* (NaRnEA with the cohort's own signature but the FSR network of a cohort from another cancer; 102 pairs).
- **Biological effect size, decomposed.** The same agreement statistics between subtypes of the same cancer, each with its own signature and FSR network (43 pairs): the magnitude of the differences the pipeline is designed to detect. Because FSR vs ARACNe3-top100 changes only the network while holding the signature fixed, the matched yardstick is the *network component* of that difference: own signature with the FSR network of another subtype of the same cancer (86 directed pairs). The complementary swap (own network, other subtype's signature) gives the signature component.
- **Sample-size requirement, measured on real data by subsampling.** For the 5 largest cohorts (n = 190–562), FSR networks were inferred from random subsets of n ∈ {20, 27, 40, 60, 80, 100, 150} tumours (5 replicates each, without replacement) and scored with NaRnEA against the *same* DESeq2 signature as the full cohort, so that only the effect of n on the network is measured. Agreement with the full-cohort result (top-20/50 overlap, ρ of NES) and per-regulator edge Jaccard were recorded; the A-vs-B, within-cancer and null agreements of Section 3.2 serve as yardsticks.
- **Non-monotonic dependencies.** In ARACNe3 networks, the fraction of edges with |Spearman| < 0.2 (the same effect-size threshold as FSR's r_min, i.e. edges that FSR could never select) at each support level, and the number of top-quartile-MI edges with |Spearman| < 0.2. The same statistics were computed on the synthetic ARACNe3 networks, where 20 % of true edges are non-monotone, as a positive control.

### 2.3 Synthetic data

Expression was simulated from a known network (`05_simulate_synthetic.R`): 200 regulators in 20 modules with within-module latent correlation 0.5 (creating indirect dependencies the DPI should remove); 5,000 targets with 1–3 regulators each; edges linear (50 %), monotone non-linear (30 %, saturating sigmoid) or non-monotone (20 %, U-shaped, Spearman ≈ 0); negative-binomial counts with variable library sizes. Networks were inferred from a "tumour" group of n ∈ {30, 60, 100, 200, 400}. Ten true MRs, chosen among regulators with ≥ 30 true targets, were shifted by ±1 SD in the tumour group relative to a reference group (n = 50); the signature is the DESeq2 Wald statistic of tumour vs reference. Two noise levels gave a median |Spearman| of true linear edges of about 0.35 ("weak") and 0.46 ("strong"); the two regimes share networks and MRs (paired design), with 5 replicates per n and regime.

Edge-level metrics: precision, recall by edge type, false positives that are *indirect* (regulator in the same module as a true regulator of the target), reported per 100 inferred edges, and average precision of each method's full edge ranking. MR recovery: the fraction of the 10 true MRs among the top 10 by |NES| (a true MR with too small a regulon to be evaluated counts as missed; a random top-10 scores 0.05), and the AUROC of |NES| over all 200 regulators, with non-evaluable regulators scored 0. The true network is included as an upper bound. With 5 replicates, the standard error of the top-10 fraction (FSR) ranges from 0.02 to 0.0775, so differences below ~0.1 are not resolved.

## 3. Results

### 3.1 corto at a fixed p-value yields regulons whose size depends on cohort size (Fig. S1)

Restricting to the 19 cohorts with n < 95, where this rule was in use: the median regulon size ranges from 2 to 102 targets and correlates with n (Spearman ρ = 0.779); a median of 518 regulators reach 30 targets; and in 17 of these cohorts fewer than 1 % of edges have negative mode (median 0.28 %), because only very strong, almost always positive, correlations pass the threshold.

Under FSR, 32 of 34 cohorts reach the 100-target cap, the minimum median is 53 targets (THCA BRAF-like, n = 358, and BRCA LumA, n = 562, where the DPI prunes more edges), and a median of 1,270 regulators can be evaluated in the small cohorts. Uniform size is, however, a property of the cap, not of the data: under FSR the CollecTRI enrichment still increases with n (Spearman ρ = 0.756), from 1.62–1.82 in the three cohorts with n < 30 to a median of 2.52 in the cohorts with n ≥ 95. Comparable size does not mean comparable quality.

**CollecTRI enrichment over the same 15 cohorts (n ≥ 95), at k = 100 where applicable:**

| Network | Enrichment |
|---|---|
| corto at fixed p | 2.1 |
| FSR | 2.52 |
| Full ARACNe3 network | 1.13 |
| ARACNe3, support ≥ 10 | 2.19 |
| ARACNe3-top100 | 2.53 |

FSR and ARACNe3-top100 are equivalent on this criterion. The full ARACNe3 network is barely above random, because most of its edges appear in a single subnetwork. (The expected fraction counts TF → TF interactions in CollecTRI, which corto cannot recover because it excludes regulators from the target set; the bias is against corto and small.) Enrichment also rises as k falls (FSR: 2.8 at k = 50, 2.3 at k = 100, 2.07 at k = 150, all cohorts), so the value for corto at fixed p, whose regulons are smaller, is not directly comparable.

### 3.2 Method effect versus biological effect (Fig. S2)

Across the 15 cohorts with n ≥ 95:

| Comparison | Median ρ of NES | Median top-20 overlap | Median top-50 overlap |
|---|---|---|---|
| FSR vs ARACNe3-top100 (between methods) | 0.926 (min 0.853) | 70 % (min 50 %) | 70 % |
| ARACNe3-top100 vs ARACNe3 support ≥ 10 (within method) | 0.985 | 70 % | 72 % |
| ARACNe3-top100 vs full ARACNe3 network (within method) | 0.928 | 60 % | 60 % |
| Own signature, FSR network of another subtype of the same cancer (network component of the biological effect) | 0.707 | 55 % | 52 % |
| Own network, signature of another subtype of the same cancer (signature component) | 0.822 | 60 % | 60 % |
| Subtypes of the same cancer, each with own signature and network (total biological effect) | 0.575 | 45 % | 40 % |
| Own signature, network from another cancer (null) | 0.54 | 35 % | 35 % |

Three observations:

- The top-k overlap between methods equals the overlap between two filterings of one ARACNe3 network, and exceeds that between ARACNe3-top100 and the full network. The genome-wide ρ of NES is systematically lower between methods than within ARACNe3 (in 15 of 15 cohorts), as expected when the two ARACNe3 versions share most edges.
- Both are far above the null, which never exceeds a top-20 overlap of 60 % (within-cancer vs null, Wilcoxon p = 0.0211; interquartile ranges 22 %–50 % and 25 %–40 %). The null is a conservative "wrong tissue" baseline rather than a random one, because co-expression modules are partly shared across cancers.
- **The method effect is smaller than, but not negligible against, the network component of the biological effect.** Switching method within a cohort (signature fixed) turns over 30 % of the top-20 MRs. Swapping only the network between two subtypes of the same cancer (signature fixed) turns over 45 %; swapping only the signature, 40 %; and the total difference between subtypes is 55 %. A design that infers networks with ARACNe3 for some subtypes and corto for others (as in COAD, where CMS2 and CMS4 have n ≥ 95 but CMS1 and CMS3 do not) would therefore add a method-driven turnover about two-thirds the size of the network differences it sets out to measure, and more than half of the total between-subtype difference.

Agreement did not degrade near the n = 95 boundary (ρ = 0.933 at n = 96). Agreement with MR results from the curated CollecTRI regulons was low for every scheme (median ρ: corto at fixed p 0.24, FSR 0.21, ARACNe3-top100 0.182, full ARACNe3 0.173), reflecting how different data-driven and literature regulons are; it does not discriminate between schemes.

### 3.3 Non-monotonic dependencies (Figs. S3, S4)

**Synthetic data.** ARACNe3 recovered non-monotone edges at large n: recall 0.82 (strong signal) and 0.401 (weak) at n = 400, but only 0.0463 at n = 100 under strong signal and 0.0256 at n = 200 under weak signal. FSR did not (0.0733 at n = 400). Its average precision over the full edge ranking was higher than ARACNe3's at every n except n = 400 under strong signal (0.689 vs 0.801).

**Real data, and why we cannot conclude that non-monotonic edges are absent.** Among ARACNe3 edges inferred from real tumours, the fraction with |Spearman| < 0.2 falls with support: 39 % over all edges, 4.3 % at support ≥ 10 (range 0.74 %–8.6 %), 0.06 % at support ≥ 50; and only 66 of 597,841 top-quartile-MI edges with support ≥ 10 have |Spearman| < 0.2. Taken alone, this would suggest that MI detects few well-supported non-monotonic edges in these cohorts. The positive control shows otherwise: on synthetic ARACNe3 networks with 20 % true non-monotone edges, the same fraction at support ≥ 10 is 5.8 % (strong signal, n = 400) and 2.0 % (weak, n = 400), and at the sizes of most real cohorts it is lower still (0.03 % and 2.8 % at n = 100 and 200, strong signal). True non-monotone edges have lower MI and lower support than linear ones (they make up 6.3 % of supported edges in the strong regime), so support-based statistics under-represent them even when they are abundant. The statistic is therefore nearly blind to non-monotone edges at n ≈ 100–200, and the real networks are compatible with a non-trivial fraction of them. (In ARACNe3-top100 regulons the fraction with |Spearman| < 0.2 is 7.5 % in real cohorts and 9.7 %/37 % in synthetic data; the weak-signal value is dominated by false-positive edges, since the precision of ARACNe3-top100 there is only 0.327, and should not be read as non-monotone content.)

What the evidence does establish is that such edges contribute little to MR scores, and it quantifies how little. The association mode of a non-monotone edge is close to zero, so it enters NaRnEA's directed enrichment with almost no weight and contributes only through the undirected component. We tested this directly with the *true* synthetic network, where no inference error intervenes: with the true regulons, the fraction of true MRs in the top 10 is 0.78 when the non-monotone edges are removed, versus 0.86 when they are present with mode 0 (as any correlation-based method would represent them) and 0.92 when they are given a perfect sign that no correlation can provide (strong signal, n = 400; weak signal: 0.8, 0.86, 0.96). Non-monotone edges thus carry some information for MR scoring (removing them costs 0.05–0.10), an unsigned representation recovers most of it, and the remainder, one to three standard errors, is the upper bound of what a scoring method able to use non-monotone dependencies could gain over the current one. Consistent with this, FSR and ARACNe3-top100 agree on MRs as closely as two versions of ARACNe3 agree with each other (Section 3.2).

### 3.4 Recovery of known master regulators in synthetic data (Fig. S5)

Fraction of the 10 true MRs among the top 10 by |NES| (minimum 30 targets; random 0.05; true-network upper bound 0.8–0.92):

| n | Strong: FSR | Strong: ARACNe3-top100 | Strong: fixed p | Weak: FSR | Weak: ARACNe3-top100 | Weak: fixed p |
|---|---|---|---|---|---|---|
| 30 | 0.38 | 0 | 0 | 0.3 | 0 | 0 |
| 60 | 0.76 | 0 | 0 | 0.74 | 0 | 0 |
| 100 | 0.84 | 0.2 | 0.12 | 0.86 | 0 | 0 |
| 200 | 0.9 | 0.76 | 0.74 | 0.88 | 0.02 | 0.32 |
| 400 | 0.82 | 0.74 | 0.8 | 0.9 | 0.76 | 0.74 |

FSR recovered MRs at every n, within 0.06 of the upper bound from n = 100 (weak signal; strong: 0.06), and well above random but well below the bound at n = 30. ARACNe3 (ARACNe3-top100) and corto at fixed p recovered none at n ≤ 60, and under weak signal none at n ≤ 200 (ARACNe3) or n ≤ 100 (fixed p): their networks were empty or their regulons too small to evaluate. At n = 400 the three schemes lie within 0.08 of each other under strong signal and within 0.16 under weak signal.

**These synthetic comparisons do not inform the choice between FSR and ARACNe3-top100 at n ≥ 95 in real data.** The synthetic ARACNe3 networks are far sparser than the real ones (at n = 100: 0.662 precision on 3,380 edges, against 2,584,773–4,140,622 edges and 1,204–1,317 evaluable regulators in real cohorts with 95 ≤ n ≤ 130), so the simulation reproduces ARACNe3's failure mode at small n but not its behaviour on real tumours with n ≥ 95. The A-vs-B comparison rests on the real-data evidence of Section 3.2, which measures consistency rather than correctness. What the synthetic data do establish is the comparison between FSR and corto at fixed p for n < 95, where the latter was in use.

FSR's edges are imprecise at small n (precision 0.186 at n = 30, 0.903 at n = 400, strong signal). Many false positives at intermediate n are indirect, that is, a regulator correlated with the true one: 33.1 per 100 edges at n = 100 (strong). Edge-level comparisons between schemes must be made at equal edge counts, because at n = 400 FSR's DPI leaves 7,288 edges while ARACNe3-top100 is at its cap of 20,000: when ARACNe3 edges are ranked by support and cut to the same count as FSR, their precision (0.893) and indirect false positives per 100 edges (10.4) are the same as FSR's (0.903 and 9.64). The two methods rank edges comparably; they differ in how many they keep. NaRnEA tolerates imprecise edges better than missing regulons, because it aggregates evidence over 100 weighted targets per regulator.

### 3.5 Sensitivity analyses (Fig. S6)

**Correlation threshold (chosen post hoc).** The first run used r_min = 0.4. On synthetic data it degraded MR recovery at large n under weak signal (n = 400: 0.2 vs 0.9 with 0.2), because once n is large the threshold removes genuine moderate correlations before the top-k selection. r_min was then re-chosen from {0, 0.1, 0.2, 0.3, 0.4} on the same data. Values of 0, 0.1 and 0.2 gave indistinguishable MR recovery (all within one standard error of each other at every n and regime), but lower thresholds reduced edge precision at large n (0.352 at r_min = 0 vs 0.903 at 0.2, n = 400, strong). r_min = 0.2 is therefore the largest value that does not cost MR recovery under weak signal (0.3 costs nothing under strong signal), and not an edge of the grid. On real data the threshold had little effect on agreement with ARACNe3-top100 (median ρ 0.926, 0.935, 0.942 for 0.2, 0.3, 0.4; top-20 overlap 70 % in all three), while 0.4 shrank the regulons of the largest cohorts well below the cap (minimum median 32 targets). At n = 27, |r| ≥ 0.2 corresponds to p ≈ 0.3, and up to n ≈ 130 to p between 0.05 and 0.3: in most cohorts the threshold is effectively inactive, the three lowest values of the grid give identical networks, and the regulon is determined by the winner-takes-all DPI and the top-k selection. The insensitivity to r_min within {0, 0.1, 0.2} on real data is therefore largely trivial; the informative comparison is the synthetic one at large n.

**Number of targets and edge weights.** MR results were robust to k (median ρ of NES against k = 100: 0.981 for k = 50, 0.994 for k = 150) and to using the raw bootstrap likelihood as `aw` (0.963).

### 3.6 Minimum cohort size, from subsampling real cohorts (Fig. S7)

Median agreement between the MRs obtained with a subsampled network and those obtained with the full cohort (same signature), over 5 cohorts × 5 replicates:

| n | Top-20 overlap | Top-50 overlap | ρ of NES | Per-regulator edge Jaccard | Evaluable regulators |
|---|---|---|---|---|---|
| 20 | 40 % | 42 % | 0.694 | 0.0355 | 1,299 |
| 27 | 50 % | 52 % | 0.761 | 0.0561 | 1,299 |
| 40 | 55 % | 58 % | 0.823 | 0.0925 | 1,299 |
| 60 | 65 % | 66 % | 0.881 | 0.136 | 1,299 |
| 80 | 70 % | 70 % | 0.909 | 0.183 | 1,299 |
| 100 | 75 % | 76 % | 0.92 | 0.227 | 1,297 |
| 150 | 80 % | 82 % | 0.954 | 0.342 | 1,266 |

For reference, two methods on the same full cohort agree at 70 % (top-20), two subtypes of the same cancer at 45 %, and the null at 35 %. The smallest n whose median reaches the ARACNe3-top100 agreement is 80, but this point estimate should not be over-read: the top-20 overlap moves in steps of 5 % (one regulator), and the interquartile ranges at n = 60 (60 %–70 %) and n = 80 (60 %–80 %) overlap almost entirely. Values at n = 150 are inflated by the overlap between subset and reference in the smaller cohorts and are shown only for completeness.

Three features of the curve matter for the choice of a minimum cohort size:

- **At n = 20–27 a subsampled network reproduces its own cohort's MRs no better than a different subtype of the same cancer does** (top-20 overlap 40 %–50 %, against 45 % between subtypes and 35 % for the null). Networks from cohorts of that size carry little subtype-specific information.
- **Between n = 60 and n = 80 the subsampled network reaches the agreement that a different inference method achieves on the full cohort** (65 % and 70 % vs 70 %), and the agreement improves only slowly beyond that (n = 100: 75 %). n = 40 is intermediate (55 %).
- **Edge-level stability is low throughout** (median per-regulator Jaccard with the full cohort 0.136 at n = 60, 0.227 at n = 100): the regulons change substantially with the sample while the top MRs do not. This is the direct, real-data counterpart of the synthetic result that NaRnEA tolerates imprecise edges, and it confirms that regulons from small cohorts should not be interpreted edge by edge.

The cohorts differ. THCA BRAF-like, whose full-cohort regulons are the smallest (median 53 targets) and whose reference is therefore itself the least stable, reaches only 60 % at n = 100; the other four reach 65 %–80 %. We adopt n ≥ 60 as the inclusion criterion for subtype-specific networks: it is the lower end of the range where the curve meets the method-to-method agreement, n = 40 falls short (55 %), and n ≤ 27 is no better than a different subtype. Cohorts below the threshold (5 of 34: three with n < 30 and two with n = 43) keep their own networks and are reported as exploratory; Section 3.7 shows that a cancer-wide network does not rescue them.

### 3.7 Cohorts below the threshold: own network or cancer-wide network (Fig. S8)

An alternative for a small cohort is to interrogate its signature with a network inferred from all tumours of its cancer, which has a large n but is not subtype-specific. We tested this on 11 large cohorts reduced to 30 random tumours (3 replicates), taking the network inferred from the full cohort as reference and keeping the full-cohort signature. Three networks were compared: the cohort's own network from the 30 tumours; a cancer-wide network (all tumours of the other subtypes plus the same 30, the situation of a real small subtype); and the network of another subtype of the same cancer. The acceptance criterion was set before the analysis: the cancer-wide network should outperform the own network (paired Wilcoxon test over cohorts, p < 0.05) and the network of another subtype.

Median top-20 overlap with the reference was 60 % for the cancer-wide network, 55 % for the own network and 43 % for the network of another subtype. The cancer-wide network was better than the own network in 7 of 11 cohorts (median difference 0.0833, p = 0.0828) and better than another subtype's network in 10 (p = 0.00506). It therefore does not meet the criterion. It also reduces the differences between subtypes: in LUSC, where three of four subtypes are below the threshold, the top-20 overlap between pairs of subtypes is 25 %–60 % with their own networks and 70 %–85 % when all are interrogated with the cancer-wide network (part of the former spread is sampling noise of small networks). In the five real small cohorts, the two networks share 40 %–75 % of the top-20 MRs.

Neither option yields a reliable subtype-specific network at n ≈ 30: the own network has internal validity (it is inferred from the subtype) but low reliability, and the cancer-wide network has reliable edges but is not specific to the subtype. Small cohorts are analysed with their own networks and labelled exploratory (`results/final/tables/analysis_tier.tsv`).

### 3.8 Subtype-specific versus cancer-wide network for large cohorts (sensitivity analysis)

The design infers one network per subtype. To measure how much this subtype specificity contributes to the TMRs, each of the 29 cohorts with n ≥ 60 was also analysed with an FSR network inferred from all tumours of its cancer (which includes the cohort itself), keeping the cohort's own signature. The top-20 TMRs shared 65 % between the two networks (median; interquartile range 55 %–75 %, range 45 %–80 %; median ρ of NES 0.872), close to the agreement between two inference methods on the same cohort (70 %). The overlap was unrelated to the cohort's share of its cancer (Spearman ρ = 0.144), so it is not explained by the cohort dominating the cancer-wide network. When all subtypes of a cancer were interrogated with the cancer-wide network, the overlap between pairs of subtypes barely changed (45 % with their own networks, 48 % with the common network; 32 pairs): the differences between subtypes are carried mostly by their signatures. Cancer-wide networks yielded fewer evaluable regulators (median 975 vs 1,233), because with more samples more regulators pass the threshold and the per-target DPI concentrates targets on fewer of them.

The subtype-specific network therefore changes the TMR ranking by about as much as the choice of inference method, and this analysis cannot separate subtype-specific rewiring from inference noise. We keep subtype-specific networks as the primary design, since it is the design that allows statements about regulatory differences between subtypes, and report the cancer-wide network as a sensitivity analysis: TMRs found with both networks form the robust core of each cohort's result.

### 3.9 Triangle DPI (sensitivity analysis)

corto applies the DPI per target: in each bootstrap only the regulator with the largest |r| keeps an edge to a given target, which discards independent co-regulators. We implemented the ARACNe rule over regulator–regulator–target triangles as an option of corto (`dpi = "triangle"`; the edge a–c is removed only if another regulator b has both |r(b,c)| and |r(a,b)| larger than |r(a,c)|) and evaluated it in FSR, with the criterion fixed in advance of a gain of at least 0.05 in true MRs recovered in the top 10 in synthetic data without loss of edge precision. It did not meet the criterion: the mean difference was -0.008 (range over regime × n -0.06 to 0.06). At n ≤ 100 the triangle rule kept the same precision and recall with fewer indirect false positives (-6.1 per 100 edges at n = 60, -13.9 at n = 100); at n = 400 edge precision fell by -0.527 (strong signal) and -0.211 (weak), because without the one-winner-per-target rule more weak edges survive before truncation. In 5 real cohorts the untruncated networks were 2.82–4 times larger, the top-20 TMRs agreed 60 %–90 % with the per-target rule, and agreement with ARACNe3-top100 rose from 50 %–75 % to 70 %–85 %, at 31.8–46.5 min per network (8 threads) instead of 5.65–7.98 min. FSR keeps the per-target rule.

### 3.10 Robustness to the choice of k (Fig. S9)

The agreement measures above use the top 20 MRs. Repeating them for 7 values of k from 5 to 50 (5, 10, 15, 20, 30, 40, 50), with 95 % bootstrap confidence intervals over cohorts (over cancers for pairs of subtypes, which share a cancer), gives the same picture (Fig. S9A). Agreement between methods stays at 60 %–70 % and agreement between the own and the cancer-wide network is similar; the ordering methods > subtypes of the same cancer > null holds for 6 of 7 values of k and fails only for k ≤ 5, where a single regulator moves the overlap by 20 %. The confidence interval of method agreement lies above that of between-subtype agreement for 5 of 7 values of k. The intervals of between-subtype agreement and of the null overlap at every k (separated for 0 of 7): with cancers as resampling units, the difference between subtypes of one cancer and the null rests on ten clusters and should not be over-read. Overlap expected by chance is negligible (1.6 % at k = 20). The minimum cohort size of Section 3.6 depends moderately on k: the smallest subsample size at which the median overlap reaches method agreement ranges from 60 to 100 across k (Fig. S9B), consistent with n ≥ 60 as a lower bound. The text keeps k = 20, which lies within the range where the conclusions are stable.

### 3.11 Computing cost

On the 15 cohorts with n ≥ 95, ARACNe3 (100 subnetworks, 18 threads, two runs at a time) took 5.53 h of wall-clock time (99.5 core-hours; up to 1.17 h for one cohort). FSR (100 bootstraps, 8 threads, four or five runs at a time) took 1.83 h (14.6 core-hours). Both were measured under contention on a shared server; per-cohort values are in `network_benchmark_runtime.tsv`.

## 4. Limitations

- **No ground truth for real data.** Agreement between schemes measures consistency, not correctness; CollecTRI is incomplete and largely cell-type agnostic. The synthetic model supplies ground truth but omits features of real tumours (pervasive co-expression modules, tumour purity, copy-number effects), and its ARACNe3 networks are much sparser than real ones.
- **Winner-takes-all DPI.** corto keeps, per target and bootstrap, only the regulator with the largest |r| among all that pass the threshold. Any target with more than one true regulator therefore loses all but the strongest, which caps recall (about 0.823 for linear edges at n = 400 in synthetic data with 1–3 regulators per target); regulators that are strongly co-expressed with each other (paralogues, co-regulated modules) are the case where the loss is systematic. ARACNe's DPI has a related tendency, mitigated by a tolerance parameter. MRs that belong to tightly co-expressed regulator families should be interpreted at the family level. corto's bootstrap likelihood reflects this competition and is not a stability estimate.
- **Comparable size, not comparable quality.** Regulons of 100 targets are produced for every cohort, including the three with n < 30, where edge precision in synthetic data was only 0.101–0.186 and where, on real data, a network from that many samples reproduces its cohort's MRs no better than a different subtype does (Section 3.6). Unless the n ≥ 60 criterion proposed there is adopted, results for STAD EBV, LUSC primitive and LIHC S1 should be treated as exploratory, and the two largest cohorts have regulons below the cap.
- **Threshold chosen on the evaluation data.** r_min was selected after inspecting results (Section 3.5). The choice is insensitive within {0, 0.1, 0.2}, largely because the threshold is inactive at the sample sizes of most cohorts, and the grids for r_min and k were narrow.
- **Non-monotonic dependencies.** Their frequency in these cohorts is unresolved (Section 3.3). Should MR scoring ever move to a method that does not weight targets by a signed correlation, the case for a correlation-based network would need to be revisited.
- **Any co-expression regulon is a proxy.** Published benchmarks find that regulons inferred from co-expression, with either measure, underperform curated or ChIP-seq-derived regulons in perturbation assays; the choice of dependence measure is second-order relative to orthogonal validation of the MRs.

## 5. Conclusion

We use a single correlation-based procedure (FSR: corto on ranks, |r| ≥ 0.2, 100 best-supported targets per regulator) for all 34 cohorts. It yields regulons of comparable size across cohort sizes from 27 to 562, and subsampling real cohorts shows that from n ≈ 60 such a network reproduces the full-cohort MRs as well as a second inference method does; where ARACNe3 can also be run, it reproduces ARACNe3's MR results as closely as alternative filterings of ARACNe3 reproduce each other, and far above a null; it recovers known MRs in synthetic data at every sample size, including sizes at which ARACNe3 yields no network; and it costs a fraction of the computing time. Its main justification, however, is comparability: the method effect on MR rankings is about two-thirds of the network component of the differences between subtypes that the study aims to characterise, and more than half of their total, so inferring networks with different methods for different subtypes would confound those comparisons. MI's capacity to capture non-monotonic dependencies is real, but it is neutralised downstream by the signed association mode that NaRnEA uses, and we found no evidence that it changes the MR results.

## Reproducibility

The analysis is `scripts/network_benchmark/` (steps 00–14b, orchestrated by `run_all.sh`; the technical appendix steps 10–11 are run separately), with parameters in `config/config.yaml: network_benchmark` and the environment in `config/env_network_benchmark.yaml`. Every number in this note is taken from `results/final/tables/network_benchmark_key_numbers.tsv`; figures are `results/final/figures/network_benchmark_fig1`–`fig9`. The CollecTRI file used is archived as `network_benchmark_collectri_edges.tsv`. Software: R 4.5.3, corto 1.3.1, NaRnEA 0.2.0 (commit 36b8fdb), DESeq2 1.50.2, ARACNe3 commit f0fa4fb (califano-lab).
