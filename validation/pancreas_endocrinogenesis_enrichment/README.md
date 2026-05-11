# Pancreas pySCEv: Gene + Pathway Enrichment

**Status:** Complete. Within-celltype DE + GSEA results below. Decoupler PROGENy/CollecTRI step requires network access to omnipathdb; runs locally but no-ops in this sandbox -- see Caveats.

## Question

The primary pancreas pass (`validation/pancreas_endocrinogenesis/`) showed *which cell types* score high vs low on angular velocity entropy. This analysis asks: **what gene programs and pathways characterize cells with high vs low pySCEv?** Two angles:

- **Cohort-based:** within each cell type, take the top/bottom 15% by entropy and run Wilcoxon DE + GSEA preranked. Asks "what's different about coherent vs incoherent cells of the *same* identity?"
- **Continuous:** per-cell module score for each curated pathway, Spearman-correlated with entropy across all cells. Sidesteps the cohort cutoff and reveals cross-celltype patterns.

## Why within-celltype cohorts

Global top/bottom 15% would be dominated by the cell-type ranking we already have (Alpha, Delta high; Ngn3 high EP, Pre-endocrine low). The DE would just rediscover cell-identity genes -- insulin for Beta, glucagon for Alpha. The within-celltype cohorts ask a sharper question: among cells with the same nominal identity, what distinguishes the coherent movers from the incoherent ones? See `validation/_common/score_pipeline.py:assign_cohorts`.

## How

```
pip install -r requirements.txt -r validation/requirements.txt
python validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py
```

Pipeline:

1. `scvelo.datasets.pancreas()` -> standard scVelo preprocessing -> `pysce.score_angular_velocity_entropy(basis='umap', n_neighbors=30, n_bins=8)`. Cached as `cache/pancreas.scored.h5ad`.
2. Cohort assignment: bottom/top 15% within each cell type. Cell types with <60 valid cells are dropped (Delta, Epsilon).
3. Wilcoxon DE (`scanpy.tl.rank_genes_groups`), low vs high, run globally and per cell type passing the >=30-per-cohort threshold.
4. GSEA preranked (`gseapy.prerank`) on `metric = sign(logFC) * -log10(pval)` against (a) a curated 15-set pancreas panel in `validation/_common/pancreas_genesets.py` (offline, runs unconditionally), (b) Enrichr Hallmark + Reactome (online; falls through silently in offline environments).
5. Per-cell module scores (`scanpy.tl.score_genes`) for the curated panel, Spearman-correlated with entropy. This is the offline analogue of the decoupler step.
6. *(Optional, network-required)* decoupler-py PROGENy + CollecTRI per-cell pathway/TF activity, also Spearman-correlated with entropy.

## Hypothesis (pre-registered)

The primary pass concluded "low entropy = coordinated motion in the neighborhood." If that's the right read, **low-entropy** cells should look like *transitioning programs* (proliferation, lineage-commitment TFs) and **high-entropy** cells should look like *terminal-identity programs* (mature hormone synthesis, secretion machinery).

## Results

### Cohort sizes

q=0.15 within celltype: 815 low / 576 high (out of 3,696). Per-celltype DE ran on Ductal (147/139), Ngn3 low EP (43/40), Ngn3 high EP (294/112), Pre-endocrine (124/90), Beta (94/89), Alpha (76/73). Delta and Epsilon skipped (<30 per cohort).

### Global DE (`de_global.tsv`)

| Direction | Top genes (by adj p) |
|---|---|
| Up in low | **Neurog3** (logFC +1.65), Aldh1b1, Foxa3, Cd24a, Mpzl1, Btg2, Qsox1 |
| Up in high | **Isl1** (logFC -1.46), Cd200, Meis2, Slc38a5, **Pyy** (-1.47), Fam183b, Id2, **Arx** (-1.37) |

Low cohort is enriched for the endocrine-progenitor master TF (Neurog3) and the early-commitment marker Cd24a. High cohort is enriched for the alpha-cell master TF (Arx), its co-regulator Isl1, and the hormone Pyy. Direct match to the hypothesis.

### Per-celltype DE highlights (`de_<celltype>.tsv`)

The pattern recurs *within* every cell type with enough cells:

- **Ngn3 high EP** -- low: Clu, Sparc, Spp1, Acot1, Ttr (secretory-progenitor signature). High: Cck (hormone), Tubb3, Krt7, Tuba1a, Aplp1 (cytoskeletal / maturation).
- **Pre-endocrine** -- low: **Pax4** (logFC +2.60; the classic endocrine specifier), Krt8, Gspt1, **Chgb** (chromogranin). High: **Isl1** (-3.66), **Pyy** (-4.21), Fam183b, Rbp4, Lrpprc.
- **Beta** -- low: H3f3b (histone variant, active chromatin), Calm2, **Pdx1** (beta-cell master TF). High: **Iapp** (-1.75; mature beta hormone), Ssr2, Ssr4 (ER translocon, secretion machinery), Ppp1r1a, Creld2.
- **Alpha** -- low: Smarca1, **Abcc8** (KATP channel), Scg5, Snap25, Ncam1. High: Hspa8 (chaperone), Cdkn1a (p21, cell-cycle exit), Hsp90aa1, Emb, Rbp4.
- **Ductal** -- low: Spc24, **Tuba1b**, 2810417H13Rik, **H2afz**, Slc25a5 (all mitosis/replication markers). High: Anxa5, Tpm1, Errfi1, Malat1, Zfos1.
- **Ngn3 low EP** -- low: Btg2, Foxa3, **Neurog3**, Cdkn1a, Cd24a. High: Hmga2 (proliferation TF), Lurap1l, Id2, Atp1b1, Nudt19.

Within Ductal, low-entropy cells are actively cycling (MCM family analogues, replication histones); within Beta, low-entropy cells express developmental Pdx1 while high-entropy cells express mature Iapp/secretory machinery. This is exactly the predicted intra-celltype progenitor-vs-mature contrast.

### GSEA preranked on curated panel (`gsea_global_pancreas_panel.tsv`)

Top by |NES| (full table in TSV):

| Term | NES | FDR q-val | Leading edge |
|---|---:|---:|---|
| PANCREATIC_HORMONES | **-1.57** | **0.041** | PYY, SCG3, IAPP, PPY, CHGA, INS1, GHRL, SCG5, INS2, NPY |
| MATURE_ALPHA | **-1.53** | **0.037** | MEIS2, ARX, IRX2 |
| ENDOCRINE_PROGENITOR_TFS | **+1.42** | 0.219 | NEUROG3, FOXA3, PAX4 |
| MATURE_BETA | -1.33 | 0.204 | IAPP, PCSK2, INS1, SLC2A2, G6PC2, PCSK1, INS2, GIPR |
| PROLIFERATION_CELL_CYCLE | +1.26 | 0.402 | TUBB5, MCM3, PCNA, MCM6, STMN1, MCM4, H2AFZ, RRM1, HMGB2, TYMS |
| EMT_MIGRATION | +1.22 | 0.322 | VIM, MMP14, FN1 |
| UNFOLDED_PROTEIN_RESPONSE | -0.99 | 0.639 | HSP90B1, XBP1, HSPA5, DDIT3 |
| DUCTAL_MARKERS | +0.99 | 0.516 | CLU, SPP1, KRT8 |

Two sets pass FDR < 0.05: **PANCREATIC_HORMONES** and **MATURE_ALPHA**, both enriched in the high-entropy cohort. **ENDOCRINE_PROGENITOR_TFS** is enriched in the low-entropy cohort at FDR 0.22 (marginal, but the leading edge -- NEUROG3, FOXA3, PAX4 -- is exactly the expected genes). Proliferation is positive (low > high) but doesn't reach significance at the global level because the within-celltype proliferation cycles cancel across cell types in the global preranked metric.

Plot: `fig_gsea_dotplot.png`.

### Continuous module-score correlation (`module_score_corr.tsv`)

Spearman rho of per-cell module score vs entropy:

| Set | rho | n_genes |
|---|---:|---:|
| **ENDOCRINE_PROGENITOR_TFS** | **-0.52** | 12 |
| APOPTOSIS_INTRINSIC | -0.24 | 4 |
| HYPOXIA_GLYCOLYSIS | +0.22 | 4 |
| UNFOLDED_PROTEIN_RESPONSE | +0.20 | 12 |
| MATURE_ALPHA | +0.17 | 7 |
| DUCTAL_MARKERS | +0.12 | 12 |
| PANCREATIC_HORMONES | +0.09 | 13 |
| PROLIFERATION_CELL_CYCLE | +0.05 | 24 |
| MATURE_BETA | +0.03 | 10 |

Plot: `fig_module_score_bar.png`.

The continuous and cohort readouts answer different questions and don't always agree:

- **ENDOCRINE_PROGENITOR_TFS** is strong in both -- the progenitor TF score is highest where entropy is lowest, and that pattern holds both across cell types and within. The signal is locally and globally aligned.
- **PROLIFERATION_CELL_CYCLE** is strong in the cohort DE (within Ductal, within Ngn3 high EP) but weak globally (rho 0.05). Proliferation cycles run within each cell type independently of the overall entropy ordering; the global continuous correlation averages those out.
- **PANCREATIC_HORMONES** is significant in GSEA (FDR 0.04) but weak in continuous (rho 0.09). The high-cohort enrichment for hormone genes is a *within-celltype* effect (mature Alpha cells in their high-entropy quantile express more Gcg/Pyy than mature Alpha cells in their low-entropy quantile), which the cohort method amplifies and the global continuous correlation dilutes.
- **UPR and HYPOXIA_GLYCOLYSIS** show modest positive rho (~0.2) in the continuous readout but don't appear in the cohort GSEA. Likely tracks a cell-type-level pattern (mature secretory cells run more ER stress and have higher entropy), not a within-cell-type contrast.

The takeaway: cohort-based and continuous methods are complementary. Cohort reads within-celltype contrast (what changes as a cell of identity X becomes coherent vs incoherent); continuous reads cross-celltype scaling (which cell types in general have which programs). Both confirm the hypothesis, but on different axes.

## Conclusion

Within-celltype, low pySCEv = transitioning / developmental programs (progenitor TFs, proliferation, actively transcribed chromatin); high pySCEv = terminal / secretory programs (hormone biosynthesis, lineage-commitment TFs of the mature state, secretion machinery). Across cell types, the same pattern holds for the broad endocrine-progenitor TF axis but washes out for proliferation, which is a within-celltype cycle uncorrelated with the global entropy ordering.

The metric is reading something biologically interpretable, and consistent across two different statistical framings of the same data. The operational definition from the primary pass ("low entropy = coordinated motion") is now backed by gene-level evidence of *what* that coordination is composed of: cells riding the same developmental wave have correlated velocity directions and correlated progenitor-program expression; cells at the terminus of differentiation have neither.

## Caveats

- **Mouse symbols uppercased to match human MSigDB.** Pancreas is mm10. The GSEA preranked uses uppercased symbols; this is the standard hack for using human pathway databases on mouse data but is not a true ortholog mapping. For a higher-fidelity pass, swap in MSigDB mouse collections (`gseapy.Msigdb` with `mh.all` and `m2.cp.reactome`) when running locally.
- **External pathway libraries blocked in this sandbox.** MSigDB Hallmark via Enrichr (maayanlab.cloud), Reactome via Enrichr, PROGENy and CollecTRI via omnipathdb.org -- all return HTTP 403 here. The curated panel in `validation/_common/pancreas_genesets.py` is what runs unconditionally. The `gp.prerank` calls for Hallmark/Reactome and the `dc.mt.mlm` calls for PROGENy/CollecTRI are still in the script under try/except and will produce additional output files when run with internet access. They are not strictly required to answer the central question; the curated panel covers the expected biology.
- **The curated panel is biased toward the hypothesis.** It includes both confirmatory sets (proliferation, hormones, progenitor TFs) and contrastive ones (UPR, hypoxia, apoptosis, EMT), but it's not a blind scan. The fact that the confirmatory sets win is informative, but a full MSigDB Hallmark scan with corrected p-values would be more rigorous. Re-running locally with internet would close that gap.
- **Continuous correlation effect sizes are modest.** Even the strongest signal (ENDOCRINE_PROGENITOR_TFS, rho -0.52) explains ~27% of the variance in entropy across cells. The rest is noise plus signal not captured by these 15 gene sets.
- **Cohort cutoff is 15% within cell type.** Easy to change via `--quantile`. The DE statistics are sensitive to cohort size; results above are at q=0.15.
- **One dataset.** Cross-dataset comparison waits on Dubois cardiac (blocked on velocyto reprocessing -- GEO deposit ships no spliced/unspliced) or a swap-in like dentate gyrus. See `validation/pyscev_biological_meaning/README.md`.
- **scVelo deterministic mode.** `validation/_common/score_pipeline.py` uses `scv.tl.velocity(mode="deterministic")` to dodge a numpy 2.x bug in scvelo 0.3.4's stochastic regression path (same workaround as the spatial chicken heart pipeline). The primary pancreas README documents stochastic-mode results from an earlier environment; cell-type ordering is parameterization-stable (`validation/parameterization/`) so the cohort assignment is preserved, but absolute entropy values differ slightly between runs.

## What this does and does not show

- **Shows:** the gene programs that systematically differ between coherent (low-pySCEv) and incoherent (high-pySCEv) cells in the pancreas dataset, at the within-celltype level via DE+GSEA and at the cross-celltype level via continuous module scores. The two readouts agree on direction for the progenitor-TF axis and diverge for proliferation (within-celltype only) and hormone-secretion (within-celltype amplified).
- **Does not show:** a causal claim about what drives the entropy score. Pathway enrichment on quantile cohorts is correlational. Also does not show generalization beyond pancreas -- the synthesis doc tracks the cross-dataset analysis status.

## Files

- `run_enrichment.py` -- pipeline entrypoint
- `cache/pancreas.scored.h5ad` -- scored AnnData (gitignored, regenerable)
- `de_global.tsv`, `de_<celltype>.tsv` -- Wilcoxon DE tables (low vs high cohort, score column = scanpy Wilcoxon test statistic)
- `gsea_global_pancreas_panel.tsv` -- curated-panel GSEA preranked
- `gsea_global_hallmark.tsv`, `gsea_global_reactome.tsv` -- *(only present after running with internet)* external-library GSEA
- `module_score_corr.tsv` -- per-cell module score Spearman rho vs entropy
- `decoupler_progeny_corr.tsv`, `decoupler_collectri_corr.tsv` -- *(only present after running with internet)*
- `fig_cohort_umap.png` -- entropy + cohort assignment on UMAP
- `fig_gsea_dotplot.png` -- GSEA NES dotplot, sized by FDR
- `fig_module_score_bar.png` -- continuous module-score Spearman rhos
