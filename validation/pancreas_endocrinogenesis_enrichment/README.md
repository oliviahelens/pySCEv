# Pancreas pySCEv: Gene + Pathway Enrichment

**Status:** Complete; re-run on the pySCEv 0.2 scorer. Within-celltype DE + GSEA results below. Decoupler PROGENy/CollecTRI step requires network access to omnipathdb; runs locally but no-ops in this sandbox -- see Caveats.

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

1. `scvelo.datasets.pancreas()` -> standard scVelo preprocessing -> `pysce.score_angular_velocity_entropy(basis='umap', n_neighbors=30, n_bins=8)`. Cached as `cache/pancreas.scored.h5ad` (also written by `../pancreas_endocrinogenesis/run_analysis.py`).
2. Cohort assignment: bottom/top 15% within each cell type. Cell types with <60 valid cells are dropped (Delta, Epsilon).
3. Wilcoxon DE (`scanpy.tl.rank_genes_groups`), low vs high, run globally and per cell type passing the >=30-per-cohort threshold.
4. GSEA preranked (`gseapy.prerank`) on `metric = sign(logFC) * -log10(pval)` against (a) a curated 15-set pancreas panel in `validation/_common/pancreas_genesets.py` (offline, runs unconditionally), (b) Enrichr Hallmark + Reactome (online; falls through silently in offline environments).
5. Per-cell module scores (`scanpy.tl.score_genes`) for the curated panel, Spearman-correlated with entropy. This is the offline analogue of the decoupler step.
6. *(Optional, network-required)* decoupler-py PROGENy + CollecTRI per-cell pathway/TF activity, also Spearman-correlated with entropy.

## Hypothesis (pre-registered)

The primary pass concluded "low entropy = coordinated motion in the neighborhood." If that's the right read, **low-entropy** cells should look like *transitioning programs* (proliferation, lineage-commitment TFs) and **high-entropy** cells should look like *terminal-identity programs* (mature hormone synthesis, secretion machinery).

## Results

Re-run on the pySCEv 0.2 scorer (pairwise angles, rotation-invariant). The previous version of this file was produced with the old absolute-angle scorer; where conclusions changed it is called out below.

### Cohort sizes

q=0.15 within celltype: 904 low / 559 high (out of 3,696). Per-celltype DE ran on Ductal (138/138), Ngn3 low EP (40/40), Ngn3 high EP (418/97), Pre-endocrine (109/89), Beta (93/89), Alpha (73/73). Delta and Epsilon skipped (<30 per cohort).

Ngn3 high EP's low cohort is oversized because 65% of those cells score exactly 0 (every neighbor pair within one 22.5° bin), so the 15% cutoff falls inside a tie.

### Global DE (`de_global.tsv`)

| Direction | Top genes (by adj p) |
|---|---|
| Up in low | **Neurog3** (logFC +2.21), Gadd45a, Btbd17, Btg2, Marcksl1, Smarcd2, Foxa3 |
| Up in high | Akr1c19, **Isl1** (-1.68), Id2, Atp1b1, Rap1b, Peg3, **Meis2** |

Same axis as before: the low cohort carries the endocrine-progenitor TF Neurog3 and its targets; the high cohort carries the mature-endocrine TFs Isl1 and Meis2. The global cohorts pool the per-celltype cohorts, so they still partly reflect cell-type composition (Ngn3 high EP dominates the low cohort).

### Per-celltype DE highlights (`de_<celltype>.tsv`)

- **Ductal** (124 genes at FDR < 0.05) -- low: Tuba1b, Spc24, H2afz, Hmgb2, **Birc5**, **Pcna** (replication / mitosis). High: Malat1, Anxa5, Rbms3, Tpm1, Errfi1. *Unchanged from the old scorer.*
- **Ngn3 high EP** (166) -- low: Clu, Sparc, Ttr, Spp1, Acot1, Mgst1. High: Cck, Tubb3, Aplp1, Igfbpl1. *Unchanged.*
- **Pre-endocrine** (130) -- low: **Nkx6-1**, Chgb, **Pax4**, Chga, Hpcal1. High: **Isl1**, **Pyy**, **Arx**, Fam183b, Cdkn1a, Scgn. *Unchanged in direction: low = beta-lineage specifiers, high = alpha-lineage / maturing.*
- **Beta** (14) -- low: **Iapp**, **Ins2**, Mapt, Nnat, Hspa5. High: Glud1, **Arx**, Fev. ***Reversed.*** The old scorer gave Pdx1 in low and Iapp/secretory machinery in high. Now the mature beta hormones are in the *low* cohort, and the high cohort carries Arx/Fev, which mark alpha-lineage or immature endocrine cells near the Beta/Alpha boundary.
- **Alpha** (25) -- low: **Gcg**, Pcsk2, **Iapp**, Ttr, Resp18. High: Hspa8, **Cdkn1a**, Tm4sf4, Usp18. ***Reversed for hormones.*** Glucagon processing (Gcg, Pcsk2) is now up in *low*-entropy Alpha cells.
- **Ngn3 low EP** (3) -- low: Neurog3, Cd24a. High: Lurap1l. Too few hits to read.

Within Ductal, Ngn3 high EP and Pre-endocrine the old within-celltype story holds: low-entropy cells are cycling or carry early specification TFs. Within the terminal types (Alpha, Beta) it does **not**: low-entropy cells there express *more* mature hormone, not less. The earlier claim that "high entropy = terminal/secretory within every cell type" depended on the old binning.

### GSEA preranked on curated panel (`gsea_global_pancreas_panel.tsv`)

| Term | NES | FDR q-val |
|---|---:|---:|
| PANCREATIC_HORMONES | **-1.57** | **0.035** |
| MATURE_ALPHA | -1.46 | 0.088 |
| ENDOCRINE_PROGENITOR_TFS | +1.39 | 0.34 |
| DUCTAL_MARKERS | -1.03 | 0.74 |
| EMT_MIGRATION | +0.96 | 1.00 |
| UNFOLDED_PROTEIN_RESPONSE | -0.89 | 0.80 |
| PROLIFERATION_CELL_CYCLE | +0.80 | 0.97 |
| MATURE_BETA | -0.75 | 0.82 |

Only **PANCREATIC_HORMONES** passes FDR < 0.05 (up in the global high cohort); MATURE_ALPHA drops to 0.09. Given the per-celltype reversal above, this global hormone enrichment reflects which cell types fill the pooled high cohort, not a within-identity maturation effect.

Plot: `fig_gsea_dotplot.png`.

### Continuous module-score correlation (`module_score_corr.tsv`)

Spearman rho of per-cell module score vs entropy:

| Set | rho | n_genes |
|---|---:|---:|
| **ENDOCRINE_PROGENITOR_TFS** | **-0.56** | 12 |
| HYPOXIA_GLYCOLYSIS | +0.27 | 4 |
| APOPTOSIS_INTRINSIC | -0.23 | 4 |
| MATURE_ALPHA | +0.20 | 7 |
| UNFOLDED_PROTEIN_RESPONSE | +0.19 | 12 |
| DUCTAL_MARKERS | +0.17 | 12 |
| PROLIFERATION_CELL_CYCLE | +0.06 | 24 |
| PANCREATIC_HORMONES | +0.05 | 13 |
| MATURE_BETA | -0.03 | 10 |

Plot: `fig_module_score_bar.png`.

- **ENDOCRINE_PROGENITOR_TFS** remains the strongest and most consistent signal (rho -0.56; was -0.52), and it agrees with the global and per-celltype DE.
- **PROLIFERATION_CELL_CYCLE** is strong within Ductal but ~0 globally, as before.
- **PANCREATIC_HORMONES / MATURE_BETA** are ~0 globally, consistent with the reversed within-celltype result: there is no simple "more hormone = higher entropy" relationship.

## Conclusion

One finding survives the scorer fix cleanly: **low pySCEv tracks the endocrine-progenitor TF program** (Neurog3, Pax4, Nkx6-1, Foxa3), across cell types (rho -0.56) and within the progenitor cell types. Low entropy also marks cycling cells within Ductal.

The second half of the old conclusion -- "high pySCEv = terminal / secretory program within each identity" -- does **not** survive. In Alpha and Beta, mature hormone genes are now higher in the *low*-entropy cohort. The safer reading is that low entropy marks cells on a coherent, well-resolved flow (the progenitor wave, or the core of a terminal cluster with a consistent residual velocity), and high entropy marks boundary and slow regions where the projected direction field is poorly determined. Velocity speed explains part of that (see `../pancreas_endocrinogenesis/speed_control.png`).

## Caveats

- **Mouse symbols uppercased to match human MSigDB.** Pancreas is mm10. The GSEA preranked uses uppercased symbols; this is the standard hack for using human pathway databases on mouse data but is not a true ortholog mapping. For a higher-fidelity pass, swap in MSigDB mouse collections (`gseapy.Msigdb` with `mh.all` and `m2.cp.reactome`) when running locally.
- **External pathway libraries blocked in this sandbox.** MSigDB Hallmark via Enrichr (maayanlab.cloud), Reactome via Enrichr, PROGENy and CollecTRI via omnipathdb.org -- all return HTTP 403 here. The curated panel in `validation/_common/pancreas_genesets.py` is what runs unconditionally. The `gp.prerank` calls for Hallmark/Reactome and the `dc.mt.mlm` calls for PROGENy/CollecTRI are still in the script under try/except and will produce additional output files when run with internet access. They are not strictly required to answer the central question; the curated panel covers the expected biology.
- **The curated panel is biased toward the hypothesis.** It includes both confirmatory sets (proliferation, hormones, progenitor TFs) and contrastive ones (UPR, hypoxia, apoptosis, EMT), but it's not a blind scan. The fact that the confirmatory sets win is informative, but a full MSigDB Hallmark scan with corrected p-values would be more rigorous. Re-running locally with internet would close that gap.
- **Continuous correlation effect sizes are modest.** Even the strongest signal (ENDOCRINE_PROGENITOR_TFS, rho -0.56) explains ~31% of the variance in entropy across cells. The rest is noise plus signal not captured by these 15 gene sets.
- **Ties at zero.** 16% of all cells (65% of Ngn3 high EP) score exactly 0, so quantile cutoffs inside that cell type are tie-dominated.
- **Speed is not controlled in the DE.** Cohorts are cut on raw entropy, which correlates with log speed (rho -0.49). A speed-matched cohort design is the obvious next refinement.
- **Cohort cutoff is 15% within cell type.** Easy to change via `--quantile`. The DE statistics are sensitive to cohort size; results above are at q=0.15.
- **One dataset.** Cross-dataset comparison waits on Dubois cardiac (blocked on velocyto reprocessing -- GEO deposit ships no spliced/unspliced) or a swap-in like dentate gyrus. See `validation/pyscev_biological_meaning/README.md`.

## What this does and does not show

- **Shows:** low pySCEv is robustly associated with the endocrine-progenitor TF program, globally and within progenitor cell types, and with proliferation within Ductal cells.
- **Does not support:** a general "high entropy = mature/secretory" reading; within Alpha and Beta the hormone genes go the other way.
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
