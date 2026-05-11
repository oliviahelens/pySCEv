# Pancreas pySCEv: Gene + Pathway Enrichment

**Status:** Pipeline complete; results pending local run.

## Question

The primary pancreas pass (`validation/pancreas_endocrinogenesis/`) showed *which cell types* score high vs low on angular velocity entropy. This analysis asks: **what gene programs and pathways characterize cells with high vs low pySCEv?** Two angles:

- **Cohort-based:** within each cell type, take the top/bottom 15% by entropy and run differential expression + GSEA. Asks "what's different about coherent vs incoherent cells of the same identity?"
- **Continuous:** per-cell PROGENy pathway and CollecTRI TF activity scores, Spearman-correlated with entropy across all cells. Sidesteps the cohort cutoff entirely.

## Why within-celltype cohorts

Global top/bottom 15% would be dominated by the cell-type ranking we already have (Alpha, Delta cells score high; Ngn3 high EP and Pre-endocrine score low). The DE would just rediscover cell-identity genes -- insulin for Beta, glucagon for Alpha, etc. -- which tells us nothing about what entropy itself is reading. Within-celltype cohorts ask a sharper question: among cells with the same nominal identity, what distinguishes the coherent movers from the incoherent ones?

The cohort assignment is in `validation/_common/score_pipeline.py:assign_cohorts` and only runs on cell types with >=60 valid cells (so each cohort has >=9 cells before the DE filter).

## How

Run:

```
pip install -r requirements.txt -r validation/requirements.txt
python validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py
```

Pipeline:

1. `scvelo.datasets.pancreas()` -> standard scVelo preprocessing -> `pysce.score_angular_velocity_entropy(basis='umap', n_neighbors=30, n_bins=8)`. Same recipe as the primary pass; cached as `cache/pancreas.scored.h5ad`.
2. Cohort assignment: bottom/top 15% of entropy *within each cell type*.
3. Wilcoxon DE (`scanpy.tl.rank_genes_groups`), low vs high, run globally and per cell type that has >=30 cells per cohort. Tables: `de_global.tsv`, `de_<celltype>.tsv`.
4. GSEA preranked (`gseapy.prerank`) on the global ranking using `metric = sign(logFC) * -log10(pval)`, against MSigDB Hallmark and Reactome (Enrichr libraries `MSigDB_Hallmark_2020`, `Reactome_2022`). Tables: `gsea_global_<library>.tsv`. Plot: `fig_gsea_dotplot.png`.
5. decoupler-py PROGENy (pathway activity) + CollecTRI (TF activity) scored per cell with multivariate linear model (`dc.run_mlm`), then Spearman-correlated with entropy. Tables: `decoupler_<source>_corr.tsv`. Plots: `fig_decoupler_<source>_bar.png`.

## Hypothesis (to keep honest)

The primary pass concluded "low entropy = coordinated motion in the neighborhood." If that's right, **low-entropy** cells should look like *transitioning programs* -- proliferation, lineage-commitment TFs (Ngn3 targets, Neurod1, Pax4/Pax6) -- because coherent neighborhood motion is the signature of a population riding a shared trajectory. **High-entropy** cells should look like *terminal-identity programs* -- mature hormone synthesis (insulin, glucagon, somatostatin), secretion machinery, mature endocrine TFs -- because once a cell stops moving along a trajectory, the local direction field is dominated by noise rather than coordinated motion.

If we instead see proliferation up at high entropy, or hormone synthesis up at low entropy, that's a real surprise about what the metric is reading.

## Results

*To be filled in after first run; expected here:*

- Top 10 globally up- and down-regulated genes (pancreas-specific commentary)
- Top 5 GSEA Hallmark pathways each direction with NES + FDR
- Top 5 PROGENy pathways and CollecTRI TFs by |Spearman rho| (continuous readout)
- Comparison: do the cohort-based and continuous methods agree on direction?

## Caveats

- **Mouse symbols upper-cased to match human MSigDB.** Pancreas is mm10. The Enrichr Hallmark/Reactome libraries are human. Most metabolic and endocrine gene names are conserved when uppercased (Ins1/Ins2 -> INS1/INS2 collapse onto the human INS family in many libraries; same for GCG, SST, NEUROD1, etc.), but this is not a true ortholog mapping. For a higher-fidelity pass, use `gseapy.Msigdb` to fetch the mouse-native collections (`mh.all` for hallmark, `m2.cp.reactome` for Reactome) and skip the symbol mangling.
- **Continuous correlation is rank-based; effect sizes are small.** Per-cell PROGENy/CollecTRI scores are noisy; a Spearman rho of 0.1 across 3,000 cells can be highly significant but biologically modest. Read direction over magnitude.
- **Cohort cutoff is 15% within cell type.** Easy to swap (`--quantile 0.10` or `0.20`) but the DE statistics are sensitive to cohort size. Re-run for different quantiles before reading too much into a marginal call.
- **One dataset.** Cross-dataset comparison (Q4 in the project plan) needs a second dataset; the original Dubois target is blocked on velocyto reprocessing (no spliced/unspliced in the GEO deposit). See `validation/pyscev_biological_meaning/README.md` for the synthesis status.
- **Wilcoxon DE on log-normalized counts** is the default scanpy recipe; not the most powerful test here, but the one that matches every other validation script in this repo.

## What this does and does not show

- **Shows:** the gene programs and pathways that systematically differ between coherent (low-pySCEv) and incoherent (high-pySCEv) cells *within the same cell type* in the pancreas dataset, by two complementary methods.
- **Does not show:** a causal claim about what drives the entropy score, or generalization beyond pancreas. Pathway enrichment on quantile cohorts is a correlational read on a 2D-projected velocity field, not a mechanistic statement about cell biology.

## Files

- `run_enrichment.py` -- pipeline entrypoint
- `cache/pancreas.scored.h5ad` -- scored AnnData (gitignored if large; regenerable)
- `de_global.tsv`, `de_<celltype>.tsv` -- Wilcoxon DE tables
- `gsea_global_<library>.tsv` -- GSEA preranked tables
- `decoupler_<source>_corr.tsv` -- PROGENy / CollecTRI Spearman correlations
- `fig_cohort_umap.png` -- entropy + cohort assignment on UMAP
- `fig_gsea_dotplot.png` -- top GSEA terms by NES, sized by FDR
- `fig_decoupler_<source>_bar.png` -- top sources by |Spearman rho|
