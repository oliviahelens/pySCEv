# pySCEv Biological Meaning -- Synthesis (in progress)

**Status:** Pancreas enrichment (Q1) complete -- see `validation/pancreas_endocrinogenesis_enrichment/`. Cross-dataset synthesis blocked on Dubois cardiac-development reprocessing on Minerva (GEO deposit GSE205950 has no spliced/unspliced -- confirmed by inspecting `GSM6235981_cc_card.rds.gz`, only `RNA` and `SCT` assays present).

This document is the cross-dataset synthesis -- it does no compute itself. Once Q2/Q3 (Dubois cell-type stratification + enrichment) or a Q5 replacement dataset is done, fill in the convergence/divergence sections below.

## Operational definition (from primary pass)

From `validation/pancreas_endocrinogenesis/`: low angular velocity entropy = coordinated motion in the cell's neighborhood; high entropy = scattered/incoherent local direction field. Terminal/mature populations score *higher* (Alpha, Delta) because once a cell stops moving along a shared trajectory, the local 2D direction field is dominated by noise. Actively transitioning populations (Ngn3 high EP, Pre-endocrine) score *lower*.

This is the operational read. Everything below either confirms, refines, or contradicts it.

## Pancreas enrichment summary (Q1)

Within-celltype cohort DE + GSEA (full results in `validation/pancreas_endocrinogenesis_enrichment/`):

- **Low pySCEv** is enriched for endocrine-progenitor TFs (Neurog3, Pax4, Foxa3) globally and for proliferation markers (MCM family, H2afz, PCNA, Tuba1b) within cell types like Ductal and Ngn3 high EP. Pdx1 also comes up in the Beta low cohort -- developmental beta TF expressed in cells still on the trajectory.
- **High pySCEv** is enriched for mature hormone biosynthesis (Pyy, Iapp, Ins, Gcg, Chga) and mature alpha identity (Isl1, Arx, Meis2). Two sets pass FDR<0.05 in the curated-panel GSEA: PANCREATIC_HORMONES and MATURE_ALPHA, both up in high.
- **Two readouts, different questions.** Cohort DE captures within-celltype contrast (does a Beta cell become Iapp+ as it loses coherent velocity?). Continuous module-score correlation captures cross-celltype scaling (does the progenitor TF score rank cell types in the same order as entropy?). They agree on the progenitor TF axis (rho -0.52) and disagree on proliferation (strong within celltype, ~0 globally).

Bottom line for the synthesis: the pancreas data is consistent with the "low entropy = cell still riding a coordinated trajectory" reading, and the gene-level signal is biologically interpretable. Now we need a second dataset to claim generality.

## Convergent signals across datasets

*To fill after Q3.* Pathways and TFs that move in the **same** direction in both pancreas (endoderm -> endocrine) and Dubois (mesoderm -> cardiac) are the strongest candidates for what the metric is reading in general, since the two datasets share no obvious lineage-specific biology.

Plotting plan: scatter of `NES_pancreas` vs `NES_dubois` for the union of GSEA-significant pathways (FDR < 0.25 in either). Diagonal = convergent; off-axis = dataset-specific.

## Divergent signals

*To fill after Q3.* Pathways significant in one dataset only, or moving opposite directions. Most plausible explanations to consider before claiming biology:

- Tissue / lineage-stage differences (early mesoderm specification vs late endocrine commitment)
- Velocity mode (deterministic vs stochastic vs dynamical)
- UMAP-projection geometry (different distortions in different datasets)
- Cohort size (Dubois has fewer cells per type than pancreas, after the cardiac subset)

## What this still doesn't show

- Generalization beyond two endo/mesoderm datasets at single-time-resolved differentiation. Q5 candidates broaden the test.
- Causal mechanism. Pathway enrichment is correlation, not "this pathway makes cells coherent."
- Single parameter set per dataset (k=30, bins=8, UMAP basis). Robustness is covered for pancreas in `validation/parameterization/`; Dubois will need its own check.
- Whether the metric adds value over `scvelo.tl.velocity_confidence` for any *specific* downstream task. The primary pass showed the two are anti-correlated at r ~ -0.48 (only ~23% shared variance), but "non-redundant" is not the same as "useful."

---

## Q5 -- Other public datasets to validate against

Selection criteria: (a) public scRNA-seq with spliced/unspliced or precomputed velocity, (b) non-trivial differentiation hierarchy, (c) author-curated cell-type annotations, (d) ideally a different germ layer or organ system from pancreas + cardiac so we're not just testing on developmental endo/mesoderm.

### Recommended (in priority order)

1. **Dentate gyrus neurogenesis** -- `scvelo.datasets.dentategyrus()` (Hochgerner et al. 2018). Mouse adult neurogenesis, ~3000 cells, native spliced/unspliced, well-annotated trajectory (RGL -> nIPC -> Neuroblast -> immature/mature granule). Different germ layer (ectoderm) and adult-stem-cell biology rather than embryonic lineage commitment. Biggest test of generality. **Estimated effort: hours, not days. No new data needed.** Same script with a one-line dataset swap.

2. **Gastrulation erythroid lineage** -- `scvelo.datasets.gastrulation_erythroid()` (Pijuan-Sala et al. 2019, subset). Mouse early hematopoiesis, ~10k cells, native spliced/unspliced. Very different temporal scale (mouse E6.5-E8.5 in vivo) from pancreas (postnatal in vitro/ex vivo) and cardiac. Useful complement -- if convergent pathways from pancreas+cardiac also appear here, that's strong evidence the metric is reading something general about coordinated differentiation rather than tissue-specific regulators.

3. **Setty CD34+ bone marrow** -- the Palantir benchmark (Setty et al. 2019). Human, ~5780 CD34+ cells, hematopoietic differentiation hierarchy (HSC -> MEP/CMP/CLP -> downstream). Standard velocity benchmark in the field, so existing comparisons are easy. Not in scvelo.datasets but the loom is on the Palantir Github.

### Considered, deprioritized

- **La Manno developing brain** (the original RNA velocity paper) -- great historical baseline but the dataset is large and the annotations are coarse; analysis effort is higher per unit of new biological signal.
- **Spatial-velocity dataset** to cross-link with `validation/spatial_chicken_heart/` -- worth doing eventually but currently no clean off-the-shelf candidate that has all four criteria. Defer until the spatial pipeline matures.
