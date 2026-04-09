# Pancreas Endocrinogenesis

**Status:** First complete real-data pass. UMAP-basis, UMAP neighbors.

## Dataset

- scVelo built-in `pancreas` dataset (Bastidas-Ponce et al. 2019, Development)
- 3,696 cells, spliced/unspliced available
- Cell types: Ductal, Ngn3 low EP, Ngn3 high EP, Pre-endocrine, Alpha, Beta, Delta, Epsilon

## Why

First test on real scRNA-seq with a well-characterized differentiation hierarchy. The purpose is to see whether angular velocity entropy picks up any structure at all, and whether that structure is (a) interpretable and (b) non-redundant with metrics that already exist.

## How

Standard scVelo preprocessing (`filter_and_normalize`, `moments`, `velocity`, `velocity_graph`) on default parameters, then `pysce.score_angular_velocity_entropy(adata, basis='umap', n_neighbors=30, n_bins=8)`.

## Results (k=30, bins=8, UMAP basis)

### Cell-type ranking

Lowest entropy: **Ngn3 high EP, Pre-endocrine** — the actively transitioning populations along the endocrine specification wave. Highest entropy: **Alpha, Delta** — terminal/mature endocrine populations. Ductal and Ngn3 low EP sit in between.

**Read:** low entropy tracks coherent velocity flow. The initial guess that "low entropy = committed" was wrong in the directional sense — committed/terminal populations here score *higher*, because once cells are no longer moving along a shared trajectory the local direction field becomes disordered. The right statement is "low entropy = coordinated motion in the neighborhood."

### Relationship to other metrics

- **Expression entropy (pySCE original):** r ~ 0.02. Effectively independent. Velocity entropy is not a stand-in for potency.
- **scVelo velocity confidence:** r ~ -0.48. Anti-correlated (as expected — confidence rewards local coherence, as does low entropy), but only ~23% shared variance. The two are not interchangeable.
- **Mean angular deviation (first-moment alternative):** r ~ 0.80. Strongly correlated. Entropy adds sensitivity to distribution shape (multimodal vs unimodal at the same mean angle), and that shape signal is where the ~20% unique variance sits.

## Figures

- `umap_panels.png` — UMAP colored by cell type + angular velocity entropy side-by-side
- `velocity_streamlines.png` — scVelo streamlines over the entropy field
- `entropy_by_celltype.png` — box plot of entropy per cell type
- `entropy_vs_confidence.png` — scatter vs scVelo velocity confidence
- `expression_vs_velocity_entropy.png` — UMAPs comparing expression entropy vs velocity entropy
- `expression_vs_velocity_scatter.png` — per-cell scatter of the two entropies
- `entropy_vs_angular_deviation.png` — entropy vs mean angular deviation scatter
- `alpha_phase_portraits.png` — phase portraits for Alpha-driver genes (sanity check on velocity)

## Caveats

- **One dataset, one parameter set.** Stability across k and bins is covered in `validation/parameterization/`.
- **UMAP-projected velocity.** Direction vectors are 2D projections; any pathology of the UMAP projection flows into the angles. Higher-dim PCA basis is a reasonable next thing to try.
- **Default scVelo stochastic mode.** No numpy 2.x issue on this run (unlike the chicken heart pipeline).
- The correlations above are from a single run; not bootstrapped.

## What this does and does not show

- Shows: on a well-studied differentiation dataset, the metric produces a cell-type ordering that is biologically readable, and it is not redundant with expression entropy or scVelo confidence.
- Does not show: that this ordering is useful for any specific biological question (TIC vs normal stem cell, disease stratification, etc.), or that it generalizes beyond this dataset. See `spatial_chicken_heart/` for one follow-up on different tissue.
