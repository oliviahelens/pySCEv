# Spatial Transcriptomics: Chicken Heart Development

**Status:** First pass complete on day-14 SIRV-imputed Visium section (1967 spots). Follow-ups noted below.

## Dataset

- Chicken Heart Development (Mantri et al. 2021). Paper: https://www.nature.com/articles/s41467-021-21892-z · GEO: GSE149457
- Visium spatial + matched scRNA-seq, days 4/7/10/14, 12 sections, ~22k cells
- Visium has no spliced/unspliced → use SIRV-imputed version from Zenodo: https://doi.org/10.5281/zenodo.6798659 (SIRV paper: https://academic.oup.com/nargab/article/6/3/lqae100/7728020)

## Why

- Tests whether velocity entropy is more informative when "neighbors" are **physically adjacent** spots rather than UMAP neighbors
- Real tissue architecture, active differentiation (epicardial EMT, cardiomyocyte maturation)

## How to run

```bash
KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  python testing/spatial_chicken_heart/run_analysis.py \
  --data /path/to/Visium_D14_adata.h5ad
```

- Env vars avoid a pynndescent/OpenMP segfault on macOS ARM
- Script aliases SIRV's non-standard keys (`X_xy_loc` → `spatial`), does PCA + neighbors via scanpy explicitly (scVelo 0.4's deprecated auto path segfaults), uses deterministic velocity mode (scVelo 0.3.4 stochastic mode has a numpy 2.x bug in `leastsq_generalized`)
- Prototype on one section first

## Results (day 14, 1967 spots, k=30, bins=8)

### Spatial vs UMAP neighbors: r = 0.39

Core payoff. Same entropy metric, same velocity vectors, only the choice of neighbor space changes. ~85% of the variance is unique to one or the other — physical adjacency is capturing something transcriptional adjacency can't see. (`spatial_vs_umap_entropy_scatter.png`, `spatial_vs_umap_neighbors.png`)

### Entropy vs mean angular deviation: r = 0.52

Down from r = 0.80 on pancreas (UMAP neighbors). In spatial context the two metrics diverge more — the scatter has a crescent shape where low deviation forces low entropy but **high deviation can pair with any entropy value**. That's exactly the "unimodal vs multimodal distribution at the same mean angle" regime where entropy's sensitivity to shape matters. ~50% unique variance here vs ~20% on pancreas. (`entropy_vs_deviation_spatial.png`)

### Tissue maps

- `spatial_entropy_tissue.png` — shows regional structure, not noise. Dark (coherent-flow) rim along one edge, broader mid-to-high entropy interior. Smooth contiguous zones, not speckle.
- `spatial_deviation_tissue.png` — mostly low with scattered bright outliers. Does **not** trace the same regional structure the entropy map shows. On this dataset the first-moment metric looks more like outlier detection than tissue zoning; entropy is the more informative of the two.

## Caveats

- **No cell-type / anatomical overlay.** Can say "entropy highlights contiguous tissue zones," can't say "that rim is epicardium undergoing EMT." Needs Mantri's regional annotations to make a biological claim.
- **Velocity direction basis is UMAP-of-Visium.** Visium spots are multi-cellular, UMAP structure is weaker than on scRNA-seq; direction vectors inherit that weakness. Consider redoing with high-dim angular entropy in PCA space.
- **Deterministic mode, not stochastic.** scVelo + numpy 2.x bug. Shouldn't affect angular structure but is a footnote.
- **One section, one developmental stage.** Day 14. Should replicate across the other 11 sections and earlier stages (day 4–10 is where the most active differentiation happens).

## Next steps

- Load Mantri cell-type / region labels and overlay on entropy map → biological interpretation
- Redo with PCA-space velocity (not UMAP-projected) to strengthen the direction signal
- Scale to all 12 sections; look at whether entropy patterns evolve across developmental time
- Sanity check: permute velocity vectors within tissue and confirm the regional pattern collapses
