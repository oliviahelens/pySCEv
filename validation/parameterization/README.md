# Parameterization Robustness

**Status:** Complete (pancreas dataset)

## Results Summary

Tested all combinations of k neighbors (5, 10, 15, 20, 30) and bin counts (8, 12, 18, 24, 36) on the pancreas endocrinogenesis dataset (3,696 cells).

**Rank correlation across settings:**
- Min: 0.70 (comparing extreme settings like k=5,b=36 vs k=30,b=8)
- Median: 0.86
- Mean: 0.86
- Max: 0.99 (adjacent settings)

**Cell-type ordering:** Stable across all tested settings. Ngn3 high EP and Pre-endocrine consistently rank lowest (most coherent flow), Alpha and Delta consistently rank highest (most disordered). The relative ordering of cell types is preserved even when absolute entropy values shift.

**Bin count vs k neighbors:** Bin count has less impact than k. Varying bins at fixed k=30 produces tight, parallel curves. Varying k at fixed bins=8 shows more spread, especially at low k where the neighborhood is small and noisy.

**Normalization:** Raw and normalized entropy are perfectly rank-correlated (Spearman r = 1.0) since normalization is a monotonic scaling. Normalization only matters for interpretability of absolute values across different embedding dimensions.

## Figures

- `rank_correlation_heatmap.png` — 25x25 Spearman rank correlation matrix across all parameter combinations
- `celltype_stability.png` — Median entropy per cell type, varying k (left) and bins (right)
- `raw_vs_normalized.png` — Raw Shannon entropy vs normalized [0,1] entropy
