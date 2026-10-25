# Parameterization Robustness

**Status:** Complete (pancreas dataset), re-run on the pySCEv 0.2 scorer. Reproduce with `python run_sweep.py` after `../pancreas_endocrinogenesis/run_analysis.py` has written the scored cache. Numbers in `summary.json`.

## Results Summary

Tested all combinations of k neighbors (5, 10, 15, 20, 30) and bin counts (8, 12, 18, 24, 36) on the pancreas endocrinogenesis dataset (3,696 cells).

**Per-cell rank correlation across the 25 settings (Spearman, 300 pairs):**
- Min: 0.82
- Median: 0.95
- Mean: 0.93
- Max: 1.00

Higher than with the old absolute-angle scorer (median 0.86), whose bin-edge splitting added setting-dependent noise.

**Cell-type ordering:**
- The low end is fully stable: Ngn3 high EP and Pre-endocrine are the two lowest cell types in all 25 settings.
- The high end is not: the two highest are Epsilon/Ductal in 17 settings, Epsilon/Delta in 5, and other pairs in 3. The rank correlation of cell-type medians against the k=30, bins=8 reference ranges from 0.73 to 1.0 (median 0.90).

This matches the speed control in `../pancreas_endocrinogenesis/`: the separation of the progenitor wave is robust, the ordering among the other cell types is weak.

**Bin count vs k neighbors:** Bin count has less impact than k. Varying bins at fixed k=30 produces tight, parallel curves; varying k at fixed bins=8 shows more spread, especially at low k where each neighborhood has few pairs.

**Normalization:** Normalized scores are raw bits divided by log2(n_bins), a fixed rescaling for a given bin count (perfect rank correlation). Absolute values are comparable across bin counts only in the sense of "fraction of the maximum possible spread".

## Figures

- `rank_correlation_heatmap.png` — 25x25 Spearman rank correlation matrix across all parameter combinations
- `celltype_stability.png` — median entropy per cell type, varying k (left) and bins (right)
- `raw_vs_normalized.png` — raw Shannon entropy vs normalized [0,1] entropy
