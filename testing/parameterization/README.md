# Parameterization Robustness

**Status:** Not yet started

## Goal

Test whether angular velocity entropy results are stable across parameter choices, or if they're sensitive to specific settings. A robust metric should give consistent biological conclusions across reasonable parameter ranges.

## Parameter sweep

| Parameter | Values to test | Default |
|---|---|---|
| `n_neighbors` | 5, 10, 15, 20, 30 | 30 |
| `n_bins` | 8, 12, 18, 24, 36 | 8 |
| `normalize` | True, False | True |

## Analysis plan

1. Run angular velocity entropy across all parameter combinations on at least two datasets:
   - Pancreas endocrinogenesis (scRNA-seq, already validated)
   - Chicken heart spatial (once spatial testing is set up)
2. For each pair of parameter settings, compute **Spearman rank correlation** between the per-cell scores — this tests whether the relative ordering of cells is preserved even if absolute values shift
3. Visualize as a heatmap: parameter setting vs parameter setting, colored by rank correlation
4. Also check: do cell-type-level conclusions (which populations are high/low entropy) change across settings?

## What a good result looks like

- Rank correlations > 0.90 across reasonable parameter ranges (the relative ordering of cells is stable)
- Cell-type-level conclusions (e.g., "Ngn3 high EP has lowest entropy") hold across all tested settings
- If certain parameter ranges produce instability, document the thresholds and add guidance to the README

## What a bad result looks like

- Rank correlations below 0.80 between adjacent parameter values
- Cell-type ordering flips depending on bin count or neighbor count
- Would suggest the metric is measuring noise rather than biology
