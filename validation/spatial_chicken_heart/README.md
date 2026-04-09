# Spatial Transcriptomics: Chicken Heart Development

**Status:** First pass complete on day-14 SIRV-imputed Visium section (1967 spots), including Mantri region + cell-type overlays. Follow-ups noted below.

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
  python validation/spatial_chicken_heart/run_analysis.py \
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

- `spatial_entropy_tissue.png` — shows dark (coherent-flow) zones along the ventricular wall and broader mid-to-high entropy interior. **Caveat: the raw "smooth contiguous patches" impression is partly mechanical.** The scoring function uses each spot's 30 nearest physical neighbors, so adjacent spots share ~28/30 of the same neighbor set and the resulting entropy field is already spatially smoothed regardless of biology. Under a random-velocity permutation control (below), the scoring kernel alone produces Moran's I ≈ 0.47 out of the observed 0.74. The biologically meaningful claim is not "patches exist" — it's "the patches align with anatomy," which is what the region ranking and the Kruskal–Wallis test below actually show.
- `spatial_deviation_tissue.png` — mostly low with scattered bright outliers. Does **not** trace the same regional structure the entropy map shows. On this dataset the first-moment metric looks more like outlier detection than tissue zoning; entropy is the more informative of the two.

### Anatomical regions (Mantri labels, n=1967)

Ranked by median entropy (low → high):

| median | n   | region                                   |
|--------|-----|------------------------------------------|
| 0.505  | 188 | Trabecular LV + endocardium              |
| 0.598  | 693 | Compact LV + inter-ventricular septum    |
| 0.649  | 286 | Right ventricle                          |
| 0.688  | 130 | Epicardium                               |
| 0.717  | 454 | Atria                                    |
| 0.745  | 216 | Valves                                   |

**Read:** the lowest-entropy regions are the **left ventricular wall** — both the trabecular + endocardial layer and the compact myocardium + septum. These are the mechanically committed cardiomyocyte zones at day 14, where local neighborhoods are dominated by a single coordinated differentiation/maturation program → coherent velocity flow → low entropy. The right ventricle is intermediate. The **highest-entropy regions are valves and atria**, both known developmental mixing zones: valves are the site of endocardial-to-mesenchymal transition and fibroblast infiltration, atria contain mixed myocardial + conductive populations. Physically adjacent spots there are heading in different transcriptional directions.

Note: **epicardium sits at the median**, not low. An earlier guess that "low entropy = EMT" was wrong; the correct statement appears to be "low entropy = coordinated differentiation wave in committed tissue; high entropy = multiple divergent trajectories in physical proximity." This is consistent with the pancreas result (Ngn3-high EP, the coordinated wave, was lowest there too).

### Cell types (Mantri predictions, only n ≥ 15 groups)

| median | n   | cell type                 |
|--------|-----|---------------------------|
| 0.524  | 101 | TMSB4X high cells         |
| 0.598  | 203 | Immature myocardial cells |
| 0.621  | 241 | Erythrocytes              |
| 0.623  |  46 | Mural cells               |
| 0.635  | 977 | Cardiomyocytes-1          |
| 0.714  |  28 | Endocardial cells         |
| 0.741  |  77 | Vascular endothelial cells|
| 0.751  |  20 | Valve cells               |
| 0.755  | 264 | Fibroblast cells          |

(Skipped for tiny n: Cardiomyocytes-2 n=1, Macrophages n=3, Epi-epithelial n=6.)

**Consistent with the region story.** Lowest: TMSB4X-high cells (thymosin β4+, a cardiac progenitor / regeneration marker — actively maturing population) and immature myocardial cells. Highest: fibroblast and valve cells (the divergent-trajectory populations). Erythrocytes are interesting at mid-low — probably noise, they have little meaningful velocity signal.

### Permutation control (N=100)

Shuffle velocity vectors across spots (vectors stay intact; only their assignment to spots is randomized), recompute the spatial-neighbor entropy, repeat 100× with fixed seeds. Measures how much of the observed structure survives when the spatial-to-velocity link is broken. Run with `--permute 100`.

| statistic          |   observed | null (mean ± std) | z-score | p      |
|--------------------|-----------:|------------------:|--------:|-------:|
| Moran's I          |     0.7415 |   0.4691 ± 0.0407 |   +6.70 | 0.0099 |
| Kruskal–Wallis H   |   329.20   |   65.76 ± 45.88   |   +5.74 | 0.0099 |

Both p-values are at the N=100 floor (1/101); the z-scores are the real signal. Moran's I is computed on a row-standardized spatial k-NN graph (k=30); Kruskal–Wallis is across the 6 Mantri regions with n≥15.

**Read:**

- **The region ranking is solid biology.** Observed KW H is ~5× the null mean and +5.74 standard deviations above it. Under random velocity→spot assignment, the 6 regions do pick up *some* apparent separation (null mean 66, not 0) because the k-NN-smoothed entropy field has enough contiguity that contiguous region labels partially align with it by chance. But the observed 329 is nowhere near that. The LV-low / valves-atria-high ordering holds up to the control by a large margin. This is the cleanest biological claim on this dataset.

- **The tissue patches are partly biology, mostly method.** Observed Moran's I is +6.70σ above null — unambiguously above, but the null is already at 0.47 out of the observed 0.74. That 0.47 is the "free" autocorrelation produced by the k-NN scoring kernel alone, regardless of what velocities it's fed. Roughly two-thirds of the raw autocorrelation in the tissue map is mechanical; the remaining ~third (plus the specific alignment with anatomy) is biology. The visual "look at those big smooth zones" framing is overselling. The correct statement: patches are real and they line up with anatomy, but the smoothness itself is largely the scoring kernel's footprint.

- **Why the two tests differ in strength:** Kruskal–Wallis directly tests whether a *specific* partition (the named Mantri regions) separates medians. The null k-NN-smoothed field has contiguity but no systematic alignment with any particular label set, so KW is cleanly powered to detect "this partition is meaningful." Moran's I just asks "is the field more clustered than chance" — and the null is already very clustered, so the signal has to fight through a high baseline.

Artifacts:

- `permutation_morans_i.png`, `permutation_kw.png` — null histograms with observed value marked
- `permutation_tissue_example.png` — observed entropy next to one shuffled (i=0) run. The shuffled panel is the critical eyeball test: it still looks patchy (because of the null ≈ 0.47), but the patches don't align with the LV wall / atria / valves.
- `permutation_summary.txt` — numeric summary with per-region null percentiles
- `permutation_null.npz` — raw arrays for re-plotting / further analysis

## Caveats

- **Velocity direction basis is UMAP-of-Visium.** Visium spots are multi-cellular, UMAP structure is weaker than on scRNA-seq; direction vectors inherit that weakness. Consider redoing with high-dim angular entropy in PCA space.
- **Deterministic mode, not stochastic.** scVelo + numpy 2.x bug. Shouldn't affect angular structure but is a footnote.
- **One section, one developmental stage.** Day 14. Should replicate across the other 11 sections and earlier stages (day 4–10 is where the most active differentiation happens).

## Next steps

- Redo with PCA-space velocity (not UMAP-projected) to strengthen the direction signal
- Scale to earlier developmental stages (day 4, 7, 10) where the epicardium is actively undergoing EMT; check whether the epicardium stratifies lower than at day 14 when EMT is largely complete
- Pairwise Mann–Whitney between regions to rank specific contrasts (the overall Kruskal–Wallis is already well above the permutation null; now ask which pairs drive it)
