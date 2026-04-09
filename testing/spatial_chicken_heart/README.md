# Spatial Transcriptomics: Chicken Heart Development

**Status:** Planned — requires downloading external data from Zenodo (https://doi.org/10.5281/zenodo.6798659)

## Dataset

**Chicken Heart Development** (Mantri et al. 2021)
- Paper: https://www.nature.com/articles/s41467-021-21892-z
- Code/data: https://github.com/madhavmantri/chicken_heart
- GEO: GSE149457
- Spatial RNA-seq (Visium) + matched scRNA-seq across 4 developmental stages (days 4, 7, 10, 14). 12 tissue sections, ~22k single-cell transcriptomes.

## Why this dataset

- Active differentiation (epicardial EMT, cardiomyocyte maturation), multiple lineages
- Spatial coordinates grounded in real tissue architecture
- Tests whether angular velocity entropy is more informative when "neighbors" are physically adjacent cells rather than UMAP neighbors

## Velocity on spatial data

Visium only gives total counts (no spliced/unspliced). We use SIRV-imputed data:
- SIRV paper: https://academic.oup.com/nargab/article/6/3/lqae100/7728020
- SIRV repo: https://github.com/tabdelaal/SIRV
- Pre-processed data with imputed spliced/unspliced: https://doi.org/10.5281/zenodo.6798659

SIRV integrates spatial + scRNA-seq into a shared latent space, then transfers spliced/unspliced ratios from k-nearest scRNA-seq neighbors to each spatial spot.

## Analysis plan

1. Download SIRV-processed chicken heart data from Zenodo (already has imputed spliced/unspliced counts)
2. Run scVelo to compute velocity vectors
3. Compute angular velocity entropy using **spatial neighbors** (physical adjacency in tissue) instead of UMAP neighbors
4. Compute mean angular deviation on the same spatial neighbors for comparison
5. Overlay both metrics on tissue coordinates (and H&E images if available)
6. Compare: do entropy and mean angular deviation diverge more in spatial context than they did in UMAP context?

## What a good result looks like

- Entropy overlaid on tissue coordinates highlights biologically interpretable zones (e.g., boundaries between cell types, regions of active EMT)
- Entropy and mean angular deviation diverge more in spatial context than in UMAP context (because physical neighbors can be in genuinely different transcriptional states)
- Bonus: entropy flags regions where SIRV imputation may be unreliable (spurious velocity disagreement from noisy imputation)

## Practical notes

- Start with a single tissue section (e.g., one day-7 slice) to prototype before scaling to all 12
- Dataset should be small enough to run locally without HPC
