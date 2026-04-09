# Spatial Transcriptomics: Chicken Heart Development

**Status:** Script ready — run locally, then commit the figures back. The SIRV-processed data is on Zenodo (https://doi.org/10.5281/zenodo.6798659) and the sandbox this repo was developed in has no internet access, so the analysis script is designed to be executed on a machine that already has the data.

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

## How to run

`run_analysis.py` implements steps 2–6 end-to-end. On a machine that has pySCEv installed and the Zenodo data on disk:

```bash
# point at a single h5ad (e.g. one tissue section to prototype)
python testing/spatial_chicken_heart/run_analysis.py \
    --data ~/Downloads/chicken_heart_sirv/day7_section1.h5ad

# or point at the whole Zenodo download directory and it will pick the first .h5ad
python testing/spatial_chicken_heart/run_analysis.py \
    --data-dir ~/Downloads/chicken_heart_sirv/
```

The script expects an AnnData object with:

- `adata.layers['spliced']` and `adata.layers['unspliced']` (SIRV output)
- `adata.obsm['spatial']` — 2D tissue coordinates

It will compute velocity with scVelo, run `pysce.score_angular_velocity_entropy` twice (once with **spatial** neighbors, once with **UMAP** neighbors for contrast), compute mean angular deviation on the spatial neighbors, and write four figures into this directory:

- `spatial_entropy_tissue.png` — entropy overlaid on tissue coordinates
- `spatial_deviation_tissue.png` — mean angular deviation on the same spots
- `entropy_vs_deviation_spatial.png` — scatter + correlation in the spatial context
- `spatial_vs_umap_neighbors.png` — side-by-side of spatial vs UMAP neighbor choice

Prototype on a single day-7 section before scaling to all 12 tissue sections.

## What a good result looks like

- Entropy overlaid on tissue coordinates highlights biologically interpretable zones (e.g., boundaries between cell types, regions of active EMT)
- Entropy and mean angular deviation diverge more in spatial context than in UMAP context (because physical neighbors can be in genuinely different transcriptional states)
- Bonus: entropy flags regions where SIRV imputation may be unreliable (spurious velocity disagreement from noisy imputation)

## Practical notes

- Start with a single tissue section (e.g., one day-7 slice) to prototype before scaling to all 12
- Dataset should be small enough to run locally without HPC
