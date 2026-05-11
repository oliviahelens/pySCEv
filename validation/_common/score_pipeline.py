"""
Shared scoring pipeline for validation analyses.

Standardizes the scVelo preprocessing + pysce.score_angular_velocity_entropy
recipe used in validation/pancreas_endocrinogenesis/, plus the cohort
assignment used by the enrichment analyses.

Two public entry points:

    score(adata, ...)           run the full scVelo -> pysce pipeline in place
    assign_cohorts(adata, ...)  add a 'pyscev_cohort' obs column with values
                                in {"low", "mid", "high"} per cell type
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

import pysce


@dataclass
class ScoreParams:
    n_top_genes: int = 2000
    n_pcs: int = 30
    n_neighbors: int = 30
    velocity_mode: str = "deterministic"
    n_bins: int = 8
    basis: str = "umap"


def score(adata, params: ScoreParams = ScoreParams()):
    """Run filter+normalize -> moments -> velocity -> velocity_graph ->
    angular_velocity_entropy. Writes adata.obs['angular_velocity_entropy'].

    Matches validation/pancreas_endocrinogenesis/velocity_vs_entropy.py
    so the enrichment analyses see the same scoring distribution that the
    primary pancreas pass produced.
    """
    import scanpy as sc
    import scvelo as scv

    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=params.n_top_genes)
    adata._inplace_subset_var(adata.var["highly_variable"].to_numpy())
    scv.pp.moments(adata, n_pcs=params.n_pcs, n_neighbors=params.n_neighbors)
    scv.tl.velocity(adata, mode=params.velocity_mode)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis=params.basis)
    scv.tl.velocity_confidence(adata)

    pysce.score_angular_velocity_entropy(
        adata,
        basis=params.basis,
        n_neighbors=params.n_neighbors,
        n_bins=params.n_bins,
    )
    return adata


def assign_cohorts(
    adata,
    *,
    score_key: str = "angular_velocity_entropy",
    celltype_key: str = "clusters",
    quantile: float = 0.15,
    min_cells_per_celltype: int = 60,
    out_key: str = "pyscev_cohort",
) -> pd.Series:
    """Tag the bottom and top `quantile` of pySCEv scores within each cell
    type as 'low' / 'high'; everything else 'mid'. NaN scores stay NaN.

    Cohorts are computed *within* cell type so the downstream DE doesn't just
    rediscover the cell-type ranking we already have from the primary pass.
    Cell types with fewer than `min_cells_per_celltype` valid cells are
    dropped (the per-cohort sample would be < ~9 cells at q=0.15, too small
    for a Wilcoxon test).

    Returns the new categorical Series and writes it to adata.obs[out_key].
    """
    if not 0 < quantile < 0.5:
        raise ValueError(f"quantile must be in (0, 0.5), got {quantile}")

    scores = adata.obs[score_key].astype(float)
    celltypes = adata.obs[celltype_key].astype("category")

    cohort = pd.Series(
        pd.Categorical([np.nan] * adata.n_obs, categories=["low", "mid", "high"]),
        index=adata.obs_names,
    )

    for ct in celltypes.cat.categories:
        mask = (celltypes == ct).to_numpy() & scores.notna().to_numpy()
        if mask.sum() < min_cells_per_celltype:
            continue
        s = scores[mask]
        lo = s.quantile(quantile)
        hi = s.quantile(1.0 - quantile)
        cohort.loc[s.index[s <= lo]] = "low"
        cohort.loc[s.index[s >= hi]] = "high"
        cohort.loc[s.index[(s > lo) & (s < hi)]] = "mid"

    adata.obs[out_key] = cohort
    return cohort
