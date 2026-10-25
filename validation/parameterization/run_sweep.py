"""
k-neighbors x bin-count sweep of angular velocity entropy on pancreas.

Reads the scored AnnData written by
validation/pancreas_endocrinogenesis/run_analysis.py (velocity already
projected), rescores on a 5x5 grid, and writes:

    rank_correlation_heatmap.png   Spearman rho between every pair of settings
    celltype_stability.png         celltype medians, varying k and varying bins
    raw_vs_normalized.png          raw bits vs normalized score
    summary.json                   numbers quoted in README.md
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad

import pysce

HERE = Path(__file__).resolve().parent
SCORED = HERE.parent / "pancreas_endocrinogenesis_enrichment" / "cache" / "pancreas.scored.h5ad"
KS = [5, 10, 15, 20, 30]
BINS = [8, 12, 18, 24, 36]
ORDER = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine",
         "Beta", "Alpha", "Delta", "Epsilon"]


def main():
    adata = ad.read_h5ad(SCORED)
    ct = adata.obs["clusters"].astype(str).to_numpy()

    scores = {}
    for k, b in itertools.product(KS, BINS):
        pysce.score_angular_velocity_entropy(
            adata, basis="umap", n_neighbors=k, n_bins=b, key_added="_s")
        scores[(k, b)] = adata.obs["_s"].to_numpy().copy()

    labels = [f"k={k},b={b}" for k, b in scores]
    df = pd.DataFrame({l: v for l, v in zip(labels, scores.values())})
    rho = df.corr(method="spearman").to_numpy()
    off = rho[np.triu_indices_from(rho, 1)]

    fig, ax = plt.subplots(figsize=(10, 8.5))
    im = ax.imshow(rho, vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    fig.colorbar(im, label="Spearman rho")
    ax.set_title("Rank correlation across k x bins settings")
    fig.tight_layout()
    fig.savefig(HERE / "rank_correlation_heatmap.png", dpi=150)
    plt.close(fig)

    def med(v):
        return [np.nanmedian(v[ct == g]) for g in ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    for g_i, g in enumerate(ORDER):
        axes[0].plot(KS, [med(scores[(k, 8)])[g_i] for k in KS], marker="o", label=g)
        axes[1].plot(BINS, [med(scores[(30, b)])[g_i] for b in BINS], marker="o", label=g)
    axes[0].set_xlabel("k neighbors (bins = 8)")
    axes[1].set_xlabel("bins (k = 30)")
    axes[0].set_ylabel("median entropy")
    axes[1].legend(fontsize=8, frameon=False, bbox_to_anchor=(1, 1))
    fig.suptitle("Cell-type medians across settings")
    fig.tight_layout()
    fig.savefig(HERE / "celltype_stability.png", dpi=150)
    plt.close(fig)

    pysce.score_angular_velocity_entropy(adata, basis="umap", normalize=False, key_added="_raw")
    raw = adata.obs["_raw"].to_numpy()
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(raw, scores[(30, 8)], s=5, alpha=0.5)
    ax.set_xlabel("raw entropy (bits)")
    ax.set_ylabel("normalized entropy (/ log2 8)")
    ax.set_title("Normalization is a fixed rescaling")
    fig.tight_layout()
    fig.savefig(HERE / "raw_vs_normalized.png", dpi=150)
    plt.close(fig)

    order_by_setting = {
        l: [ORDER[i] for i in np.argsort(med(v))] for l, v in zip(labels, scores.values())
    }
    ref = np.array(med(scores[(30, 8)]))
    celltype_rank_rho = [
        float(pd.Series(med(v)).corr(pd.Series(ref), method="spearman"))
        for v in scores.values()
    ]
    summary = {
        "spearman_offdiag": {"min": float(off.min()), "median": float(np.median(off)),
                             "mean": float(off.mean()), "max": float(off.max())},
        "celltype_median_rank_rho_vs_k30_b8": {
            "min": float(np.min(celltype_rank_rho)),
            "median": float(np.median(celltype_rank_rho))},
        "lowest_two_by_setting": {l: o[:2] for l, o in order_by_setting.items()},
        "highest_two_by_setting": {l: o[-2:] for l, o in order_by_setting.items()},
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["spearman_offdiag"], indent=2))
    print(json.dumps(summary["celltype_median_rank_rho_vs_k30_b8"], indent=2))


if __name__ == "__main__":
    main()
