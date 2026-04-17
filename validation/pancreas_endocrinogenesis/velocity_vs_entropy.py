"""
Scatter per-cell velocity magnitude against angular velocity entropy on
the scVelo pancreas dataset.

"Velocity values" here means the per-cell velocity length (L2 norm of the
UMAP-projected velocity vector, as also reported by
scv.tl.velocity_confidence as 'velocity_length'). Each point is one cell.

Usage
-----
    python velocity_vs_entropy.py
    python velocity_vs_entropy.py --out velocity_vs_entropy.png

Output
------
    velocity_vs_entropy.png  — per-cell scatter, colored by cell type
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scvelo as scv

import pysce


def prepare(adata):
    import scanpy as sc
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata._inplace_subset_var(adata.var["highly_variable"].to_numpy())
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    # deterministic mode avoids a numpy 2.x incompatibility in scvelo 0.3.4
    # stochastic regression path
    scv.tl.velocity(adata, mode="deterministic")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="umap")
    scv.tl.velocity_confidence(adata)
    return adata


def plot(adata, out_path: Path) -> None:
    V = np.asarray(adata.obsm["velocity_umap"])
    length = np.linalg.norm(V, axis=1)
    entropy = np.asarray(adata.obs["angular_velocity_entropy"])

    finite = np.isfinite(length) & np.isfinite(entropy) & (length > 0)
    r_raw = np.corrcoef(length[finite], entropy[finite])[0, 1]
    r_log = np.corrcoef(np.log10(length[finite]), entropy[finite])[0, 1]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    def _scatter(ax, x, mask_fn):
        if "clusters" in adata.obs:
            cats = adata.obs["clusters"].astype("category")
            palette = plt.get_cmap("tab10")
            for i, cat in enumerate(cats.cat.categories):
                m = (cats == cat).to_numpy() & mask_fn
                ax.scatter(
                    x[m], entropy[m],
                    s=8, alpha=0.6, color=palette(i % 10), label=str(cat),
                )
        else:
            ax.scatter(x[mask_fn], entropy[mask_fn], s=8, alpha=0.6)

    _scatter(axes[0], length, finite)
    axes[0].set_xlabel("velocity magnitude  (||velocity_umap||)")
    axes[0].set_ylabel("angular velocity entropy  (normalized)")
    axes[0].set_title(f"linear x  (r = {r_raw:.2f})")

    _scatter(axes[1], length, finite)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("velocity magnitude  (log scale)")
    axes[1].set_ylabel("angular velocity entropy  (normalized)")
    axes[1].set_title(f"log x  (r with log-length = {r_log:.2f})")
    axes[1].legend(loc="best", fontsize=8, frameon=False, markerscale=1.5)

    fig.suptitle("Pancreas: per-cell velocity values vs velocity entropy")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    print(f"[plot] wrote {out_path}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--out", type=Path,
        default=Path(__file__).parent / "velocity_vs_entropy.png",
    )
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--n-bins", type=int, default=8)
    args = p.parse_args()

    print("[load] scv.datasets.pancreas()")
    adata = scv.datasets.pancreas()
    print(f"[load] {adata.n_obs} cells x {adata.n_vars} genes")

    print("[prep] scVelo preprocessing + velocity")
    prepare(adata)

    print("[score] angular velocity entropy")
    pysce.score_angular_velocity_entropy(
        adata, basis="umap",
        n_neighbors=args.n_neighbors, n_bins=args.n_bins,
    )

    plot(adata, args.out)


if __name__ == "__main__":
    main()
