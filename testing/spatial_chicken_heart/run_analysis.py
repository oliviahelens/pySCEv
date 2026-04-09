"""
Spatial angular velocity entropy on chicken heart Visium data.

Runs locally — expects SIRV-imputed chicken heart data (Mantri et al. 2021,
pre-processed by SIRV: https://doi.org/10.5281/zenodo.6798659) to already
exist on disk. Produces the figures for testing/spatial_chicken_heart/.

Usage
-----
    python run_analysis.py --data /path/to/chicken_heart_sirv.h5ad
    python run_analysis.py --data-dir /path/to/sirv_zenodo_download/

The Zenodo archive contains several .h5ad files (one per developmental
stage / tissue section). Either point --data at a single file or point
--data-dir at the directory and the script will pick the first .h5ad.

Output
------
    spatial_entropy_tissue.png       — entropy overlaid on tissue coordinates
    spatial_deviation_tissue.png     — mean angular deviation on the same spots
    entropy_vs_deviation_spatial.png — scatter, spatial context
    spatial_vs_umap_neighbors.png    — entropy computed with spatial vs UMAP
                                       neighbors, side by side
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

import pysce


def load_adata(data: Path | None, data_dir: Path | None):
    if data is not None:
        path = data
    elif data_dir is not None:
        candidates = sorted(data_dir.glob("*.h5ad"))
        if not candidates:
            sys.exit(f"No .h5ad files found in {data_dir}")
        path = candidates[0]
        print(f"[load] picking {path.name} from {data_dir}")
    else:
        sys.exit("Must pass --data FILE or --data-dir DIR")

    print(f"[load] reading {path}")
    adata = sc.read_h5ad(path)
    print(f"[load] {adata.n_obs} spots x {adata.n_vars} genes")
    print(f"[load] layers: {list(adata.layers.keys())}")
    print(f"[load] obsm:   {list(adata.obsm.keys())}")

    # SIRV-imputed output should have spliced/unspliced layers.
    # Accept a few known alternate names.
    layer_aliases = {
        "spliced": ["spliced", "Ms", "S", "s"],
        "unspliced": ["unspliced", "Mu", "U", "u"],
    }
    for canonical, options in layer_aliases.items():
        if canonical in adata.layers:
            continue
        for alt in options:
            if alt in adata.layers:
                print(f"[load] aliasing layer '{alt}' -> '{canonical}'")
                adata.layers[canonical] = adata.layers[alt]
                break
        else:
            sys.exit(
                f"No '{canonical}' layer found (tried {options}). "
                f"Layers present: {list(adata.layers.keys())}"
            )

    # Spatial coordinates — accept common alternate obsm keys
    spatial_aliases = ["spatial", "X_spatial", "X_xy_loc", "xy_loc"]
    if "spatial" not in adata.obsm:
        for alt in spatial_aliases:
            if alt in adata.obsm:
                print(f"[load] aliasing obsm '{alt}' -> 'spatial'")
                adata.obsm["spatial"] = adata.obsm[alt]
                break
        else:
            sys.exit(
                "No spatial coordinates found. Tried "
                f"{spatial_aliases}. Present: {list(adata.obsm.keys())}"
            )

    return adata


def compute_velocity(adata):
    """Run scVelo pipeline and project onto UMAP.

    Pre-computes PCA and neighbors via scanpy so scVelo's moments step
    reuses them — the deprecated auto path in scvelo 0.4+ segfaults on
    macOS ARM via pynndescent.
    """
    import scvelo as scv

    print("[velocity] filter genes")
    scv.pp.filter_genes(adata, min_shared_counts=10)

    # SIRV output is typically already log-normalized. Re-normalizing is
    # harmless (scvelo's normalize_per_cell no-ops when it detects that),
    # but log1p on already-logged data is NOT. Only log if scanpy hasn't
    # recorded a log1p step and X values look like raw/normalized counts.
    already_logged = "log1p" in adata.uns
    if not already_logged:
        max_val = float(adata.X.max())
        if max_val > 50:  # heuristic: raw or cp10k, not logged
            print("[velocity] normalize + log1p")
            scv.pp.normalize_per_cell(adata)
            sc.pp.log1p(adata)
        else:
            print(f"[velocity] X looks log-scaled (max={max_val:.2f}), skipping log1p")
    else:
        print("[velocity] adata.uns['log1p'] present, skipping log1p")

    if adata.n_vars > 2000:
        print("[velocity] HVG selection (2000)")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
        adata._inplace_subset_var(adata.var["highly_variable"].to_numpy())

    print("[velocity] PCA (scanpy)")
    sc.pp.pca(adata, n_comps=30)

    print("[velocity] neighbors (scanpy)")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    print("[velocity] scvelo moments (reusing PCA + neighbors)")
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

    print("[velocity] deterministic velocity")
    # Stochastic mode has a numpy 2.x incompatibility in scvelo 0.3.4
    # (ragged array assignment in leastsq_generalized). Deterministic
    # mode uses a simpler regression path and is sufficient for angular
    # structure analysis.
    scv.tl.velocity(adata, mode="deterministic")

    print("[velocity] velocity graph")
    scv.tl.velocity_graph(adata, n_jobs=1)

    if "X_umap" not in adata.obsm:
        print("[velocity] UMAP (direction basis)")
        sc.tl.umap(adata)

    print("[velocity] projecting velocity onto UMAP")
    scv.tl.velocity_embedding(adata, basis="umap")


def mean_angular_deviation(adata, coords_key: str, n_neighbors: int = 30):
    """Per-cell mean angle between its velocity and its spatial neighbors'.

    Uses the 2D UMAP velocity projection (a fixed direction representation)
    but picks neighbors via the given coordinate key.
    """
    from sklearn.neighbors import NearestNeighbors

    V = np.asarray(adata.obsm["velocity_umap"])
    coords = np.asarray(adata.obsm[coords_key])

    norms = np.linalg.norm(V, axis=1)
    valid = norms > 1e-10
    V_unit = np.zeros_like(V)
    V_unit[valid] = V[valid] / norms[valid, np.newaxis]

    nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, len(coords)))
    nn.fit(coords)
    idx = nn.kneighbors(return_distance=False)[:, 1:]

    out = np.full(len(coords), np.nan)
    for i in range(len(coords)):
        if not valid[i]:
            continue
        nbrs = idx[i][valid[idx[i]]]
        if len(nbrs) == 0:
            continue
        cos = V_unit[nbrs] @ V_unit[i]
        cos = np.clip(cos, -1.0, 1.0)
        out[i] = np.arccos(cos).mean()
    return out


def score_with_spatial_neighbors(adata, n_neighbors=30, n_bins=8):
    """Score angular velocity entropy using physical adjacency as neighbors.

    Tricks score_angular_velocity_entropy by aliasing the spatial
    coordinates into X_spatial and the UMAP velocity projection into
    velocity_spatial, then calls it with basis='spatial'.
    """
    # The scoring function treats X_{basis} as the neighbor space and
    # velocity_{basis} as the direction vectors. We want neighbors by
    # physical location but direction from the 2D UMAP projection.
    adata.obsm["X_spatial_basis"] = np.asarray(adata.obsm["spatial"], dtype=float)
    adata.obsm["velocity_spatial_basis"] = np.asarray(adata.obsm["velocity_umap"])

    pysce.score_angular_velocity_entropy(
        adata,
        basis="spatial_basis",
        n_neighbors=n_neighbors,
        n_bins=n_bins,
        key_added="angular_velocity_entropy_spatial",
    )


def plot_tissue(adata, color_key: str, title: str, out: Path, cmap="viridis"):
    coords = np.asarray(adata.obsm["spatial"])
    values = adata.obs[color_key].to_numpy()

    fig, ax = plt.subplots(figsize=(6, 6))
    finite = np.isfinite(values)
    ax.scatter(
        coords[~finite, 0], coords[~finite, 1],
        s=12, c="lightgray", label="NaN", linewidths=0,
    )
    sc_plot = ax.scatter(
        coords[finite, 0], coords[finite, 1],
        c=values[finite], s=14, cmap=cmap, linewidths=0,
    )
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image-space convention
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.solids.set_alpha(1.0)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_entropy_vs_deviation(adata, out: Path):
    ent = adata.obs["angular_velocity_entropy_spatial"].to_numpy()
    dev = adata.obs["mean_angular_deviation_spatial"].to_numpy()
    mask = np.isfinite(ent) & np.isfinite(dev)
    r = np.corrcoef(ent[mask], dev[mask])[0, 1]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(dev[mask], ent[mask], s=8, alpha=0.5, linewidths=0)
    ax.set_xlabel("Mean angular deviation (spatial neighbors)")
    ax.set_ylabel("Angular velocity entropy (spatial neighbors)")
    ax.set_title(f"Pearson r = {r:.2f}  (n = {mask.sum()})")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out}  (r={r:.3f})")


def plot_spatial_vs_umap(adata, out: Path):
    coords = np.asarray(adata.obsm["spatial"])
    sp = adata.obs["angular_velocity_entropy_spatial"].to_numpy()
    um = adata.obs["angular_velocity_entropy_umap"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, vals, title in [
        (axes[0], sp, "Entropy (spatial neighbors)"),
        (axes[1], um, "Entropy (UMAP neighbors)"),
    ]:
        finite = np.isfinite(vals)
        ax.scatter(
            coords[~finite, 0], coords[~finite, 1],
            s=12, c="lightgray", linewidths=0,
        )
        pc = ax.scatter(
            coords[finite, 0], coords[finite, 1],
            c=vals[finite], s=14, cmap="viridis", vmin=0, vmax=1, linewidths=0,
        )
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title)
        cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
        cb.solids.set_alpha(1.0)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=None, help="Single .h5ad file")
    p.add_argument("--data-dir", type=Path, default=None, help="Directory of .h5ad files")
    p.add_argument("--outdir", type=Path, default=Path(__file__).parent)
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--n-bins", type=int, default=8)
    args = p.parse_args()

    adata = load_adata(args.data, args.data_dir)
    compute_velocity(adata)

    print("[score] angular velocity entropy with spatial neighbors")
    score_with_spatial_neighbors(adata, args.n_neighbors, args.n_bins)
    adata.obs["angular_velocity_entropy_spatial"] = adata.obs[
        "angular_velocity_entropy_spatial_basis"
    ]

    print("[score] angular velocity entropy with UMAP neighbors (for contrast)")
    pysce.score_angular_velocity_entropy(
        adata,
        basis="umap",
        n_neighbors=args.n_neighbors,
        n_bins=args.n_bins,
        key_added="angular_velocity_entropy_umap",
    )

    print("[score] mean angular deviation with spatial neighbors")
    adata.obs["mean_angular_deviation_spatial"] = mean_angular_deviation(
        adata, coords_key="spatial", n_neighbors=args.n_neighbors
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    plot_tissue(
        adata,
        "angular_velocity_entropy_spatial",
        "Angular velocity entropy (spatial neighbors)",
        args.outdir / "spatial_entropy_tissue.png",
    )
    plot_tissue(
        adata,
        "mean_angular_deviation_spatial",
        "Mean angular deviation (spatial neighbors)",
        args.outdir / "spatial_deviation_tissue.png",
        cmap="magma",
    )
    plot_entropy_vs_deviation(adata, args.outdir / "entropy_vs_deviation_spatial.png")
    plot_spatial_vs_umap(adata, args.outdir / "spatial_vs_umap_neighbors.png")

    print("[done]")


if __name__ == "__main__":
    main()
