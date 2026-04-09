"""
Spatial angular velocity entropy on chicken heart Visium data.

Runs locally — expects SIRV-imputed chicken heart data (Mantri et al. 2021,
pre-processed by SIRV: https://doi.org/10.5281/zenodo.6798659) to already
exist on disk. Produces the figures for validation/spatial_chicken_heart/.

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


def plot_spatial_vs_umap_scatter(adata, out: Path):
    sp = adata.obs["angular_velocity_entropy_spatial"].to_numpy()
    um = adata.obs["angular_velocity_entropy_umap"].to_numpy()
    mask = np.isfinite(sp) & np.isfinite(um)
    r = np.corrcoef(sp[mask], um[mask])[0, 1]

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(um[mask], sp[mask], s=8, alpha=0.5, linewidths=0)
    lo = min(sp[mask].min(), um[mask].min())
    hi = max(sp[mask].max(), um[mask].max())
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=1, alpha=0.6, label="y = x")
    ax.set_xlabel("Angular velocity entropy (UMAP neighbors)")
    ax.set_ylabel("Angular velocity entropy (spatial neighbors)")
    ax.set_title(f"Pearson r = {r:.2f}  (n = {mask.sum()})")
    ax.legend(loc="upper left", frameon=False)
    ax.set_aspect("equal")
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


def plot_annotation_overlay(adata, annotation_col: str, out: Path):
    """Two-panel tissue map: annotation labels (left) + entropy (right)."""
    coords = np.asarray(adata.obsm["spatial"])
    labels = adata.obs[annotation_col].astype(str).to_numpy()
    ent = adata.obs["angular_velocity_entropy_spatial"].to_numpy()
    unique_labels = sorted(np.unique(labels).tolist())

    # Pick a categorical colormap with enough colors
    n = len(unique_labels)
    if n <= 10:
        cmap = plt.get_cmap("tab10")
    elif n <= 20:
        cmap = plt.get_cmap("tab20")
    else:
        cmap = plt.get_cmap("gist_ncar")
    label_to_color = {lab: cmap(i % cmap.N) for i, lab in enumerate(unique_labels)}
    point_colors = np.array([label_to_color[l] for l in labels])

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Left: annotation
    axes[0].scatter(
        coords[:, 0], coords[:, 1], c=point_colors, s=14, linewidths=0
    )
    axes[0].set_aspect("equal")
    axes[0].invert_yaxis()
    axes[0].set_xticks([]); axes[0].set_yticks([])
    axes[0].set_title(annotation_col.replace("_", " ").title())
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=label_to_color[l], markersize=8,
                   label=l.replace("\n", " "))
        for l in unique_labels
    ]
    axes[0].legend(
        handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
        fontsize=8, frameon=False,
    )

    # Right: entropy
    finite = np.isfinite(ent)
    axes[1].scatter(
        coords[~finite, 0], coords[~finite, 1],
        s=12, c="lightgray", linewidths=0,
    )
    pc = axes[1].scatter(
        coords[finite, 0], coords[finite, 1],
        c=ent[finite], s=14, cmap="viridis", vmin=0, vmax=1, linewidths=0,
    )
    axes[1].set_aspect("equal")
    axes[1].invert_yaxis()
    axes[1].set_xticks([]); axes[1].set_yticks([])
    axes[1].set_title("Angular velocity entropy (spatial neighbors)")
    cb = fig.colorbar(pc, ax=axes[1], fraction=0.046, pad=0.04)
    cb.solids.set_alpha(1.0)

    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_entropy_by_group(
    adata, group_col: str, out: Path,
    title: str | None = None, min_n: int = 15,
):
    """Box plot of spatial entropy distribution per group, ordered by median.

    Groups with fewer than `min_n` spots are dropped from the plot (their
    medians are noise) but still printed in the stdout ranking with a
    '[skipped: small n]' note.
    """
    ent = adata.obs["angular_velocity_entropy_spatial"].to_numpy()
    groups = adata.obs[group_col].astype(str).to_numpy()
    mask = np.isfinite(ent)
    ent = ent[mask]
    groups = groups[mask]

    unique_groups = np.unique(groups)
    all_values = [ent[groups == g] for g in unique_groups]
    all_medians = [np.median(v) for v in all_values]

    # Filter out small groups for the plot, but remember them for the printout
    plot_groups, plot_values, plot_medians = [], [], []
    skipped = []
    for g, v, m in zip(unique_groups, all_values, all_medians):
        if len(v) >= min_n:
            plot_groups.append(g); plot_values.append(v); plot_medians.append(m)
        else:
            skipped.append((g, len(v), m))

    # Sort plotted groups by median
    order = np.argsort(plot_medians)
    sorted_groups = [plot_groups[i] for i in order]
    sorted_values = [plot_values[i] for i in order]
    sorted_medians = [plot_medians[i] for i in order]

    fig, ax = plt.subplots(figsize=(max(7, 0.6 * len(sorted_groups) + 4), 5))
    bp = ax.boxplot(
        sorted_values, vert=True, patch_artist=True,
        showfliers=False, widths=0.6,
    )
    # Color boxes by median entropy (same cmap as tissue map)
    cmap = plt.get_cmap("viridis")
    for patch, m in zip(bp["boxes"], sorted_medians):
        patch.set_facecolor(cmap(np.clip(m, 0, 1)))
        patch.set_alpha(0.85)
    for median_line in bp["medians"]:
        median_line.set_color("black")

    labels = [g.replace("\n", " ") for g in sorted_groups]
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Angular velocity entropy (spatial neighbors)")
    ax.set_ylim(0, 1)
    ax.axhline(np.median(ent), color="gray", linestyle="--", linewidth=1,
               label=f"overall median = {np.median(ent):.2f}")
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    ax.set_title(title or f"Entropy by {group_col}")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    # Print ranked medians for quick reference in the terminal
    print(f"[plot] wrote {out}")
    print(f"  {group_col} medians (low -> high entropy, n >= {min_n}):")
    for lab, m, v in zip(labels, sorted_medians, sorted_values):
        print(f"    {m:.3f}  n={len(v):4d}  {lab}")
    if skipped:
        print(f"  {group_col} skipped (n < {min_n}):")
        for g, n, m in sorted(skipped, key=lambda x: x[2]):
            print(f"    {m:.3f}  n={n:4d}  {g.replace(chr(10), ' ')}")


# ----------------------------------------------------------------------------
# Permutation control
# ----------------------------------------------------------------------------


def build_spatial_knn_weights(coords, k):
    """Row-standardized k-NN weight matrix over spatial coordinates.

    Returns a scipy.sparse CSR matrix W with W[i, j] = 1/k when j is one of
    i's k nearest spatial neighbors (self excluded), 0 otherwise.
    """
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import csr_matrix

    coords = np.asarray(coords)
    n = len(coords)
    k = min(k, n - 1)

    nn = NearestNeighbors(n_neighbors=k + 1)
    nn.fit(coords)
    idx = nn.kneighbors(return_distance=False)[:, 1:]  # drop self

    rows = np.repeat(np.arange(n), k)
    cols = idx.reshape(-1)
    data = np.full(n * k, 1.0 / k)
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def morans_i(values, W):
    """Moran's I of a scalar field under a precomputed weight matrix.

    NaN entries in `values` are masked out per call; the submatrix is
    re-row-standardized after masking. Returns NaN if fewer than 2 valid
    entries or if the variance is zero.
    """
    from scipy.sparse import diags

    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    if mask.sum() < 2:
        return float("nan")

    W_m = W[mask][:, mask].tocsr()
    row_sum = np.asarray(W_m.sum(axis=1)).ravel()
    # Avoid divide-by-zero for rows whose neighbors all dropped out.
    inv = np.where(row_sum > 0, 1.0 / np.where(row_sum > 0, row_sum, 1.0), 0.0)
    W_m = diags(inv) @ W_m

    x = values[mask] - values[mask].mean()
    denom = float(x @ x)
    if denom == 0.0:
        return float("nan")
    S0 = float(W_m.sum())
    if S0 == 0.0:
        return float("nan")
    num = float(x @ (W_m @ x))
    return (len(x) / S0) * (num / denom)


def kruskal_by_region(values, groups, min_n=15):
    """Kruskal-Wallis H statistic across groups with at least min_n members.

    Returns (H, n_groups). (NaN, 0) if fewer than 2 groups survive the filter.
    """
    from scipy.stats import kruskal

    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    mask = np.isfinite(values)
    values = values[mask]
    groups = groups[mask]

    arrays = []
    for g in np.unique(groups):
        v = values[groups == g]
        if len(v) >= min_n:
            arrays.append(v)

    if len(arrays) < 2:
        return float("nan"), len(arrays)
    return float(kruskal(*arrays).statistic), len(arrays)


def region_medians(values, groups, region_order, min_n=15):
    """Fixed-order vector of per-region medians (NaN if group drops below min_n)."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    mask = np.isfinite(values)
    values = values[mask]
    groups = groups[mask]

    out = np.full(len(region_order), np.nan)
    for i, g in enumerate(region_order):
        v = values[groups == g]
        if len(v) >= min_n:
            out[i] = float(np.median(v))
    return out


def score_entropy_from_velocity(adata, V_perm, n_neighbors, n_bins):
    """Score angular velocity entropy on a permuted velocity field.

    Overwrites adata.obsm['velocity_spatial_basis'] with V_perm, runs
    pysce.score_angular_velocity_entropy under a temporary key, returns the
    resulting score array, and cleans up the temporary keys. Never touches
    adata.obsm['velocity_umap'].
    """
    adata.obsm["velocity_spatial_basis"] = np.asarray(V_perm)
    pysce.score_angular_velocity_entropy(
        adata,
        basis="spatial_basis",
        n_neighbors=n_neighbors,
        n_bins=n_bins,
        key_added="_perm_ent",
    )
    ent = adata.obs["_perm_ent"].to_numpy().copy()
    # Clean up temporary obs / uns entries
    del adata.obs["_perm_ent"]
    adata.uns.pop("_perm_ent_params", None)
    return ent


def plot_null_histogram(null, observed, title, xlabel, out: Path):
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(null, bins=30, color="gray", alpha=0.8, edgecolor="white")
    if np.isfinite(observed):
        ax.axvline(
            observed, color="crimson", linewidth=2,
            label=f"observed = {observed:.3g}",
        )
        ax.legend(loc="upper left", frameon=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("permutations")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {out}")


def plot_tissue_example(adata, observed_ent, shuffled_ent, out: Path):
    """Two-panel tissue map: observed entropy vs one shuffled entropy."""
    coords = np.asarray(adata.obsm["spatial"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, vals, title in [
        (axes[0], observed_ent, "Observed entropy"),
        (axes[1], shuffled_ent, "Shuffled velocities (i=0)"),
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


def run_permutation_test(adata, args):
    """Shuffle velocity vectors across spots and recompute the spatial entropy.

    Records Moran's I and (if Mantri region labels are present) Kruskal-Wallis
    H for each permutation, writes a null-distribution npz, a text summary,
    two histograms, and a side-by-side tissue example.
    """
    import time

    N = int(args.permute)
    base_seed = int(args.permute_seed)
    min_n = 15  # matches plot_entropy_by_group filter

    print(f"[perm] running {N} permutations (seed base = {base_seed})")

    V_observed = np.asarray(adata.obsm["velocity_umap"]).copy()
    V_basis_observed = np.asarray(adata.obsm["velocity_spatial_basis"]).copy()
    observed_ent = adata.obs["angular_velocity_entropy_spatial"].to_numpy().copy()
    coords = np.asarray(adata.obsm["spatial"])
    n = len(V_observed)

    # Build the spatial weight matrix once (reused for every Moran's I call)
    W = build_spatial_knn_weights(coords, args.n_neighbors)

    has_regions = "region" in adata.obs.columns
    if has_regions:
        groups = adata.obs["region"].astype(str).to_numpy()
        finite_obs = np.isfinite(observed_ent)
        # Regions passing the min_n filter on the observed (finite) data
        region_counts = {}
        for g in np.unique(groups[finite_obs]):
            region_counts[g] = int(np.sum(groups[finite_obs] == g))
        region_order = sorted(
            [g for g, c in region_counts.items() if c >= min_n]
        )
        print(f"[perm] {len(region_order)} regions passing n >= {min_n}")
    else:
        groups = None
        region_order = []
        print("[perm] no 'region' column — skipping Kruskal / per-region stats")

    observed_I = morans_i(observed_ent, W)
    if has_regions:
        observed_kw, n_obs_groups = kruskal_by_region(observed_ent, groups, min_n)
        observed_med = region_medians(observed_ent, groups, region_order, min_n)
    else:
        observed_kw, n_obs_groups = float("nan"), 0
        observed_med = np.zeros(0)

    print(f"[perm] observed Moran's I = {observed_I:.4f}")
    if has_regions:
        print(f"[perm] observed Kruskal-Wallis H = {observed_kw:.2f} "
              f"({n_obs_groups} groups)")

    null_I = np.full(N, np.nan)
    null_kw = np.full(N, np.nan)
    null_med = (
        np.full((N, len(region_order)), np.nan) if has_regions else None
    )
    example_shuffled_ent = None

    try:
        loop_start = time.time()
        for i in range(N):
            t0 = time.time()
            rng = np.random.default_rng(base_seed + i)
            perm = rng.permutation(n)
            V_perm = V_observed[perm]

            ent_perm = score_entropy_from_velocity(
                adata, V_perm, args.n_neighbors, args.n_bins
            )

            null_I[i] = morans_i(ent_perm, W)
            if has_regions:
                null_kw[i], _ = kruskal_by_region(ent_perm, groups, min_n)
                null_med[i] = region_medians(
                    ent_perm, groups, region_order, min_n
                )

            if i == 0:
                example_shuffled_ent = ent_perm.copy()

            if i in (0, 9) or (i + 1) % 25 == 0 or i == N - 1:
                dt = time.time() - t0
                elapsed = time.time() - loop_start
                print(
                    f"[perm] {i + 1}/{N}  last={dt:.2f}s  elapsed={elapsed:.1f}s  "
                    f"I={null_I[i]:.3f}"
                    + (f"  KW={null_kw[i]:.1f}" if has_regions else "")
                )
    finally:
        # Restore observed state so subsequent code / interactive inspection
        # sees the real run, not the last shuffle.
        adata.obsm["velocity_spatial_basis"] = V_basis_observed
        adata.obs["angular_velocity_entropy_spatial"] = observed_ent

    # One-sided empirical p-values (observed >= null)
    finite_I = np.isfinite(null_I)
    p_I = (1 + int(np.sum(null_I[finite_I] >= observed_I))) / (1 + int(finite_I.sum()))
    if has_regions:
        finite_kw = np.isfinite(null_kw)
        p_kw = (1 + int(np.sum(null_kw[finite_kw] >= observed_kw))) / (
            1 + int(finite_kw.sum())
        )
    else:
        p_kw = float("nan")

    # Expected value of Moran's I under spatial randomness with row-standardized
    # weights: E[I] = -1 / (n - 1). Print for sanity-checking the null.
    expected_I = -1.0 / (n - 1)
    print(
        f"[perm] null Moran's I mean = {np.nanmean(null_I):.4f}  "
        f"(expected ~{expected_I:.4f})"
    )
    print(f"[perm] p(Moran's I >= observed) = {p_I:.4f}")
    if has_regions:
        print(f"[perm] null KW H mean = {np.nanmean(null_kw):.2f}")
        print(f"[perm] p(KW H >= observed) = {p_kw:.4f}")

    # Save artifacts
    npz_path = args.outdir / "permutation_null.npz"
    save_dict = {
        "observed_morans_i": np.array(observed_I),
        "null_morans_i": null_I,
        "observed_kw_h": np.array(observed_kw),
        "null_kw_h": null_kw,
        "p_morans_i": np.array(p_I),
        "p_kw_h": np.array(p_kw),
        "n_permutations": np.array(N),
        "base_seed": np.array(base_seed),
    }
    if has_regions:
        save_dict["observed_medians"] = observed_med
        save_dict["null_medians"] = null_med
        save_dict["region_order"] = np.array(region_order, dtype=object)
    np.savez(npz_path, **save_dict)
    print(f"[perm] wrote {npz_path}")

    # Text summary
    summary_path = args.outdir / "permutation_summary.txt"
    with open(summary_path, "w") as f:
        f.write("Permutation control for spatial chicken heart entropy\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"n permutations : {N}\n")
        f.write(f"base seed      : {base_seed}\n")
        f.write(f"n spots        : {n}\n")
        f.write(f"k neighbors    : {args.n_neighbors}\n")
        f.write(f"n bins         : {args.n_bins}\n")
        f.write("p-values are one-sided (observed >= null).\n\n")

        f.write("Moran's I\n")
        f.write("-" * 40 + "\n")
        f.write(f"observed : {observed_I:.6f}\n")
        f.write(f"expected (spatial randomness, -1/(n-1)) : {expected_I:.6f}\n")
        if finite_I.any():
            ni = null_I[finite_I]
            f.write(f"null mean   : {ni.mean():.6f}\n")
            f.write(f"null std    : {ni.std():.6f}\n")
            f.write(
                "null 5 / 50 / 95 percentile : "
                f"{np.percentile(ni, 5):.6f}  "
                f"{np.percentile(ni, 50):.6f}  "
                f"{np.percentile(ni, 95):.6f}\n"
            )
            f.write(f"null max    : {ni.max():.6f}\n")
        f.write(f"p-value     : {p_I:.4f}\n\n")

        if has_regions:
            f.write("Kruskal-Wallis H across regions (n >= 15)\n")
            f.write("-" * 40 + "\n")
            f.write(f"observed H  : {observed_kw:.4f}\n")
            f.write(f"n groups    : {n_obs_groups}\n")
            if finite_kw.any():
                nk = null_kw[finite_kw]
                f.write(f"null mean   : {nk.mean():.4f}\n")
                f.write(f"null std    : {nk.std():.4f}\n")
                f.write(
                    "null 5 / 50 / 95 percentile : "
                    f"{np.percentile(nk, 5):.4f}  "
                    f"{np.percentile(nk, 50):.4f}  "
                    f"{np.percentile(nk, 95):.4f}\n"
                )
                f.write(f"null max    : {nk.max():.4f}\n")
            f.write(f"p-value     : {p_kw:.4f}\n\n")

            f.write("Per-region medians (observed vs null mean / 5 / 95)\n")
            f.write("-" * 40 + "\n")
            f.write(
                f"{'region':<40} {'obs':>8} {'null_mean':>10} "
                f"{'null_5':>8} {'null_95':>8}\n"
            )
            for i, g in enumerate(region_order):
                col = null_med[:, i]
                col = col[np.isfinite(col)]
                if col.size > 0:
                    nm = np.mean(col)
                    n5 = np.percentile(col, 5)
                    n95 = np.percentile(col, 95)
                else:
                    nm = n5 = n95 = float("nan")
                label = g.replace("\n", " ")[:40]
                f.write(
                    f"{label:<40} {observed_med[i]:>8.3f} "
                    f"{nm:>10.3f} {n5:>8.3f} {n95:>8.3f}\n"
                )
    print(f"[perm] wrote {summary_path}")

    # Figures
    plot_null_histogram(
        null_I, observed_I,
        title=f"Moran's I null distribution (N={N})",
        xlabel="Moran's I",
        out=args.outdir / "permutation_morans_i.png",
    )
    if has_regions:
        plot_null_histogram(
            null_kw, observed_kw,
            title=f"Kruskal-Wallis H null distribution (N={N})",
            xlabel="Kruskal-Wallis H",
            out=args.outdir / "permutation_kw.png",
        )
    if example_shuffled_ent is not None:
        plot_tissue_example(
            adata, observed_ent, example_shuffled_ent,
            out=args.outdir / "permutation_tissue_example.png",
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=None, help="Single .h5ad file")
    p.add_argument("--data-dir", type=Path, default=None, help="Directory of .h5ad files")
    p.add_argument("--outdir", type=Path, default=Path(__file__).parent)
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--n-bins", type=int, default=8)
    p.add_argument(
        "--permute", type=int, default=0,
        help="Number of velocity-shuffle permutations (0 = skip).",
    )
    p.add_argument(
        "--permute-seed", type=int, default=0,
        help="Base seed; permutation i uses default_rng(base + i).",
    )
    args = p.parse_args()

    adata = load_adata(args.data, args.data_dir)
    compute_velocity(adata)

    print("[score] angular velocity entropy with spatial neighbors")
    score_with_spatial_neighbors(adata, args.n_neighbors, args.n_bins)

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
    plot_spatial_vs_umap_scatter(
        adata, args.outdir / "spatial_vs_umap_entropy_scatter.png"
    )

    # Anatomical region / cell-type overlays if annotations are present.
    # SIRV-imputed Mantri files ship with 'region' and 'celltype_prediction'.
    if "region" in adata.obs.columns:
        plot_annotation_overlay(
            adata, "region", args.outdir / "region_vs_entropy_tissue.png"
        )
        plot_entropy_by_group(
            adata, "region", args.outdir / "entropy_by_region.png",
            title="Angular velocity entropy by anatomical region",
        )
    else:
        print("[plot] skipping region overlay (no 'region' column in obs)")

    if "celltype_prediction" in adata.obs.columns:
        plot_entropy_by_group(
            adata, "celltype_prediction",
            args.outdir / "entropy_by_celltype.png",
            title="Angular velocity entropy by predicted cell type",
        )
    else:
        print("[plot] skipping celltype overlay (no 'celltype_prediction' column)")

    if args.permute > 0:
        run_permutation_test(adata, args)

    print("[done]")


if __name__ == "__main__":
    main()
