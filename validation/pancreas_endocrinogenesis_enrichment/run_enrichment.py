"""
Pancreas pySCEv enrichment: what genes/pathways characterize cells with
high vs low angular velocity entropy?

Pipeline:
    1. Load scvelo.datasets.pancreas() and run the standard scoring recipe
       (validation._common.score_pipeline.score) -> adata.obs['angular_velocity_entropy'].
    2. Cache the scored AnnData under cache/pancreas.scored.h5ad.
    3. Define low/high cohorts as bottom/top 15% of pySCEv *within each cell type*
       (validation._common.score_pipeline.assign_cohorts). This avoids the DE
       just rediscovering the cell-type-level ranking we already have from the
       primary pancreas pass.
    4. Wilcoxon DE (low vs high) with scanpy.tl.rank_genes_groups, run once
       globally and once per cell type with enough cells.
    5. GSEA preranked (gseapy.prerank) on the global signed -log10(p)*sign(logFC)
       ranking, against MSigDB Hallmark + Reactome via Enrichr libraries
       (mouse genes uppercased to match human symbols -- a standard hack;
       see README caveats).
    6. Per-cell pathway and TF activity via decoupler-py (PROGENy + CollecTRI,
       both mouse), then Spearman correlation of each pathway/TF score vs
       angular_velocity_entropy across all cells. This sidesteps the cohort
       cutoff entirely.

Outputs (under this folder):
    cache/pancreas.scored.h5ad
    de_global.tsv
    de_<celltype>.tsv  (one per cell type passing min-cell threshold)
    gsea_global_hallmark.tsv
    gsea_global_reactome.tsv
    decoupler_progeny_corr.tsv
    decoupler_collectri_corr.tsv
    fig_cohort_umap.png
    fig_gsea_dotplot.png
    fig_decoupler_progeny_bar.png
    fig_decoupler_collectri_bar.png

Usage:
    python validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py
    python validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py --skip-gsea
    python validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py --skip-decoupler
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO))

from validation._common.score_pipeline import ScoreParams, assign_cohorts, score


CACHE = HERE / "cache"
MIN_CELLS_PER_COHORT = 30  # Wilcoxon DE needs both groups non-trivial


# ---------------------------------------------------------------------------
# Step 1-2: load + score (with cache)
# ---------------------------------------------------------------------------

def load_scored(force: bool = False):
    import anndata as ad
    import scvelo as scv

    cache_path = CACHE / "pancreas.scored.h5ad"
    if cache_path.exists() and not force:
        print(f"[load] cached: {cache_path}")
        return ad.read_h5ad(cache_path)

    print("[load] scv.datasets.pancreas()")
    adata = scv.datasets.pancreas()
    print(f"[load] {adata.n_obs} cells x {adata.n_vars} genes")

    print("[score] scVelo preprocess + pysce.score_angular_velocity_entropy")
    score(adata, ScoreParams())

    CACHE.mkdir(exist_ok=True)
    adata.write_h5ad(cache_path)
    print(f"[load] wrote {cache_path}")
    return adata


# ---------------------------------------------------------------------------
# Step 3: cohort UMAP
# ---------------------------------------------------------------------------

def plot_cohort_umap(adata, out: Path) -> None:
    import scanpy as sc

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sc.pl.umap(
        adata, color="angular_velocity_entropy",
        ax=axes[0], show=False, frameon=False, color_map="viridis",
    )
    axes[0].set_title("angular_velocity_entropy")
    sc.pl.umap(
        adata, color="pyscev_cohort",
        ax=axes[1], show=False, frameon=False,
        palette={"low": "#2166ac", "mid": "#d3d3d3", "high": "#b2182b"},
    )
    axes[1].set_title("pySCEv cohort (within cell type)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[fig] {out}")


# ---------------------------------------------------------------------------
# Step 4: Wilcoxon DE (low vs high)
# ---------------------------------------------------------------------------

def _rank_genes_table(adata, group: str = "low") -> pd.DataFrame:
    """Pull the rank_genes_groups result for one group into a tidy DataFrame."""
    rg = adata.uns["rank_genes_groups"]
    return pd.DataFrame({
        "gene": [r[group] for r in rg["names"]],
        "logfoldchange": [r[group] for r in rg["logfoldchanges"]],
        "pval": [r[group] for r in rg["pvals"]],
        "pval_adj": [r[group] for r in rg["pvals_adj"]],
        "score": [r[group] for r in rg["scores"]],
    })


def run_de(adata, out_dir: Path) -> dict[str, pd.DataFrame]:
    import scanpy as sc

    out_dir.mkdir(exist_ok=True)
    tables: dict[str, pd.DataFrame] = {}

    # Global low vs high
    mask = adata.obs["pyscev_cohort"].isin(["low", "high"])
    sub = adata[mask].copy()
    sc.tl.rank_genes_groups(
        sub, groupby="pyscev_cohort", groups=["low"], reference="high",
        method="wilcoxon",
    )
    tables["global"] = _rank_genes_table(sub, "low")
    tables["global"].to_csv(out_dir / "de_global.tsv", sep="\t", index=False)
    print(f"[de] global: {len(tables['global'])} genes")

    # Per cell type
    for ct in sub.obs["clusters"].cat.categories:
        ct_mask = (sub.obs["clusters"] == ct).to_numpy()
        n_lo = ((sub.obs["pyscev_cohort"] == "low") & ct_mask).sum()
        n_hi = ((sub.obs["pyscev_cohort"] == "high") & ct_mask).sum()
        if n_lo < MIN_CELLS_PER_COHORT or n_hi < MIN_CELLS_PER_COHORT:
            print(f"[de] skip {ct}: low={n_lo}, high={n_hi}")
            continue
        ct_sub = sub[ct_mask].copy()
        sc.tl.rank_genes_groups(
            ct_sub, groupby="pyscev_cohort", groups=["low"], reference="high",
            method="wilcoxon",
        )
        tables[ct] = _rank_genes_table(ct_sub, "low")
        safe = ct.replace(" ", "_").replace("/", "_")
        tables[ct].to_csv(out_dir / f"de_{safe}.tsv", sep="\t", index=False)
        print(f"[de] {ct}: low={n_lo}, high={n_hi}")

    return tables


# ---------------------------------------------------------------------------
# Step 5: GSEA preranked
# ---------------------------------------------------------------------------

def _signed_rank(de: pd.DataFrame) -> pd.DataFrame:
    """Rank metric = sign(logFC) * -log10(pval). Pancreas symbols are mouse;
    we uppercase to match human MSigDB symbols (standard hack -- documented
    in README caveats)."""
    df = de.dropna(subset=["pval", "logfoldchange"]).copy()
    df["pval"] = df["pval"].clip(lower=1e-300)
    df["metric"] = np.sign(df["logfoldchange"]) * -np.log10(df["pval"])
    df["gene"] = df["gene"].str.upper()
    df = df.drop_duplicates("gene").sort_values("metric", ascending=False)
    return df[["gene", "metric"]]


def run_gsea(de_tables: dict[str, pd.DataFrame], out_dir: Path) -> dict:
    import gseapy as gp

    out_dir.mkdir(exist_ok=True)
    rnk = _signed_rank(de_tables["global"])
    results = {}
    for label, lib in [
        ("hallmark", "MSigDB_Hallmark_2020"),
        ("reactome", "Reactome_2022"),
    ]:
        try:
            res = gp.prerank(
                rnk=rnk, gene_sets=lib, outdir=None, threads=4,
                min_size=15, max_size=500, permutation_num=1000, seed=0,
            )
            df = res.res2d.sort_values("NES", key=abs, ascending=False)
            df.to_csv(out_dir / f"gsea_global_{label}.tsv", sep="\t", index=False)
            results[label] = df
            print(f"[gsea] {label}: {len(df)} terms")
        except Exception as exc:
            print(f"[gsea] {label} FAILED: {exc!r}")
    return results


def plot_gsea_dotplot(gsea_results: dict, out: Path, top_n: int = 12) -> None:
    if not gsea_results:
        return
    import seaborn as sns

    rows = []
    for lib, df in gsea_results.items():
        d = df.copy()
        d["library"] = lib
        # Top up- and down-regulated by NES
        d_sig = d[d["FDR q-val"] < 0.25]
        if d_sig.empty:
            d_sig = d
        d_sig = pd.concat([
            d_sig.nlargest(top_n // 2, "NES"),
            d_sig.nsmallest(top_n // 2, "NES"),
        ])
        rows.append(d_sig)
    plot_df = pd.concat(rows)
    plot_df["term_short"] = plot_df["Term"].str.slice(0, 50)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(plot_df))))
    sns.scatterplot(
        data=plot_df, x="NES", y="term_short", hue="library",
        size="FDR q-val", sizes=(200, 30), ax=ax,
    )
    ax.axvline(0, color="grey", linewidth=0.8)
    ax.set_title("GSEA preranked: low pySCEv (NES > 0) vs high pySCEv (NES < 0)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[fig] {out}")


# ---------------------------------------------------------------------------
# Step 6: decoupler-py PROGENy + CollecTRI
# ---------------------------------------------------------------------------

def run_decoupler(adata, out_dir: Path) -> dict:
    import decoupler as dc
    from scipy.stats import spearmanr

    out_dir.mkdir(exist_ok=True)
    out: dict[str, pd.DataFrame] = {}
    entropy = adata.obs["angular_velocity_entropy"].astype(float).to_numpy()
    finite = np.isfinite(entropy)

    for label, fetch in [
        ("progeny", lambda: dc.get_progeny(organism="mouse", top=500)),
        ("collectri", lambda: dc.get_collectri(organism="mouse")),
    ]:
        try:
            net = fetch()
            dc.run_mlm(
                mat=adata, net=net, source="source", target="target",
                weight="weight", verbose=False, use_raw=False,
            )
            est = adata.obsm["mlm_estimate"]
            corrs = []
            for src in est.columns:
                vals = est[src].to_numpy()
                m = finite & np.isfinite(vals)
                if m.sum() < 50:
                    continue
                rho, p = spearmanr(vals[m], entropy[m])
                corrs.append((src, rho, p, m.sum()))
            df = pd.DataFrame(corrs, columns=["source", "spearman_rho", "pval", "n"])
            df = df.sort_values("spearman_rho", key=abs, ascending=False)
            df.to_csv(out_dir / f"decoupler_{label}_corr.tsv", sep="\t", index=False)
            out[label] = df
            print(f"[decoupler] {label}: {len(df)} sources scored")
        except Exception as exc:
            print(f"[decoupler] {label} FAILED: {exc!r}")
    return out


def plot_decoupler_bars(decoupler_results: dict, out_dir: Path, top_n: int = 15) -> None:
    import seaborn as sns

    for label, df in decoupler_results.items():
        d = df.head(top_n).copy()
        if d.empty:
            continue
        d["sig"] = d["pval"].clip(lower=1e-300).apply(lambda p: -np.log10(p))
        fig, ax = plt.subplots(figsize=(7, max(3, 0.3 * len(d))))
        colors = ["#2166ac" if r < 0 else "#b2182b" for r in d["spearman_rho"]]
        sns.barplot(data=d, y="source", x="spearman_rho", palette=colors, ax=ax)
        ax.axvline(0, color="grey", linewidth=0.8)
        ax.set_xlabel(f"Spearman rho ({label} score vs angular_velocity_entropy)")
        ax.set_title(f"{label}: top {top_n} sources by |rho|")
        fig.tight_layout()
        out_path = out_dir / f"fig_decoupler_{label}_bar.png"
        fig.savefig(out_path, dpi=160)
        plt.close(fig)
        print(f"[fig] {out_path}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--force-rescore", action="store_true",
                   help="ignore cache/pancreas.scored.h5ad and recompute")
    p.add_argument("--skip-gsea", action="store_true")
    p.add_argument("--skip-decoupler", action="store_true")
    p.add_argument("--quantile", type=float, default=0.15,
                   help="cohort cutoff (top/bottom q within each cell type)")
    args = p.parse_args()

    adata = load_scored(force=args.force_rescore)

    cohort = assign_cohorts(adata, quantile=args.quantile)
    n_lo = (cohort == "low").sum()
    n_hi = (cohort == "high").sum()
    print(f"[cohort] q={args.quantile}: low={n_lo}, high={n_hi}")

    plot_cohort_umap(adata, HERE / "fig_cohort_umap.png")

    de_tables = run_de(adata, HERE)

    if not args.skip_gsea:
        gsea_results = run_gsea(de_tables, HERE)
        plot_gsea_dotplot(gsea_results, HERE / "fig_gsea_dotplot.png")

    if not args.skip_decoupler:
        decoupler_results = run_decoupler(adata, HERE)
        plot_decoupler_bars(decoupler_results, HERE)


if __name__ == "__main__":
    main()
