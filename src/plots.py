"""Write the four required figures (design doc 5.2, 5.4, 5.5).

Usage: python -m src.plots   (requires `make de` and `make compare` first)
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")  # headless: works in CI and over SSH

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from src.common import FDR, FIGURES_DIR, MEMORY_GENES_CSV, SUMMARY_JSON
from src.compare import join_results, load_de_results
from src.load import load_counts, load_edger, load_meta

DPI = 200
EXPOSURE_COLORS: dict[str, str] = {"NSD": "#4C72B0", "SD": "#DD8452"}
BATCH_MARKERS: dict[str, str] = {"batch1": "o", "batch2": "s"}


def _save(fig: plt.Figure, name: str) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES_DIR / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote figures/{name}")


def plot_pca(counts: pd.DataFrame, meta: pd.DataFrame, n_top: int = 500) -> None:
    """PCA on log2(CPM + 1) of the n_top most variable genes; color = exposure, marker = batch."""
    lib = counts.sum(axis=1)
    logcpm = np.log2(counts.div(lib, axis=0) * 1e6 + 1)
    top = logcpm.var(axis=0).nlargest(n_top).index
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(logcpm[top])  # PCA centers each gene
    ev = pca.explained_variance_ratio_ * 100
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for (exposure, batch), idx in meta.groupby(["exposure", "batch"]).groups.items():
        rows = meta.index.get_indexer(idx)
        ax.scatter(pcs[rows, 0], pcs[rows, 1], c=EXPOSURE_COLORS.get(exposure, "gray"),
                   marker=BATCH_MARKERS.get(batch, "^"), s=60, edgecolor="black", linewidth=0.5,
                   label=f"{exposure}, {batch}")
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}% variance)")
    ax.set_title(f"PCA, log2(CPM+1), top {n_top} variable genes")
    ax.legend(fontsize=8, frameon=False)
    _save(fig, "01_pca.png")


def plot_overlap(summary: dict) -> None:
    labels = ["Shared", "PyDESeq2 only", "edgeR only"]
    values = [summary["n_shared"], summary["n_pydeseq2_only"], summary["n_edger_only"]]
    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(labels, values, color=["#55A868", "#DD8452", "#4C72B0"])
    ax.bar_label(bars, labels=[f"{v:,}" for v in values])
    ax.set_ylabel(f"Genes at FDR < {FDR}")
    ax.set_title(f"Significant genes: PyDESeq2 ({summary['n_sig_pydeseq2']:,}) vs "
                 f"published edgeR ({summary['n_sig_edger']:,})", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    _save(fig, "02_overlap.png")


def plot_logfc_scatter(j: pd.DataFrame, summary: dict) -> None:
    both = j.dropna(subset=["log2FoldChange", "logFC"])
    sig_both = (both["padj"] < FDR) & (both["FDR"] < FDR)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(both.loc[~sig_both, "logFC"], both.loc[~sig_both, "log2FoldChange"],
               s=3, alpha=0.3, color="gray", rasterized=True, label="other genes")
    ax.scatter(both.loc[sig_both, "logFC"], both.loc[sig_both, "log2FoldChange"],
               s=5, alpha=0.7, color="#55A868", rasterized=True, label="significant in both")
    lim = float(np.nanmax(np.abs(both[["logFC", "log2FoldChange"]].to_numpy())))
    ax.plot([-lim, lim], [-lim, lim], color="black", linewidth=0.7, linestyle="--")
    ax.set_xlabel("edgeR logFC (published)")
    ax.set_ylabel("PyDESeq2 log2 fold-change")
    ax.set_title(f"SD vs NSD fold-changes, Spearman rho = {summary['spearman_rho']:.3f} "
                 f"(n = {summary['n_spearman']:,})", fontsize=10)
    ax.legend(fontsize=8, frameon=False, markerscale=3)
    _save(fig, "03_logfc_scatter.png")


def plot_memory_genes(mem: pd.DataFrame) -> None:
    mem = mem.sort_values("log2FoldChange")
    y = np.arange(len(mem))
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.hlines(y, mem["ci_low"], mem["ci_high"], color="black", linewidth=1)
    sig = mem["significant"].to_numpy(dtype=bool)
    ax.scatter(mem["log2FoldChange"][sig], y[sig], s=50, color="black", zorder=3,
               label=f"padj < {FDR}")
    ax.scatter(mem["log2FoldChange"][~sig], y[~sig], s=50, facecolor="white",
               edgecolor="black", zorder=3, label=f"padj >= {FDR}")
    ax.axvline(0, color="gray", linewidth=0.7, linestyle="--")
    ax.set_yticks(y, mem.index, fontstyle="italic")
    ax.set_xlabel("PyDESeq2 log2 fold-change, SD vs NSD (95% CI)")
    ax.set_title("Selected memory-consolidation genes (descriptive)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    _save(fig, "04_memory_genes.png")


def main() -> None:
    if not SUMMARY_JSON.is_file() or not MEMORY_GENES_CSV.is_file():
        raise FileNotFoundError("results/summary.json or memory_genes.csv missing; run `make compare` first")
    summary = json.loads(SUMMARY_JSON.read_text())
    edger = load_edger()
    counts = load_counts(universe=edger.index)
    plot_pca(counts, load_meta(counts.index))
    plot_overlap(summary)
    plot_logfc_scatter(join_results(load_de_results(), edger), summary)
    plot_memory_genes(pd.read_csv(MEMORY_GENES_CSV, index_col=0))


if __name__ == "__main__":
    main()
