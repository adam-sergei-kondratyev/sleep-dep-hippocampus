"""Compare PyDESeq2 results with the published edgeR results (design doc 5.4) and tabulate memory genes (5.5).

Usage: python -m src.compare   (requires results/pydeseq2_results.csv from src.de)
Writes results/summary.json (source of every README number), results/comparison_summary.csv,
and results/memory_genes.csv.
"""
from __future__ import annotations

import json
import platform
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src.common import (
    DE_RESULTS_FILE, FDR, MEMORY_GENES, MEMORY_GENES_CSV, RESULTS_DIR, SUMMARY_CSV, SUMMARY_JSON,
)
from src.load import load_edger


def load_de_results() -> pd.DataFrame:
    if not DE_RESULTS_FILE.is_file():
        raise FileNotFoundError(f"{DE_RESULTS_FILE} not found; run `make de` first")
    return pd.read_csv(DE_RESULTS_FILE, index_col=0)


def join_results(de: pd.DataFrame, edger: pd.DataFrame) -> pd.DataFrame:
    """Join on the versioned Ensembl index (never on gene symbol). Keeps NaN padj rows."""
    joined = de.join(edger[["gene_name", "logFC", "FDR"]], how="inner")
    if len(joined) != len(edger):
        raise ValueError(f"join kept {len(joined)} of {len(edger)} edgeR genes; gene universe mismatch")
    return joined


def comparison_summary(j: pd.DataFrame) -> dict[str, Any]:
    """Every number the README quotes. NaN padj counts as not significant."""
    py_sig = (j["padj"] < FDR).to_numpy()  # NaN < FDR is False
    ed_sig = (j["FDR"] < FDR).to_numpy()
    shared = py_sig & ed_sig
    edger_only = ed_sig & ~py_sig
    padj_nan = j["padj"].isna().to_numpy()

    both = j.dropna(subset=["log2FoldChange", "logFC"])
    rho = float(spearmanr(both["log2FoldChange"], both["logFC"]).statistic)

    sh = j.loc[shared]
    same_dir = int((np.sign(sh["log2FoldChange"]) == np.sign(sh["logFC"])).sum())
    n_ed_sig = int(ed_sig.sum())
    return {
        "n_tested": int(len(j)),
        "n_sig_edger": n_ed_sig,
        "n_sig_edger_up": int((ed_sig & (j["logFC"] > 0).to_numpy()).sum()),
        "n_sig_edger_down": int((ed_sig & (j["logFC"] < 0).to_numpy()).sum()),
        "n_sig_pydeseq2": int(py_sig.sum()),
        "n_padj_nan": int(padj_nan.sum()),
        "n_shared": int(shared.sum()),
        "n_pydeseq2_only": int((py_sig & ~ed_sig).sum()),
        "n_edger_only": int(edger_only.sum()),
        "recovery_pct": round(100 * int(shared.sum()) / n_ed_sig, 1),
        "n_edger_only_padj_nan": int((edger_only & padj_nan).sum()),
        "n_edger_only_tested": int((edger_only & ~padj_nan).sum()),
        "n_padj_nan_cooks": int(j["pvalue"].isna().sum()),  # Cook's outliers get pvalue NaN
        "direction_agreement_pct": round(100 * same_dir / max(len(sh), 1), 1),
        "spearman_rho": round(rho, 3),
        "n_spearman": int(len(both)),
        "fdr": FDR,
        "python": platform.python_version(),
        "pydeseq2": version("pydeseq2"),
        "pandas": version("pandas"),
        "numpy": version("numpy"),
        "scipy": version("scipy"),
    }


def memory_gene_table(j: pd.DataFrame) -> pd.DataFrame:
    """Section 5.5 genes; each symbol is unique in the edgeR table (checked here, not assumed)."""
    rows = []
    for gene in MEMORY_GENES:
        hit = j[j["gene_name"] == gene]
        if len(hit) != 1:
            raise ValueError(f"{gene}: expected exactly 1 row by gene_name, found {len(hit)}")
        r = hit.iloc[0]
        rows.append({
            "gene": gene,
            "ensembl_id": hit.index[0],
            "log2FoldChange": r["log2FoldChange"],
            "lfcSE": r["lfcSE"],
            "ci_low": r["log2FoldChange"] - 1.96 * r["lfcSE"],
            "ci_high": r["log2FoldChange"] + 1.96 * r["lfcSE"],
            "padj": r["padj"],
            "significant": bool(r["padj"] < FDR),
            "edger_logFC": r["logFC"],
            "edger_FDR": r["FDR"],
        })
    return pd.DataFrame(rows).set_index("gene")


def main() -> None:
    j = join_results(load_de_results(), load_edger())
    summary = comparison_summary(j)
    mem = memory_gene_table(j)
    for gene in ("Homer1", "Arc", "Cirbp"):
        summary[f"{gene.lower()}_lfc"] = round(float(mem.loc[gene, "log2FoldChange"]), 3)
        summary[f"{gene.lower()}_padj"] = float(f"{mem.loc[gene, 'padj']:.2g}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n")
    pd.Series(summary, name="value").rename_axis("metric").to_csv(SUMMARY_CSV)
    mem.to_csv(MEMORY_GENES_CSV)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
