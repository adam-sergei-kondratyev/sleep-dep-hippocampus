"""Compare PyDESeq2 with the published edgeR results and tabulate the memory genes (design doc 5.4, 5.5).

Usage: python -m src.compare
"""
from __future__ import annotations

import json
import platform
from importlib.metadata import version
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src.common import FDR, MEMORY_GENES, MEMORY_GENES_CSV, RESULTS_DIR, SUMMARY_CSV, SUMMARY_JSON
from src.load import load_de_results, load_edger

Z_95 = 1.959963984540054


def join_results(de: pd.DataFrame, edger: pd.DataFrame) -> pd.DataFrame:
    joined = de.join(edger[["gene_name", "logFC", "FDR"]], how="inner")
    if len(joined) != len(edger):
        raise ValueError(f"join kept {len(joined)} of {len(edger)} edgeR genes; gene universe mismatch")
    return joined


def _pct(part: int, whole: int) -> float:
    return round(100 * part / whole, 1) if whole else float("nan")


def comparison_summary(j: pd.DataFrame) -> dict[str, Any]:
    """NaN padj counts as not significant; rows are never dropped."""
    py_sig = j["padj"].lt(FDR).to_numpy()
    ed_sig = j["FDR"].lt(FDR).to_numpy()
    padj_nan = j["padj"].isna().to_numpy()
    lfc_py = j["log2FoldChange"].to_numpy()
    lfc_ed = j["logFC"].to_numpy()

    shared = py_sig & ed_sig
    edger_only = ed_sig & ~py_sig
    n_shared = int(shared.sum())
    n_ed_sig = int(ed_sig.sum())
    n_edger_only = int(edger_only.sum())
    n_edger_only_nan = int((edger_only & padj_nan).sum())

    finite = np.isfinite(lfc_py) & np.isfinite(lfc_ed)
    rho = float(spearmanr(lfc_py[finite], lfc_ed[finite]).statistic)
    same_dir = int((np.sign(lfc_py[shared]) == np.sign(lfc_ed[shared])).sum())

    return {
        "n_tested": int(len(j)),
        "n_sig_edger": n_ed_sig,
        "n_sig_edger_up": int((ed_sig & (lfc_ed > 0)).sum()),
        "n_sig_edger_down": int((ed_sig & (lfc_ed < 0)).sum()),
        "n_sig_pydeseq2": int(py_sig.sum()),
        "n_padj_nan": int(padj_nan.sum()),
        "n_shared": n_shared,
        "n_pydeseq2_only": int((py_sig & ~ed_sig).sum()),
        "n_edger_only": n_edger_only,
        "recovery_pct": _pct(n_shared, n_ed_sig),
        "n_edger_only_padj_nan": n_edger_only_nan,
        "n_edger_only_tested": n_edger_only - n_edger_only_nan,
        "edger_only_padj_nan_pct": _pct(n_edger_only_nan, n_edger_only),
        "direction_agreement_pct": _pct(same_dir, n_shared),
        "spearman_rho": round(rho, 3),
        "n_spearman": int(finite.sum()),
        "fdr": FDR,
        "python": platform.python_version(),
        "pydeseq2": version("pydeseq2"),
        "pandas": version("pandas"),
        "numpy": version("numpy"),
        "scipy": version("scipy"),
    }


def memory_gene_table(j: pd.DataFrame, genes: tuple[str, ...] = MEMORY_GENES) -> pd.DataFrame:
    sub = j[j["gene_name"].isin(genes)]
    n_rows = sub["gene_name"].value_counts()
    bad = {g: int(n_rows.get(g, 0)) for g in genes if n_rows.get(g, 0) != 1}
    if bad:
        raise ValueError(f"genes not matching exactly one row by gene_name: {bad}")
    t = sub.rename_axis("ensembl_id").reset_index().set_index("gene_name").loc[list(genes)]
    half = Z_95 * t["lfcSE"]
    return pd.DataFrame({
        "ensembl_id": t["ensembl_id"],
        "log2FoldChange": t["log2FoldChange"],
        "lfcSE": t["lfcSE"],
        "ci_low": t["log2FoldChange"] - half,
        "ci_high": t["log2FoldChange"] + half,
        "padj": t["padj"],
        "significant": t["padj"].lt(FDR),
        "edger_logFC": t["logFC"],
        "edger_FDR": t["FDR"],
    }).rename_axis("gene")


def main() -> None:
    j = join_results(load_de_results(), load_edger())
    summary = comparison_summary(j)
    mem = memory_gene_table(j)
    for gene in ("Homer1", "Arc", "Cirbp"):
        summary[f"{gene.lower()}_lfc"] = round(float(mem.at[gene, "log2FoldChange"]), 3)
        summary[f"{gene.lower()}_padj"] = float(f"{mem.at[gene, 'padj']:.2g}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    pd.Series(summary, name="value").rename_axis("metric").to_csv(SUMMARY_CSV)
    mem.to_csv(MEMORY_GENES_CSV)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
