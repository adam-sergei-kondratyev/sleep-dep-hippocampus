"""Logic tests on small synthetic tables; no download needed."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.compare import comparison_summary, join_results, memory_gene_table
from src.render_readme import render


@pytest.fixture
def toy() -> pd.DataFrame:
    idx = pd.Index([f"G{i}.1" for i in range(6)])
    de = pd.DataFrame(
        {
            "log2FoldChange": [1.0, -1.0, 0.5, 0.2, -0.3, np.nan],
            "lfcSE": [0.1, 0.1, 0.2, 0.2, 0.2, np.nan],
            "pvalue": [1e-6, 1e-6, 0.01, 0.5, 0.2, np.nan],
            "padj": [1e-5, 1e-5, 0.05, np.nan, 0.3, np.nan],
        },
        index=idx,
    )
    edger = pd.DataFrame(
        {
            "gene_name": ["A", "B", "C", "D", "E", "F"],
            "logFC": [0.9, -1.1, 0.4, 0.3, -0.2, 0.0],
            "FDR": [1e-4, 1e-4, 0.5, 0.01, 0.05, 0.9],
        },
        index=idx,
    )
    return join_results(de, edger)


def test_nan_padj_is_not_significant_and_counted(toy: pd.DataFrame) -> None:
    s = comparison_summary(toy)
    assert s["n_sig_pydeseq2"] == 3
    assert s["n_sig_edger"] == 4
    assert s["n_shared"] == 2
    assert s["n_edger_only"] == 2
    assert s["n_edger_only_padj_nan"] == 1
    assert s["n_edger_only_tested"] == 1
    assert s["n_padj_nan"] == 2
    assert s["recovery_pct"] == 50.0
    assert s["direction_agreement_pct"] == 100.0


def test_spearman_skips_missing_fold_changes(toy: pd.DataFrame) -> None:
    s = comparison_summary(toy)
    assert s["n_spearman"] == 5
    assert s["spearman_rho"] == 1.0


def test_join_rejects_universe_mismatch(toy: pd.DataFrame) -> None:
    de = toy[["log2FoldChange", "lfcSE", "pvalue", "padj"]].iloc[:-1]
    edger = toy[["gene_name", "logFC", "FDR"]]
    with pytest.raises(ValueError, match="universe mismatch"):
        join_results(de, edger)


def test_memory_table_ci_and_order(toy: pd.DataFrame) -> None:
    mem = memory_gene_table(toy, genes=("C", "A"))
    assert list(mem.index) == ["C", "A"]
    assert mem.at["A", "ci_high"] - mem.at["A", "ci_low"] == pytest.approx(2 * 1.959963984540054 * 0.1)
    assert bool(mem.at["C", "significant"]) is True


def test_memory_table_rejects_missing_gene(toy: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match="Zzz"):
        memory_gene_table(toy, genes=("A", "Zzz"))


def test_render_formats_and_rejects_unknown_keys() -> None:
    assert render("{{n}} genes, rho {{ r }}", {"n": 22582, "r": 0.974}) == "22,582 genes, rho 0.974"
    with pytest.raises(KeyError, match="missing"):
        render("{{missing}}", {})
