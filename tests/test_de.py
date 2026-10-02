"""Section 8 tolerances on the real PyDESeq2 fit (pipeline output if present, otherwise a fresh fit)."""
from __future__ import annotations

import pandas as pd
import pytest

from src.common import DE_RESULTS_FILE, RAW_DIR
from src.compare import comparison_summary, join_results, memory_gene_table
from src.download import EXPECTED
from src.load import load_counts, load_de_results, load_edger, load_meta

if not all((RAW_DIR / name).is_file() for name in EXPECTED):
    pytest.skip("raw data not downloaded; run `make download`", allow_module_level=True)


@pytest.fixture(scope="module")
def joined() -> pd.DataFrame:
    edger = load_edger()
    if DE_RESULTS_FILE.is_file():
        return join_results(load_de_results(), edger)
    from src.de import run_de

    counts = load_counts(universe=edger.index)
    return join_results(run_de(counts, load_meta(counts.index)), edger)


@pytest.fixture(scope="module")
def summary(joined: pd.DataFrame) -> dict:
    return comparison_summary(joined)


def test_significant_count_within_5pct(summary: dict) -> None:
    assert 1520 * 0.95 <= summary["n_sig_pydeseq2"] <= 1520 * 1.05, summary["n_sig_pydeseq2"]


def test_spearman_within_tolerance(summary: dict) -> None:
    assert abs(summary["spearman_rho"] - 0.974) <= 0.02, summary["spearman_rho"]


def test_shared_genes_same_direction(summary: dict) -> None:
    assert summary["direction_agreement_pct"] >= 99.0


def test_nan_padj_kept_not_dropped(joined: pd.DataFrame) -> None:
    assert len(joined) == 22582
    assert joined["padj"].isna().any()


def test_sign_convention(joined: pd.DataFrame) -> None:
    mem = memory_gene_table(joined)
    assert mem.at["Cirbp", "log2FoldChange"] < 0
    assert mem.at["Arc", "log2FoldChange"] > 0
