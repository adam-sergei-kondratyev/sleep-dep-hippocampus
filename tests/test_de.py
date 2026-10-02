"""Runs the real PyDESeq2 model (about 1 minute) and checks design doc Section 8 tolerances."""
from __future__ import annotations

import pandas as pd
import pytest

from src.common import FDR, RAW_DIR
from src.compare import comparison_summary, join_results, memory_gene_table
from src.download import EXPECTED
from src.load import load_counts, load_edger, load_meta

if not all((RAW_DIR / name).is_file() for name in EXPECTED):
    pytest.skip("raw data not downloaded; run `make download`", allow_module_level=True)


@pytest.fixture(scope="module")
def joined() -> pd.DataFrame:
    from src.de import run_de

    edger = load_edger()
    counts = load_counts(universe=edger.index)
    return join_results(run_de(counts, load_meta(counts.index)), edger)


def test_significant_count_within_5pct(joined: pd.DataFrame) -> None:
    n_sig = int((joined["padj"] < FDR).sum())
    assert 1520 * 0.95 <= n_sig <= 1520 * 1.05, n_sig


def test_spearman_within_tolerance(joined: pd.DataFrame) -> None:
    rho = comparison_summary(joined)["spearman_rho"]
    assert abs(rho - 0.974) <= 0.02, rho


def test_shared_genes_same_direction(joined: pd.DataFrame) -> None:
    assert comparison_summary(joined)["direction_agreement_pct"] >= 99.0


def test_nan_padj_kept_not_dropped(joined: pd.DataFrame) -> None:
    assert len(joined) == 22582
    assert joined["padj"].isna().sum() > 0


def test_sign_convention(joined: pd.DataFrame) -> None:
    mem = memory_gene_table(joined)
    assert mem.loc["Cirbp", "log2FoldChange"] < 0  # published: lower after sleep deprivation
    assert mem.loc["Arc", "log2FoldChange"] > 0
