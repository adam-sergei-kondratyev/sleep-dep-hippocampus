"""Data loading and alignment (design doc Section 7: shapes, sample order, group counts)."""
from __future__ import annotations

import pandas as pd
import pytest

from src.common import RAW_DIR
from src.download import EXPECTED, check_cross_file, check_file
from src.load import load_counts, load_edger, load_meta

if not all((RAW_DIR / name).is_file() for name in EXPECTED):
    pytest.skip("raw data not downloaded; run `make download`", allow_module_level=True)


@pytest.mark.parametrize("name", list(EXPECTED))
def test_file_shapes(name: str) -> None:
    check_file(name)


def test_cross_file_alignment() -> None:
    check_cross_file()


def test_counts_universe_and_order() -> None:
    edger = load_edger()
    counts = load_counts(universe=edger.index)
    assert counts.shape == (18, 22582)
    assert list(counts.columns) == list(edger.index)
    assert (counts.dtypes == "int64").all()
    assert (counts.to_numpy() >= 0).all()


def test_group_and_batch_counts() -> None:
    counts = load_counts(universe=load_edger().index)
    meta = load_meta(counts.index)
    assert meta.index.equals(counts.index)
    assert meta["exposure"].value_counts().to_dict() == {"NSD": 9, "SD": 9}
    assert meta["batch"].value_counts().to_dict() == {"batch2": 10, "batch1": 8}
    assert pd.crosstab(meta["batch"], meta["exposure"]).loc["batch1"].tolist() == [4, 4]


def test_edger_join_key_is_versioned_id() -> None:
    edger = load_edger()
    assert edger.index.is_unique
    assert edger.index.str.contains(r"\.\d+$").all()  # versioned, unlike the ensembl_gene_id column
