"""Data loading and alignment: shapes, sample order, group and batch counts."""
from __future__ import annotations

import pandas as pd
import pytest

from src.common import RAW_DIR
from src.download import EXPECTED, DataCheckError, check_alignment, check_frame, read_raw
from src.load import load_counts, load_edger, load_meta

if not all((RAW_DIR / name).is_file() for name in EXPECTED):
    pytest.skip("raw data not downloaded; run `make download`", allow_module_level=True)


@pytest.fixture(scope="module")
def frames() -> dict[str, pd.DataFrame]:
    return {name: read_raw(name) for name in EXPECTED}


@pytest.mark.parametrize("name", list(EXPECTED))
def test_file_shapes(frames: dict[str, pd.DataFrame], name: str) -> None:
    check_frame(name, frames[name])


def test_cross_file_alignment(frames: dict[str, pd.DataFrame]) -> None:
    check_alignment(frames)


def test_counts_universe_and_order() -> None:
    edger = load_edger()
    counts = load_counts(universe=edger.index)
    assert counts.shape == (18, 22582)
    assert counts.columns.equals(edger.index)
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
    assert edger.index.str.contains(r"\.\d+$").all()


def test_truncated_file_is_rejected(frames: dict[str, pd.DataFrame]) -> None:
    with pytest.raises(DataCheckError, match="sample_info.csv: expected shape"):
        check_frame("sample_info.csv", frames["sample_info.csv"].head(5))
