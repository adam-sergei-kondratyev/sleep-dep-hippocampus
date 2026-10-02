"""Aligned loaders for counts, sample metadata, and the published edgeR table."""
from __future__ import annotations

import pandas as pd

from src.common import COUNTS_FILE, DE_RESULTS_FILE, EDGER_FILE, SAMPLE_INFO_FILE


def load_edger() -> pd.DataFrame:
    """Index is the versioned Ensembl ID, the only safe join key."""
    return pd.read_csv(EDGER_FILE, index_col=0)


def load_counts(universe: pd.Index | None = None) -> pd.DataFrame:
    """Samples x genes int64 matrix, restricted to and ordered by `universe` when given."""
    counts = pd.read_csv(COUNTS_FILE, index_col=0)
    if universe is not None:
        missing = universe.difference(counts.index)
        if len(missing):
            raise KeyError(f"{len(missing)} universe genes absent from counts, e.g. {list(missing[:3])}")
        counts = counts.loc[universe]
    return counts.T.astype("int64")


def load_meta(samples: pd.Index) -> pd.DataFrame:
    return pd.read_csv(SAMPLE_INFO_FILE, index_col=0).loc[samples, ["batch", "exposure"]]


def load_de_results() -> pd.DataFrame:
    if not DE_RESULTS_FILE.is_file():
        raise FileNotFoundError(f"{DE_RESULTS_FILE} not found; run `make de` first")
    return pd.read_csv(DE_RESULTS_FILE, index_col=0, float_precision="round_trip")
