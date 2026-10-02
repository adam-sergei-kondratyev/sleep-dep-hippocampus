"""Load counts, sample metadata, and the authors' edgeR table, aligned (design doc 5.2)."""
from __future__ import annotations

import pandas as pd

from src.common import COUNTS_FILE, EDGER_FILE, LIBRARY_SIZES_CSV, RESULTS_DIR, SAMPLE_INFO_FILE


def load_edger() -> pd.DataFrame:
    """Authors' results. Index = versioned Ensembl ID; never join on the unversioned ensembl_gene_id column."""
    return pd.read_csv(EDGER_FILE, index_col=0)


def load_counts(universe: pd.Index | None = None) -> pd.DataFrame:
    """Samples x genes integer matrix, optionally restricted to (and ordered by) a gene universe."""
    counts = pd.read_csv(COUNTS_FILE, index_col=0)
    if universe is not None:
        missing = universe.difference(counts.index)
        if len(missing) > 0:
            raise KeyError(f"{len(missing)} universe genes absent from counts, e.g. {list(missing[:3])}")
        counts = counts.loc[universe]
    return counts.T.astype(int)


def load_meta(samples: pd.Index) -> pd.DataFrame:
    """batch and exposure per sample, rows in the same order as `samples`."""
    meta = pd.read_csv(SAMPLE_INFO_FILE, index_col=0)[["batch", "exposure"]]
    return meta.loc[samples]


def write_library_sizes(counts: pd.DataFrame) -> pd.Series:
    """QC table of total counts per sample (a CSV, so figures/ holds exactly the four required PNGs)."""
    sizes = counts.sum(axis=1).rename("library_size")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    sizes.to_csv(LIBRARY_SIZES_CSV)
    return sizes
