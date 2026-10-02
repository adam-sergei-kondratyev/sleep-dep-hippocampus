"""Differential expression with PyDESeq2 on the authors' tested gene universe (design doc 5.3).

Usage: python -m src.de
"""
from __future__ import annotations

import time

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from src.common import DE_RESULTS_FILE, FDR, RESULTS_DIR
from src.load import load_counts, load_edger, load_meta, write_library_sizes


def run_de(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Design ~ batch + exposure; contrast SD vs NSD, so positive log2FoldChange = higher after sleep deprivation."""
    if not counts.index.equals(meta.index):
        raise ValueError("counts and metadata sample order differ")
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~ batch + exposure", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["exposure", "SD", "NSD"], quiet=True)
    stats.summary()
    return stats.results_df


def main() -> None:
    edger = load_edger()
    counts = load_counts(universe=edger.index)
    meta = load_meta(counts.index)
    write_library_sizes(counts)
    start = time.time()
    results = run_de(counts, meta)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results.to_csv(DE_RESULTS_FILE)
    n_sig = int((results["padj"] < FDR).sum())
    n_nan = int(results["padj"].isna().sum())
    print(f"PyDESeq2 done in {time.time() - start:.0f}s: {n_sig} genes at padj<{FDR}, {n_nan} NaN padj")


if __name__ == "__main__":
    main()
