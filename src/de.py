"""PyDESeq2 on the authors' tested gene universe (design doc 5.3).

Usage: python -m src.de
"""
from __future__ import annotations

import time

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from src.common import DE_RESULTS_FILE, FDR, RESULTS_DIR
from src.load import load_counts, load_edger, load_meta


def run_de(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Design ~ batch + exposure; positive log2FoldChange means higher in SD than NSD."""
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
    start = time.perf_counter()
    results = run_de(counts, load_meta(counts.index))
    elapsed = time.perf_counter() - start
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tmp = DE_RESULTS_FILE.with_name(DE_RESULTS_FILE.name + ".part")
    results.to_csv(tmp)
    tmp.replace(DE_RESULTS_FILE)
    print(f"PyDESeq2 done in {elapsed:.0f}s: {int((results['padj'] < FDR).sum())} genes at padj<{FDR}, "
          f"{int(results['padj'].isna().sum())} NaN padj")


if __name__ == "__main__":
    main()
