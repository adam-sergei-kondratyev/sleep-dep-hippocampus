"""Shared paths and constants."""
from __future__ import annotations

from pathlib import Path

ROOT: Path = Path(__file__).resolve().parent.parent
RAW_DIR: Path = ROOT / "data" / "raw"
RESULTS_DIR: Path = ROOT / "results"
FIGURES_DIR: Path = ROOT / "figures"

COUNTS_FILE: Path = RAW_DIR / "counts_all.csv"
SAMPLE_INFO_FILE: Path = RAW_DIR / "sample_info.csv"
EDGER_FILE: Path = RAW_DIR / "edger_diffexp_results.csv"
MANIFEST_FILE: Path = RAW_DIR / "manifest.json"

DE_RESULTS_FILE: Path = RESULTS_DIR / "pydeseq2_results.csv"
SUMMARY_JSON: Path = RESULTS_DIR / "summary.json"
SUMMARY_CSV: Path = RESULTS_DIR / "comparison_summary.csv"
MEMORY_GENES_CSV: Path = RESULTS_DIR / "memory_genes.csv"
LIBRARY_SIZES_CSV: Path = RESULTS_DIR / "library_sizes.csv"

FDR: float = 0.1
MEMORY_GENES: tuple[str, ...] = (
    "Arc", "Nr4a1", "Bdnf", "Creb1", "Camk2a", "Prkaca", "Homer1", "Cirbp", "Upf2",
)
