"""Download the authors' five processed data files and verify them (design doc 5.1).

Usage: python -m src.download
Files already present and passing every check are not downloaded again.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from src.common import FDR, RAW_DIR

BASE_URL = "https://raw.githubusercontent.com/ethanbahl/gaine2021_sleepdeprivation/main/data/"
REPO_URL = "https://github.com/ethanbahl/gaine2021_sleepdeprivation"
SAMPLES: tuple[str, ...] = tuple(f"NSD{i}" for i in range(1, 10)) + tuple(f"SD{i}" for i in range(1, 10))

# file name -> (separator, shape when read with index_col=0, columns that must be present)
EXPECTED: dict[str, tuple[str, tuple[int, int], tuple[str, ...]]] = {
    "counts_all.csv": (",", (51826, 18), SAMPLES),
    "sample_info.csv": (",", (18, 16), ("batch", "exposure")),
    "edger_diffexp_results.csv": (
        ",",
        (22582, 14),
        ("ensembl_gene_id", "gene_name", "gene_type", "gc", "logFC", "logCPM", "F", "PValue", "FDR", "effect.size"),
    ),
    "gene_annotation.csv": (",", (51826, 8), ("ensembl_gene_id", "gene_name", "gene_type", "gc")),
    "pmid_22930738_geo2r_results.txt": (
        "\t",
        (45101, 7),
        ("adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title"),
    ),
}


class DataCheckError(RuntimeError):
    """A data file does not match the verified shape, columns, or alignment."""


def read_raw(name: str, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    sep = EXPECTED[name][0]
    return pd.read_csv(raw_dir / name, sep=sep, index_col=0)


def check_file(name: str, raw_dir: Path = RAW_DIR) -> None:
    """Raise DataCheckError naming the file if it is missing, unparseable, or the wrong shape."""
    path = raw_dir / name
    if not path.is_file():
        raise DataCheckError(f"{name}: not found at {path}")
    _, shape, required = EXPECTED[name]
    try:
        df = read_raw(name, raw_dir)
    except Exception as exc:  # truncated or malformed download
        raise DataCheckError(f"{name}: could not be parsed ({exc!r})") from exc
    if df.shape != shape:
        raise DataCheckError(f"{name}: expected shape {shape}, got {df.shape}")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise DataCheckError(f"{name}: missing columns {missing}")
    if name == "counts_all.csv" and tuple(df.columns) != SAMPLES:
        raise DataCheckError(f"{name}: columns are not NSD1-NSD9, SD1-SD9 in that order")


def check_cross_file(raw_dir: Path = RAW_DIR) -> None:
    """Alignment facts the analysis depends on (design doc Section 2)."""
    counts = read_raw("counts_all.csv", raw_dir)
    meta = read_raw("sample_info.csv", raw_dir)
    edger = read_raw("edger_diffexp_results.csv", raw_dir)
    if set(meta.index) != set(counts.columns):
        raise DataCheckError("sample_info.csv: sample IDs do not match counts_all.csv columns")
    if not edger.index.is_unique:
        raise DataCheckError("edger_diffexp_results.csv: duplicated gene IDs in the index")
    if not edger.index.isin(counts.index).all():
        raise DataCheckError("edger_diffexp_results.csv: some gene IDs are absent from counts_all.csv")
    n_sig = int((edger["FDR"] < FDR).sum())
    if n_sig != 1146:
        raise DataCheckError(f"edger_diffexp_results.csv: expected 1146 genes at FDR<{FDR}, got {n_sig}")


def _download(url: str, dest: Path) -> None:
    """HTTPS GET to a .part file, then rename, so a failed download never leaves a half-written file."""
    tmp = dest.with_name(dest.name + ".part")
    with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    tmp.replace(dest)


def _clone_fallback(names: list[str], raw_dir: Path) -> None:
    """Design doc 5.1 fallback when a raw URL returns 404."""
    git = shutil.which("git")
    if git is None:
        raise DataCheckError("a raw URL returned 404 and the git clone fallback needs git, which is not installed")
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([git, "clone", "--depth", "1", REPO_URL, tmp], check=True)
        for name in names:
            src = Path(tmp) / "data" / name
            if not src.is_file():
                raise DataCheckError(f"{name}: not in the cloned repository's data/ folder")
            shutil.copy2(src, raw_dir / name)


def main() -> int:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    need_clone: list[str] = []
    for name in EXPECTED:
        try:
            check_file(name)
            print(f"ok (cached)  {name}")
            continue
        except DataCheckError:
            pass
        try:
            _download(BASE_URL + name, RAW_DIR / name)
            print(f"downloaded   {name}")
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                raise
            print(f"404          {name} (will use git clone fallback)")
            need_clone.append(name)
    try:
        if need_clone:
            _clone_fallback(need_clone, RAW_DIR)
        for name in EXPECTED:
            check_file(name)
        check_cross_file()
    except DataCheckError as exc:
        print(f"DATA CHECK FAILED: {exc}", file=sys.stderr)
        return 1
    print("all five files pass shape, column, and alignment checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
