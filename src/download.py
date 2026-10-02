"""Download and verify the five processed data files (design doc 5.1).

Usage: python -m src.download
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from src.common import FDR, MANIFEST_FILE, RAW_DIR

BASE_URL = "https://raw.githubusercontent.com/ethanbahl/gaine2021_sleepdeprivation/main/data/"
REPO_URL = "https://github.com/ethanbahl/gaine2021_sleepdeprivation"
SAMPLES: tuple[str, ...] = tuple(f"NSD{i}" for i in range(1, 10)) + tuple(f"SD{i}" for i in range(1, 10))
PUBLISHED_SIG = 1146
RETRIES = 3


@dataclass(frozen=True)
class FileSpec:
    sep: str
    shape: tuple[int, int]
    columns: tuple[str, ...]


EXPECTED: dict[str, FileSpec] = {
    "counts_all.csv": FileSpec(",", (51826, 18), SAMPLES),
    "sample_info.csv": FileSpec(",", (18, 16), ("batch", "exposure")),
    "edger_diffexp_results.csv": FileSpec(
        ",", (22582, 14),
        ("ensembl_gene_id", "gene_name", "gene_type", "gc", "logFC", "logCPM", "F", "PValue", "FDR", "effect.size"),
    ),
    "gene_annotation.csv": FileSpec(",", (51826, 8), ("ensembl_gene_id", "gene_name", "gene_type", "gc")),
    "pmid_22930738_geo2r_results.txt": FileSpec(
        "\t", (45101, 7), ("adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title"),
    ),
}


class DataCheckError(RuntimeError):
    pass


def read_raw(name: str, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    path = raw_dir / name
    if not path.is_file():
        raise DataCheckError(f"{name}: not found at {path}")
    try:
        return pd.read_csv(path, sep=EXPECTED[name].sep, index_col=0)
    except Exception as exc:
        raise DataCheckError(f"{name}: could not be parsed ({exc!r})") from exc


def check_frame(name: str, df: pd.DataFrame) -> None:
    spec = EXPECTED[name]
    if df.shape != spec.shape:
        raise DataCheckError(f"{name}: expected shape {spec.shape}, got {df.shape}")
    missing = [c for c in spec.columns if c not in df.columns]
    if missing:
        raise DataCheckError(f"{name}: missing columns {missing}")
    if name == "counts_all.csv" and tuple(df.columns) != SAMPLES:
        raise DataCheckError(f"{name}: columns are not NSD1-NSD9, SD1-SD9 in that order")


def check_alignment(frames: dict[str, pd.DataFrame]) -> None:
    counts = frames["counts_all.csv"]
    meta = frames["sample_info.csv"]
    edger = frames["edger_diffexp_results.csv"]
    if set(meta.index) != set(counts.columns):
        raise DataCheckError("sample_info.csv: sample IDs do not match counts_all.csv columns")
    if not edger.index.is_unique:
        raise DataCheckError("edger_diffexp_results.csv: duplicated gene IDs in the index")
    if not edger.index.isin(counts.index).all():
        raise DataCheckError("edger_diffexp_results.csv: some gene IDs are absent from counts_all.csv")
    n_sig = int((edger["FDR"] < FDR).sum())
    if n_sig != PUBLISHED_SIG:
        raise DataCheckError(f"edger_diffexp_results.csv: expected {PUBLISHED_SIG} genes at FDR<{FDR}, got {n_sig}")


def verify(raw_dir: Path = RAW_DIR) -> dict[str, pd.DataFrame]:
    frames = {name: read_raw(name, raw_dir) for name in EXPECTED}
    for name, df in frames.items():
        check_frame(name, df)
    check_alignment(frames)
    return frames


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp, open(tmp, "wb") as fh:
                shutil.copyfileobj(resp, fh)
            tmp.replace(dest)
            return
        except urllib.error.HTTPError:
            raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == RETRIES:
                raise
            time.sleep(2 ** attempt)
        finally:
            tmp.unlink(missing_ok=True)


def _clone_fallback(names: list[str], raw_dir: Path) -> None:
    git = shutil.which("git")
    if git is None:
        raise DataCheckError("a raw URL returned 404 and the git clone fallback needs git")
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([git, "clone", "--depth", "1", REPO_URL, tmp], check=True)
        for name in names:
            src = Path(tmp) / "data" / name
            if not src.is_file():
                raise DataCheckError(f"{name}: not in the cloned repository's data/ folder")
            shutil.copy2(src, raw_dir / name)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(raw_dir: Path = RAW_DIR) -> None:
    manifest = {
        name: {"url": BASE_URL + name, "bytes": (raw_dir / name).stat().st_size, "sha256": _sha256(raw_dir / name)}
        for name in EXPECTED
    }
    MANIFEST_FILE.write_text(json.dumps(manifest, indent=2) + "\n")


def fetch_missing(raw_dir: Path = RAW_DIR) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    need_clone: list[str] = []
    for name in EXPECTED:
        if (raw_dir / name).is_file():
            print(f"present      {name}")
            continue
        try:
            _download(BASE_URL + name, raw_dir / name)
            print(f"downloaded   {name}")
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                raise
            print(f"404          {name} (git clone fallback)")
            need_clone.append(name)
    if need_clone:
        _clone_fallback(need_clone, raw_dir)


def main() -> int:
    try:
        fetch_missing()
        try:
            verify()
        except DataCheckError as first:
            print(f"cached data failed check ({first}); downloading again", file=sys.stderr)
            for name in EXPECTED:
                (RAW_DIR / name).unlink(missing_ok=True)
            fetch_missing()
            verify()
    except (DataCheckError, urllib.error.URLError) as exc:
        print(f"DATA CHECK FAILED: {exc}", file=sys.stderr)
        return 1
    write_manifest()
    print("all five files pass shape, column, and alignment checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
