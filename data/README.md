# Data

`make download` (or `python -m src.download`) fetches five processed files from https://github.com/ethanbahl/gaine2021_sleepdeprivation (`data/` folder, branch `main`) into `data/raw/`, then checks each file's shape and columns and fails loudly on any mismatch. Sizes and SHA-256 hashes of the verified files are written to `data/raw/manifest.json`. That repository has no license, so the data are not committed here (`data/raw/` is gitignored).

Source study: Gaine et al., Molecular Brain 2021;14:125, doi:10.1186/s13041-021-00835-1. Raw reads: GEO GSE166831.
