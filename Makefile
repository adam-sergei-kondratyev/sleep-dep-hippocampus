.RECIPEPREFIX = >
.DELETE_ON_ERROR:
PY ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

COMMON  := src/__init__.py src/common.py src/load.py
RAW     := data/raw/manifest.json
DE      := results/pydeseq2_results.csv
SUMMARY := results/summary.json
FIGS    := figures/01_pca.png figures/02_overlap.png figures/03_logfc_scatter.png figures/04_memory_genes.png

.PHONY: all download de compare figures readme test clean distclean

all: $(FIGS) README.md test

download: $(RAW)
de: $(DE)
compare: $(SUMMARY)
figures: $(FIGS)
readme: README.md

$(RAW): src/download.py src/common.py
> $(PY) -m src.download

$(DE): $(RAW) src/de.py $(COMMON)
> $(PY) -m src.de

$(SUMMARY) results/memory_genes.csv results/comparison_summary.csv &: $(DE) src/compare.py $(COMMON)
> $(PY) -m src.compare

$(FIGS) results/library_sizes.csv &: $(SUMMARY) src/plots.py src/compare.py $(COMMON)
> $(PY) -m src.plots

README.md: README.template.md $(SUMMARY) src/render_readme.py
> $(PY) -m src.render_readme

test: $(SUMMARY)
> $(PY) -m pytest -q

clean:
> rm -rf results figures

distclean: clean
> rm -rf data/raw
