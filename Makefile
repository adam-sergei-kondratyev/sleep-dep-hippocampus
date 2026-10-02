# GNU Make 3.82+; ">" replaces the tab recipe prefix so copy-paste cannot break the file.
.RECIPEPREFIX = >
PY ?= python
.PHONY: all download de compare figures readme test clean

all: download de compare figures readme test

download:
> $(PY) -m src.download

de:
> $(PY) -m src.de

compare:
> $(PY) -m src.compare

figures:
> $(PY) -m src.plots

readme:
> $(PY) -m src.render_readme

test:
> $(PY) -m pytest -q

clean:
> rm -rf results figures
