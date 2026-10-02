# Decisions (one line each)

- Library sizes are written to results/library_sizes.csv instead of a figure, because acceptance criterion 3 requires exactly four PNGs in figures/.
- Figures are PNG only (Section 7 mentions SVG, Section 8 requires PNG).
- PCA uses the authors' 22,582-gene tested universe, top 500 by variance of log2(CPM+1), sklearn PCA (gene-wise centering, no scaling).
- Notebooks are deferred: they are not in the acceptance criteria, and the 1-day build prioritizes a tested pipeline; the scripts in src/ are the single source of logic.
- README.md is generated from README.template.md with numbers from results/summary.json (acceptance criterion 4).
- Test thresholds follow Section 8: significant genes within 1,520 +/- 5%, Spearman within 0.974 +/- 0.02.
- NaN padj is treated as not significant and rows are never dropped before comparison.
- Large per-gene results CSV is gitignored (derived from unlicensed data); small summaries and figures are committed.
- environment.yml delegates to requirements.txt so versions are pinned in one place.
- CI does not diff README.md or figures, because last-digit drift across machines could cause false failures; tests enforce the Section 8 tolerances instead.
- Requires Python 3.12+, not the spec's 3.11+: the pip-freeze pins (e.g. anndata 0.13.4) have no Python 3.11 builds; full pipeline verified identical on Python 3.12.3 and 3.14.4.
