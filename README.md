<!-- README.md is generated from this template by `make readme`, so edit this file instead of README.md. -->
# Sleep deprivation and the mouse hippocampal transcriptome: an independent PyDESeq2 reanalysis

[![ci](https://github.com/adam-sergei-kondratyev/sleep-dep-hippocampus/actions/workflows/ci.yml/badge.svg)](https://github.com/adam-sergei-kondratyev/sleep-dep-hippocampus/actions/workflows/ci.yml)

**Question:** if I take the published data from Gaine et al. (2021), who kept mice awake for 5 hours and then measured gene expression in the hippocampus, and analyze it again with a different statistical method (PyDESeq2 in Python instead of the authors’ edgeR in R), will I find the same genes changing?

## Results

Using the same 22,582 genes the authors tested and the same model, which accounts for the two sequencing batches before comparing the 9 sleep-deprived (SD) mice with the 9 undisturbed (NSD) mice, PyDESeq2 found 1,521 genes at FDR < 0.1, compared with the 1,146 genes the authors reported. Of those 1,146 published genes, PyDESeq2 recovered 886 (77.3%), with 100.0% of them changing in the same direction in both analyses, and it also flagged 635 genes the authors did not. Across all 22,582 genes, the fold changes from the two methods were almost identical in rank (Spearman rho = 0.974), so both methods see the same effect of sleep loss but disagree on which genes clear the significance cutoff.

![Overlap of significant genes](figures/02_overlap.png)
![Fold-change agreement](figures/03_logfc_scatter.png)

## Why the two analyses differ

The two tools test each gene differently: PyDESeq2 uses the DESeq2 Wald test, while the authors used edgeR’s quasi-likelihood F-test (`glmQLFit(..., robust=TRUE)`), which adjusts for uncertainty in how much each gene varies between mice and is generally more conservative. They normalize for sequencing depth differently as well, with PyDESeq2 using median-of-ratios size factors and the authors using EDASeq GC-content normalization. PyDESeq2 also skips genes with too few reads to test reliably, leaving the adjusted p-value missing for 7,005 genes, which I counted as not significant. That filtering explains 73 (28.1%) of the 260 published genes PyDESeq2 did not recover, but the remaining 187 were tested and still missed the FDR < 0.1 cutoff because of the different test and normalization.

## Memory-consolidation genes

![Memory-consolidation genes](figures/04_memory_genes.png)

Since my senior capstone was a literature review on sleep and memory consolidation, I also looked at seven genes from that review (Arc, Nr4a1, Bdnf, Creb1, Camk2a, Prkaca, and Homer1), adding the paper’s top hits Cirbp and Upf2 for comparison, and plotted each gene’s log2 fold change with its 95% confidence interval. Because I picked these genes after the data were already published, I treat the figure as a description of what happened to them and not as a new hypothesis test. Arc, for example, went up after sleep deprivation (log2 fold change 0.418, padj 0.00047), while Cirbp went down (-0.52, padj 7.1e-20).

**Open question (Homer1).** The same lab’s earlier microarray study (PMID 22930738, GEO GSE33302) found its most significant change in probe 1436387_at, which decreased after sleep deprivation, but that probe is annotated to two genes at once (`C330006P03Rik///Homer1`), and in this RNA-seq data Homer1 does not change significantly at the gene level (log2 fold change 0.116, padj 0.28), leaving the reason for the mismatch as an open question.

## Limitations

I used the processed counts the authors shared without realigning the raw reads myself, and with the data coming from a single experiment on male mice using the whole hippocampus, the results cannot say anything about sex differences or specific cell types. Agreement between PyDESeq2 and edgeR shows that the analysis reproduces computationally, but confirming the biology would take new experiments. Exact counts can shift slightly between package versions, which is why `requirements.txt` pins every package; this README was generated with Python 3.12.3, PyDESeq2 0.5.4, pandas 3.0.6, numpy 2.5.3, and scipy 1.18.1.

## Reproduce

You need Python 3.12 or newer (the pinned requirements do not install on older versions), along with GNU Make 4.3 or newer and an internet connection (on Windows, use WSL). After cloning this repository, run:

~~
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
make all
~~

`make all` downloads the data, runs the analysis, writes the tables to `results/` and the plots to `figures/`, regenerates this README, and runs the tests. The data are downloaded from the authors’ repository each time and never stored here, as that repository has no license.

## Citations

- Gaine ME, Bahl E, Chatterjee S, Michaelson JJ, Abel T, Lyons LC. Altered hippocampal transcriptome dynamics following sleep deprivation. Molecular Brain. 2021;14:125. doi:10.1186/s13041-021-00835-1
- Processed data and original R code: https://github.com/ethanbahl/gaine2021_sleepdeprivation
- Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023. doi:10.1093/bioinformatics/btad547
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550. doi:10.1186/s13059-014-0550-8

## AI assistance

I used Claude (Anthropic) to help plan the project, write and debug the code, and draft this README, but I chose the study and the question, picked the memory-consolidation genes from my own literature review, ran every step on my own computer, and reviewed the results, so I am responsible for everything in this repository.

## License

The code is released under the MIT license (see `LICENSE`), but the data belong to the original authors.
