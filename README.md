# FIRE-mouse-2021-paper

This repository contains all of the code used to process and analyze the SPLiT-seq
data for Shabestari et al 2021.


## pre-processing

We used the `split-pipe` pipeline to quantify gene expression in single-nuclei,
check the `split-pipe.sh` script for the exact code.


## primary processing

Next we performed quality control filtering, doublet detection, clustering, and
data visualizations using Seurat and Scanpy in R and Python respectively. Check out the
`Processing.Rmd` script for the code. Furthermore, for the analysis without the
transplant conditions, check a very similar script `Processing-4conditions.Rmd`.

## differential expression

We compared different conditions and identified marker genes using Seurat's differential
gene expression platform, check out the following scripts: `DEGs.Rmd` and `parallel_DEGs.R`.

## cell communication analysis

We performed cell-cell communications analysis using CellChat: `cellchat.Rmd`.

## Co-expression network analysis

We used WGCNA to perform gene co-expression network analysis: `scWGCNA.Rmd`

## Additional plotting:

Additional plotting code can be found here: `plotting-for-paper.Rmd`
