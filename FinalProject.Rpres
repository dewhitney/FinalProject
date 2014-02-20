GC-content normalization for RNA-Seq data
========================================================

**Paper:** Risso, D., Schwartz, K., Sherlock, G. & Dudoit, S. GC-content normalization for RNA-Seq data. BMC Bioinformatics 12, 480 (2011).

**Group:** Jun Hwang, Caitlin McHugh, David Whitney.


```{r, echo=FALSE, cache=TRUE, background=TRUE}
source("http://bioconductor.org/biocLite.R")
biocLite(c("SRAdb","EDASeq","edgeR"))
library(SRAdb)
library(EDASeq)
library(edgeR)
```

Introduction
=============
_Discuss RNA-seq methodology briefly. Redefine lanes._

_Why do we need to do within/between-lane normalizations?_

The authors present three approaches to within-lane normalization. The goal of each method of normalization is to adjust for dependence of read counts on GC-content. The three approaches discuss are as follows.

1. **Regression normalization.** Regress the counts on GC-content and subtract the loess fit from the counts to remove dependence.

2. Global-scaling normalization

3. Full-quantile normalization

Implementations of the within-lane and between-lane methods used in the paper are available from the EDASeq package in R/Bioconductor.

Data
=============
Yeast (SRA048710) and MAQC (SRA010153) datasets were analyzed using the new approaches. _Some discussion of the data background._ We seek to reproduce the analyses of the authors as described in the paper.

Analyses Using New Methods
==========================
To account for GC-content bias within lanes and technical differences between lanes, there are three steps to each analysis.

1. Apply _within-lane_ normalization.

2. Apply _between-lane_ normalization.

3. Perform differential expression analysis (package: edgeR).

Evaluation of New Methods
==========================
Bias, MSE, Type I error rate, p-value distributions.