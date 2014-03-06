GC-content normalization for RNA-Seq data
========================================================

**Paper:** Risso, D., Schwartz, K., Sherlock, G. & Dudoit, S. GC-content normalization for RNA-Seq data. BMC Bioinformatics 12, 480 (2011).

**Group:** Jun Hwang, Caitlin McHugh, David Whitney.


Introduction
=============
_RNA-Seq_ is a 'recent' high-throughput sequencing assay technology.

- Wang, Gerstein, Snyder (2009) claimed that RNA-Seq "_can capture transcriptome dynamics across different tissues or conditions without sophisticated normalization of data sets._"
- A bold claim!

RNA-Seq: Overview
============
1. mRNA $\rightarrow$ cDNA fragments.

2. cDNA sequenced into millions of short (25-100 bases) reads.

3. Reads mapped (aligned) to reference genome.

4. Read count for given gene ~ abundance of transcript in sample.

RNA-Seq: Expected Sources of Read Count Error
============
- Gene length
- Variation between replicate lanes from differing sequencing depth.
- GC-content, which tends to be lane specific
- Mappability (how common sequence pattern is within genome)

RNA-Seq: Sophisticated Normalization
==============
We can identify 2 broad types of effect on read counts:

1. Within-lane (gene length, GC-content, mappability)

2. Between-lane (sequencing depth)

Thus do normalization for both types of effect.

Methods for Within-lane Normalization
=====================================
The authors present three approaches to within-lane normalization. The goal of each method is to adjust for dependence of read counts on GC-content.

- **Regression normalization.** Regress the counts on GC-content and subtract the loess fit from the counts to remove dependence.
- **Global-scaling normalization** Counts are rescaled through values of a summary statistic within a GC-content stratum

- **Full-quantile normalization.** Force the distribution of counts to be the same across bins of genes based on GC-content.

Implementations of within-lane and between-lane methods: `EDASeq` package in Bioconductor.

Methods for Between-lane Normalization
======================================
The authors opt for **full-quantile** normalization method covered by Bullard, et al (2010).


Data: What Samples?
=============
14 Saccharomyces cerevisiae **yeast** samples were used to test the normalization techniques.
Among the 14 samples, there are
- Three different growth conditions: standard YP Glucose (YPD), Delft Glucose (Del) and YP Glycerol (Gly)
- Two library preparation protocols: 1 and 2
- Five flow cells: 428R1, 4328B, 61MKN, 247L, 62OAY

Factor | Biological or Technical?
 ---|---
Growth conditions | biological
Library prep | technical
Flow cell | technical

Data: How Did We Prepare Them?
=============
- Data were _downloaded_ from NCBI's Sequence Read Archive (SRA), accession number SRA048710
- Then data were _aligned_ to the yeast reference genome using the `RBowtie` package

<small>**Note:** the authors did not specify which reference genome they used for alignment</small>

- The aligned SAM files were converted to BAM format using the `Rsamtools` package
- Finally, genes were _filtered_ out if

max(average read count within growth conditions)<10 

**We were unable to mirror exactly the counts the authors present**


Data: Download Again
=============
After completing the SRA downloads and running the alignment, we discovered the authors provide the read counts for all 14 yeast samples online.

All remaining analyses use read counts from **5,690 filtered genes**, as provided by the authors.


Analyses Using New Methods
==========================
To account for GC-content bias within lanes and technical differences between lanes, there are three steps to each analysis.

1. Apply _within-lane_ normalization.

2. Apply _between-lane_ normalization.

3. Perform differential expression analysis (package: `edgeR`).

Regression Normalization
========================
For gene $j = 1, \ldots, J$

0. Log-read counts $y_j$ are regressed on GC-content $x_j$ via loess robust local regression, obtaining $\hat{y}_j$

1. Calculate the residuals $y_j - \hat{y}_j$

2. Calculate a summary statistic such (e.g. the median), $T(y_1,\ldots,y_J)$

3. The normalized values are given by $y_j' = y_j - \hat{y}_j + T(y_1,\ldots,y_J)$

**Note:** In `EDASeq` it is not possible to specify bandwidth for the loess regression or the statistic $T$, which is presumed to be the median (based on the paper).

Global-scaling Normalization
==========================
 **Intuition:** Bias due to GC-content can be corrected by comparison of expression levels within similar GC-content ranges

 **Specifics:** On log-scale: $y_j' = y_j - T(y_{j'} : j' \in k(j))) + T(y_1,...,y_j)$ for gene j and summary statistic T. 

0. Bin J genes into K bins based on GC-content

1. Determine T among genes y_j' within GC-bin k(j) within lane

2. Determine T among all genes within lane

3. Original scale: $\exp(y_j') = \frac{exp(y_j)}{T(\exp(\overrightarrow{y_{j'}}))/T(\exp(\overrightarrow{y}))}$ 

Full-quantile Normalization
==========================
 **Goal:** make _distribution_ of log(read counts) the same across bins

0. Bin J genes into K bins based on GC-content

1. Form a matrix of log(read counts) $X_{\frac{J}{K} \times K}$ where each column is a bin
 
2. Sort each column (bin) of $X_{\frac{J}{K} \times K}$ = $X_{sort}$

3. Take mean by row, create vector $\mathbf{m}$

4. Make $X'_{sort}$ that has entries of row means

5. Put columns back in original order of $X$, create $X_{norm}$

Full-quantile Normalization: An Example
==========================
$X$ =

  |  |  |
 ---|---|---
**3** | **7** | **2**
_5_ | _1_ | _6_
1 | 4 | 3

$X_{sort}$ = 

 | | |
 ---|---|---
 _5_ | **7** | _6_
**3** | 4 | 3
1 | _1_ | **2**
***
$\mathbf{m}$ = (6,3.3,1.3)

$X'_{sort}$ = 

 | | |
 ---|---|---
 _6_ | **6** | _6_
**3.3** | 3.3 | 3.3
1.3 | _1.3_ | **1.3**

$X_{norm}$ = 

 | | |
 ---|---|---
**3.3**  | **6** |**1.3**
_6_ | _1.3_ | _6_ 
1.3 |  3.3 | 3.3


Evaluation of New Methods
==========================
Bias, MSE, Type I error rate are considered for differential expression (DE) analysis.
 - We choose the 8 yeast cultures grown via YPD
 - Group them into two groups of 4 lanes each
 - Take the log-ratio of the normalized counts
 - Test whether the counts vary between the two groups 
 - **There should be no difference**

Results of Normalization Methods: Read Count
==========================
GC-content differs across cultures but is similar for the same cultures.

![Fig1](Images/figure_1.png)


Results of Normalization Methods: Log fold change
==========================
Fold change vs GC-content in the same culture with different flow cells (left) and
different cultures from the same flow cell (right).
![Fig2](Images/figure_2.png)

Results of Normalization Methods: Normalized log fold change
==========================
Normalized fold change vs GC-content using proposed normalization techniques.

Loess-normalized log fold change
==========================
![Fig3a](Images/figure_3a.png)

GS-normalized log fold change
==========================
![Fig3b](Images/figure_3b.png)

FQ-normalized log fold change
==========================
![Fig3c](Images/figure_3c.png)


Results of Normalization Methods: Bias
==========================
![FigS9](Images/figure_S9.png)

Results of Normalization Methods: MSE
==========================
![FigS10](Images/figure_S10.png)

Results of Normalization Methods: Type I error
==========================
![FigS12b](Images/figure_S12b.png)

Results of Normalization Methods: Type I error
==========================
![FigS12c](Images/figure_S12c.png)

Results of Normalization Methods: Type I error
==========================
![FigS12d](Images/figure_S12d.png)

Conclusions
==========================
Why should we use the methods developed by Risso et al.?
 - Easy to implement: `EDASeq` package for R available

Better than the state-of-the-art? Better than Hansen, Brenner, Dudoit (2011)?
 - Lower bias and MSE of expression fold-changed estimates
 - Lower Type I error and better p-values when testing differential expression

Bonus:
 - Design and results indicate library preparation as source of GC bias



Remarks on Reproducibility
==========================
Authors providing their counts is a mixed blessing:
 - Able to reproduce many analyses exactly
 - But what does this verify, given author-supplied data and authors' package?
 - Authors' documentation should be sufficient to replicate results exactly
 - Reference genome should be given explicitly
 - Do we know if there were issues beyond reference genome?

