# S2N R package

# Introduction
Existing methods have been proposed to identify sample-specific networks (a network for a sample) from bulk, single-cell and spatial transcriptomics data with different programming languages (e.g. Matlab, Python, R). To facilitate comparison and the selection of a more suitable method for exploring sample-specific networks in transcriptomics data, we implement S2N R package.

A schematic illustration of **S2N** is shown in the folowing.

<p align="center">
  <img src="https://github.com/zhangjunpeng411/S2N/blob/main/S2N_schematic_illustration.png" alt="Schematic illustration of S2N" border="0.1">
</p>

S2N implements 10 in silico methods to infer sample-specific networks from transcriptomics data (including bulk, single-cell and spatial gene expression data). Based on the identified sample-specific networks, S2N performs downstream analyses, including network validation, identification of hubs and modules, and sample similarity calculation.

# Install S2N R package

```{r, eval=TRUE, include=TRUE}
install.packages("remotes")
library(remotes)
remotes::install_github("zhangjunpeng411/S2N")
```
# Load BRCA sample data
The BRCA sample data contains the following:
1. Bulk RNA-seq data, single-cell RNA-seq (scRNA-seq) data, and spatial transcriptomic (ST) data. 
2. Mutation data.
3. Groundtruth of miRNA-target interactions.
4. Groundtruth of lncRNA-target interactions.
5. Groundtruth of TF-target interactions.

```{r, eval=TRUE, include=TRUE}
data(bulk_ncR_tar_normal_cancer)
data(scRNA_ncR_tar_normal_cancer)
data(scRNA_ncR_tar_count_cancer)
data(ST_ncR_tar_cancer)
data(ST_ncR_tar_count_cancer)
data(mut_ncR)
data(miRTarget_Groundtruth)
data(lncRTarget_Groundtruth)
data(TFTarget_Groundtruth)
```

# A quick example for Inferring sample-specific networks
S2N provides 10 computational methods to identify sample-specific network, including SSN (Single-Sample Network), Paired-SSN (Paired-Single-Sample Network), LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples), SCS (single-sample controller strategy), CSN (cell-specific network), Scan.perturb (Sample-speCific miRNA regulAtioN with statistical perturbation strategy), Scan.interp (Sample-speCific miRNA regulAtioN with linear interpolation strategy), SWEET (Sample-specific weighted correlation network), SINUM (SIngle-cell Network Using Mutual information), and BONOBO (Bayesian Optimized Networks Obtained By assimilating Omics data).

```{r echo=FALSE, results='hide', message=FALSE}
library(S2N)

# bulk data
bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')

# scRNA data
scRNA_LIONESS_res <- LIONESS(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = 'pearson')

# ST data
ST_LIONESS_res <- LIONESS(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], cormethod = 'pearson')
```

# License
GPL-3

