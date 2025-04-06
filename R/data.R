#' bulk ncRNA expression data in breast cancer
#'
#' @docType data
#' @name bulk_ncR_cancer
#' @aliases bulk_ncR_cancer
#' @format bulk_ncR_cancer: A data frame object with six BRCA samples (columns) and 4794 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in breast cancer is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' bulk ncRNA expression data in normal breast
#'
#' @docType data
#' @name bulk_ncR_normal
#' @aliases bulk_ncR_normal
#' @format bulk_ncR_normal: A data frame object with six normal samples (columns) and 4929 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in normal breast is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' bulk target gene expression data in breast cancer
#'
#' @docType data
#' @name bulk_tar_cancer
#' @aliases bulk_tar_cancer
#' @format bulk_tar_cancer: A data frame object with six BRCA samples (columns) and 22427 target genes (rows).
#' @details The matched lncRNA and mRNA expression data in breast cancer is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' bulk target gene expression data in normal breast
#'
#' @docType data
#' @name bulk_tar_normal
#' @aliases bulk_tar_normal
#' @format bulk_tar_normal: A dataframe object with six normal samples (columns) and 22427 target genes (rows).
#' @details The matched lncRNA and mRNA expression data in normal breast is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' single-cell ncRNA expression data in breast cancer
#'
#' @docType data
#' @name scRNA_ncR_cancer
#' @aliases scRNA_ncR_cancer
#' @format scRNA_ncR_cancer: A data frame object with 736 barcodes of CID3946 breast cancer cells (columns) and 147 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in breast cancer cells is obtained from Wu et al.
#' The data is pre-processed by using Seurat R package.
#' The matched miRNA and lncRNA expression data is regarded as single-cell ncRNA expression data in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' single-cell ncRNA expression data in normal breast
#'
#' @docType data
#' @name scRNA_ncR_normal
#' @aliases scRNA_ncR_normal
#' @format scRNA_ncR_normal: A data frame object with 2430 barcodes of normal breast cells (columns) and 343 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in normal breast cells is obtained from GEO (GSE153889).
#' The data is pre-processed by using Seurat R package.
#' The matched miRNA and lncRNA expression data is regarded as single-cell ncRNA expression data in normal breast.
#' @references Martin Carli JF, Trahan GD, Jones KL, Hirsch N, Rolloff KP, Dunn EZ, Friedman JE, Barbour LA, Hernandez TL, MacLean PS, Monks J,
#' McManaman JL, Rudolph MC. Single Cell RNA Sequencing of Human Milk-Derived Cells Reveals Sub-Populations of Mammary Epithelial Cells with Molecular
#' Signatures of Progenitor and Mature States: a Novel, Non-invasive Framework for Investigating Human Lactation Physiology. J Mammary Gland Biol Neoplasia.
#' 2020;25(4):367-387.
NULL

#' single-cell target gene expression data in breast cancer
#'
#' @docType data
#' @name scRNA_tar_cancer
#' @aliases scRNA_tar_cancer
#' @format scRNA_tar_cancer: A data frame object with 736 barcodes of CID3946 breast cancer cells (columns) and 9469 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in breast cancer cells is obtained from Wu et al.
#' The data is pre-processed by using Seurat R package.
#' The matched lncRNA and mRNA expression data is regarded as single-cell target gene expression data in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' single-cell target gene expression data in normal breast
#'
#' @docType data
#' @name scRNA_tar_normal
#' @aliases scRNA_tar_normal
#' @format scRNA_tar_normal: A data frame object with 2430 barcodes of normal breast cells (columns) and 11226 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in normal breast cells is obtained from GEO (GSE15388).
#' The data is pre-processed by using Seurat R package.
#' The matched lncRNA and mRNA expression data is regarded as single-cell target gene expression data in normal breast.
#' @references Martin Carli JF, Trahan GD, Jones KL, Hirsch N, Rolloff KP, Dunn EZ, Friedman JE, Barbour LA, Hernandez TL, MacLean PS, Monks J, McManaman JL,
#' Rudolph MC. Single Cell RNA Sequencing of Human Milk-Derived Cells Reveals Sub-Populations of Mammary Epithelial Cells with Molecular Signatures of Progenitor
#' and Mature States: a Novel, Non-invasive Framework for Investigating Human Lactation Physiology. J Mammary Gland Biol Neoplasia. 2020;25(4):367-387.
NULL

#' single-cell target gene count data in breast cancer
#'
#' @docType data
#' @name scRNA_tar_count_cancer
#' @aliases scRNA_tar_count_cancer
#' @format scRNA_tar_count_cancer: A data frame object with 736 barcodes of CID3946 breast cancer cells (columns) and 9323 target genes (rows).
#' @details The matched lncRNA and mRNA count data is obtained from Wu et al.
#' The count data is pre-processed by using Seurat R package.
#' The matched lncRNA and mRNA expression data is regarded as single-cell target gene count data in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' single-cell ncRNA count data in breast cancer
#'
#' @docType data
#' @name scRNA_ncR_count_cancer
#' @aliases scRNA_ncR_count_cancer
#' @format scRNA_ncR_count_cancer: A data frame object with 736 barcodes of CID3946 breast cancer cells (columns) and 147 ncRNAs (rows).
#' @details The matched miRNA and lncRNA count data is obtained from Wu et al.
#' The count data is pre-processed by using Seurat R package.
#' The matched miRNA and lncRNA expression data is regarded as single-cell ncRNA count data in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' spatial transcriptomics data of lncRNA in breast cancer
#'
#' @docType data
#' @name ST_ncR_cancer
#' @aliases ST_ncR_cancer
#' @format ST_ncR_cancer: A data frame object with 1057 barcodes of CID44971 breast cancer samples (columns) and 2109 lncRNAs (rows).
#' @details The spatial transcriptomics data of lncRNAs in breast cancer is obtained from Wu et al.
#' The spatial transcriptomics data is pre-processed by using Seurat R package.
#' The lncRNA expression data is regarded as spatial transcriptomics data of lncRNA in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' spatial transcriptomics data of mRNA in breast cancer
#'
#' @docType data
#' @name ST_tar_cancer
#' @aliases ST_tar_cancer
#' @format ST_tar_cancer: A data frame object with 1057 barcodes of CID44971 breast cancer samples (columns) and 14801 mRNAs (rows).
#' @details The spatial transcriptomics data of mRNAs in breast cancer is obtained from Wu et al.
#' The spatial transcriptomics data is pre-processed by using Seurat R package.
#' The mRNA expression data is regarded as spatial transcriptomics data of mRNA in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' spatial transcriptomics count data of lncRNA in breast cancer
#'
#' @docType data
#' @name ST_ncR_count_cancer
#' @aliases ST_ncR_count_cancer
#' @format ST_ncR_cancer: A data frame object with 1057 barcodes of CID44971 breast cancer samples (columns) and 2110 lncRNAs (rows).
#' @details The spatial transcriptomics count data of lncRNAs in breast cancer is obtained from Wu et al.
#' The spatial transcriptomics data is pre-processed by using Seurat R package.
#' The lncRNA count data is regarded as spatial transcriptomics count data of lncRNA in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' spatial transcriptomics count data of mRNA in breast cancer
#'
#' @docType data
#' @name ST_tar_count_cancer
#' @aliases ST_tar_count_cancer
#' @format ST_tar_count_cancer: A data frame object with 1057 barcodes of CID44971 breast cancer samples (columns) and 14804 mRNAs (rows).
#' @details The spatial transcriptomics count data of mRNAs in breast cancer is obtained from Wu et al.
#' The spatial transcriptomics data is pre-processed by using Seurat R package.
#' The mRNA count data is regarded as spatial transcriptomics count data of mRNA in breast cancer.
#' @references Wu SZ, Roden DL, Junankar S et al. A single-cell and spatially resolved atlas of human breast cancers. Nat Genet, 2021;53(9):1334-1347.
NULL

#' mutation data
#'
#' @docType data
#' @name mut_ncR
#' @aliases mut_ncR
#' @format mut_ncR: A character object with 158 SNV genes.
#' @details The mutation genes are collected from miRNASNP v3.0 (http://bioinfo.life.hust.edu.cn/miRNASNP/#!/) and lncRNASNP v3.0 (http://bioinfo.life.hust.edu.cn/lncRNASNP#!/) databases.
#' @references Liu CJ, Fu X, Xia M, Zhang Q, Gu Z, Guo AY. miRNASNP-v3: a comprehensive database for SNPs and disease-related variations in miRNAs and miRNA targets. Nucleic Acids Res. 2021;49(D1):D1276-D1281.
#' @references Yang Y, Wang D, Miao YR, Wu X, Luo H, Cao W, Yang W, Yang J, Guo AY, Gong J. lncRNASNP v3: an updated database for functional variants in long non-coding RNAs. Nucleic Acids Res. 2023;51(D1):D192-D198.
#'
NULL

#' BRCA genes
#'
#' @docType data
#' @name BRCA_genes
#' @aliases BRCA_genes
#' @format BRCA_genes: A matrix object with 4819 BRCA related genes (including lncRNAs and mRNAs).
#' @details The BRCA related lncRNAs are from LncRNADisease v2.0, Lnc2Cancer v2.0
#' and MNDR v2.0. The BRCA related mRNAs are from DisGeNET v5.0 and COSMIC v86.
#' @references Bao Z, Yang Z, Huang Z, Zhou Y, Cui Q, Dong D. LncRNADisease 2.0: an updated database of long non-coding RNA-associated diseases. Nucleic Acids Res. 2019;47(D1):D1034-D1037.
#' @references Cui T, Zhang L, Huang Y, Yi Y, Tan P, Zhao Y, Hu Y, Xu L, Li E, Wang D. MNDR v2.0: an updated resource of ncRNA-disease associations in mammals. Nucleic Acids Res. 2018;46(D1):D371-D374.
#' @references Gao Y, Wang P, Wang Y, Ma X, Zhi H, Zhou D, Li X, Fang Y, Shen W, Xu Y, Shang S, Wang L, Wang L, Ning S, Li X. Lnc2Cancer v2.0: updated database of experimentally supported
#' long non-coding RNAs in human cancers. Nucleic Acids Res. 2019;47(D1):D1028-D1033.
#' @references Forbes SA, Beare D, Boutselakis H, Bamford S, Bindal N, Tate J, Cole CG, Ward S, Dawson E, Ponting L, Stefancsik R, Harsha B, Kok CY, Jia M, Jubb H, Sondka Z, Thompson S, De T,
#' Campbell PJ. COSMIC: somatic cancer genetics at high-resolution. Nucleic Acids Res. 2017;45(D1):D777-D783.
#' @references Piñero J, Bravo À, Queralt-Rosinach N, Gutiérrez-Sacristán A, Deu-Pons J, Centeno E, García-García J, Sanz F, Furlong LI. DisGeNET: a comprehensive platform integrating information
#' on human disease-associated genes and variants. Nucleic Acids Res. 2017;45(D1):D833-D839.
#' @references Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. miRSM: an R package
#' to infer and analyse miRNA sponge modules in heterogeneous data. RNA Biol. 2021;18(12):2308-2320.
NULL

#' network
#'
#' @docType data
#' @name network
#' @aliases network
#' @format network: A list object with five ncRNA-target regulatory networks predicted by the SCS method
#' @details A list of five ncRNA-target regulatory networks as a demo
NULL

#' hub_ncRNAs
#'
#' @docType data
#' @name hub_ncRNAs
#' @aliases hub_ncRNAs
#' @format hub_ncRNAs: A list object with hub ncRNAs existed in five ncRNA-target regulatory networks predicted by the SCS method
#' @details A list of hub ncRNAs existed in five ncRNA-target regulatory networks as a demo
NULL


