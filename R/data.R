#' Bulk ncRNA expression data in breast cancer
#'
#' @docType data
#' @name bulk_ncR_cancer
#' @aliases bulk_ncR_cancer
#' @format bulk_ncR_cancer: A data frame object with six BRCA samples (columns) and 4794 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in breast cancer is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' Bulk ncRNA expression data in normal breast
#'
#' @docType data
#' @name bulk_ncR_normal
#' @aliases bulk_ncR_normal
#' @format bulk_ncR_normal: A data frame object with six normal samples (columns) and 4929 ncRNAs (rows).
#' @details The matched miRNA and lncRNA expression data in normal breast is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' Bulk target gene expression data in breast cancer
#'
#' @docType data
#' @name bulk_tar_cancer
#' @aliases bulk_tar_cancer
#' @format bulk_tar_cancer: A data frame object with six BRCA samples (columns) and 22427 target genes (rows).
#' @details The matched lncRNA and mRNA expression data in breast cancer is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' Bulk target gene expression data in normal breast
#'
#' @docType data
#' @name bulk_tar_normal
#' @aliases bulk_tar_normal
#' @format bulk_tar_normal: A dataframe object with six normal samples (columns) and 22427 target genes (rows).
#' @details The matched lncRNA and mRNA expression data in normal breast is obtained from TCGA (http://cancergenome.nih.gov/).
NULL

#' Single-cell ncRNA expression data in breast cancer
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

#' Single-cell ncRNA expression data in normal breast
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

#' Single-cell target gene expression data in breast cancer
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

#' Single-cell target gene expression data in normal breast
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

#' Single-cell target gene count data in breast cancer
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

#' Single-cell ncRNA count data in breast cancer
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

#' Spatial transcriptomics data of lncRNA in breast cancer
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

#' Spatial transcriptomics data of mRNA in breast cancer
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

#' Spatial transcriptomics count data of lncRNA in breast cancer
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

#' Spatial transcriptomics count data of mRNA in breast cancer
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

#' Mutation data
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

#' Groundtruth of miRNA-target interactions
#'
#' @docType data
#' @name miRTarget_Groundtruth
#' @aliases miRTarget_Groundtruth
#' @format miRTarget_Groundtruth: A data frame object with 2,956,530 miRNA-target interactions.
#' @details The groundtruth of miRNA-target interactions are from miRTarBase v10.0 and TarBase v9.0.
#' @references Cui S, Yu S, Huang HY, Lin YC, Huang Y, Zhang B, Xiao J, Zuo H, Wang J, Li Z, Li G, Ma J, Chen B, Zhang H, Fu J, Wang L,
#' Huang HD. miRTarBase 2025: updates to the collection of experimentally validated microRNA-target interactions. Nucleic Acids Res. 2025;53(D1):D147-D156.
#' @references Skoufos G, Kakoulidis P, Tastsoglou S, Zacharopoulou E, Kotsira V, Miliotis M, Mavromati G, Grigoriadis D, Zioga M, Velli A, Koutou I,
#' Karagkouni D, Stavropoulos S, Kardaras FS, Lifousi A, Vavalou E, Ovsepian A, Skoulakis A, Tasoulis SK, Georgakopoulos SV, Plagianakos VP, Hatzigeorgiou AG.
#' TarBase-v9.0 extends experimentally supported miRNA-gene interactions to cell-types and virally encoded miRNAs. Nucleic Acids Res. 2024;52(D1):D304-D310.
NULL

#' Groundtruth of lncRNA-target interactions
#'
#' @docType data
#' @name lncRTarget_Groundtruth
#' @aliases lncRTarget_Groundtruth
#' @format lncRTarget_Groundtruth: A data frame object with 3,039,533 lncRNA-target interactions.
#' @details The groundtruth of lncRNA-target interactions are from RegNetwork v2025.
#' @references Li B, Wang C, Wang Y, Li P, Liu ZP. RegNetwork 2025: an integrative data repository for gene regulatory networks in human and mouse.
#' Nucleic Acids Res. 2025:gkaf779. doi: 10.1093/nar/gkaf779.
NULL

#' Groundtruth of TF-target interactions
#'
#' @docType data
#' @name TFTarget_Groundtruth
#' @aliases TFTarget_Groundtruth
#' @format TFTarget_Groundtruth: A data frame object with 43,178 TF-target interactions.
#' @details The groundtruth of TF-target interactions are from CollecTRI.
#' @references Müller-Dott S, Tsirvouli E, Vazquez M, Ramirez Flores RO, Badia-I-Mompel P, Fallegger R, Türei D, Lægreid A, Saez-Rodriguez J.
#' Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. Nucleic Acids Res. 2023 Nov 10;51(20):10934-10949.
NULL
