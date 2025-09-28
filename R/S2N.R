#' @useDynLib S2N, .registration = TRUE
NULL

#' Sample-specific network inferred with the SSN (single-sample network) method
#'
#' @title SSN
#' @param reg_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network
#' @import utils
#' @import stats
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#' @examples
#'
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SSN_res <- SSN(bulk_ncR_normal[1:100, 1:5], bulk_tar_normal[1:500, 1:5], bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SSN_res <- SSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5], scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' @references
#' Liu X, Wang Y, Ji H, Aihara K, Chen L. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016;44(22):e164.
#'
#' @export
SSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer, cormethod = "pearson",
    p.value.cutoff = 0.05) {

    if(is(reg_normal, "SummarizedExperiment") ||
       is(reg_normal, "SingleCellExperiment") ||
       is(reg_normal, "SpatialExperiment")){
      reg_normal <- as.matrix(assay(reg_normal))
      tar_normal <- as.matrix(assay(tar_normal))
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    valid_reg <- intersect(rownames(reg_normal), rownames(reg_cancer))
    valid_reg <- valid_reg[!is.na(valid_reg)]
    valid_tar <- intersect(rownames(tar_normal), rownames(tar_cancer))
    valid_tar <- valid_tar[!is.na(valid_tar)]
    all_pairs <- expand.grid(regulator = valid_reg, target = valid_tar, stringsAsFactors = FALSE)

    pre_r1_info <- lapply(seq_len(nrow(all_pairs)), function(i) {

        reg <- all_pairs$regulator[i]
        tar <- all_pairs$target[i]

        x <- as.numeric(reg_normal[reg, ])
        y <- as.numeric(tar_normal[tar, ])
        r1 <- cor(x, y, method = cormethod)

        list(regulator = reg, target = tar, r1 = r1, x_norm = x, y_norm = y)

    })

    samples_num <- ncol(reg_cancer)
    nn <- ncol(reg_normal)

    SSN_net <- lapply(seq_len(samples_num), function(k) {
        cat("Processing sample", k, "/", samples_num, "\n")

        res <- do.call(rbind, lapply(pre_r1_info, function(info) {
            x_cancer <- as.numeric(reg_cancer[info$regulator, k])
            y_cancer <- as.numeric(tar_cancer[info$target, k])

            combined_x <- c(info$x_norm, x_cancer)
            combined_y <- c(info$y_norm, y_cancer)
            r2 <- cor(combined_x, combined_y, method = cormethod)

            delta <- r2 - info$r1
            z <- ssn_zscore(delta, info$r1, nn)

            if (is.na(z)) {
                p <- NA_real_
            } else if (z < 0) {
                p <- pnorm(z)
            } else {
                p <- 1 - pnorm(z)
            }

            c(info$regulator, info$target, info$r1, r2, delta, z, p)

        }))
    })

    SSN_net <- lapply(seq(SSN_net), function(i) SSN_net[[i]][which(as.numeric(SSN_net[[i]][,
        7]) < p.value.cutoff), seq(2)])

    SSN_net_graph <- lapply(seq(SSN_net), function(i) make_graph(c(t(SSN_net[[i]])), directed = FALSE))
    names(SSN_net_graph) <- colnames(reg_cancer)

    return(SSN_net_graph)
}

#' Sample-specific network inferred with the Paired-SSN (paired single sample network) method
#'
#' @title PairedSSN
#' @param reg_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network
#' @import utils
#' @import stats
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_PaiedSSN_res <- PairedSSN(bulk_ncR_normal[1:100, 1:5], bulk_tar_normal[1:500, 1:5], bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_PaiedSSN_res <- PairedSSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5], scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' @references
#' Guo WF, Zhang SW, Zeng T, Li Y, Gao J, Chen L. A novel network control model for identifying personalized driver genes in cancer. PLoS Comput Biol. 2019;15(11):e1007520.
#'
#' @export
PairedSSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer, cormethod = "pearson",
    p.value.cutoff = 0.05) {

   if(is(reg_normal, "SummarizedExperiment") ||
      is(reg_normal, "SingleCellExperiment") ||
      is(reg_normal, "SpatialExperiment")){
     reg_normal <- as.matrix(assay(reg_normal))
     tar_normal <- as.matrix(assay(tar_normal))
     reg_cancer <- as.matrix(assay(reg_cancer))
     tar_cancer <- as.matrix(assay(tar_cancer))
   }

    reg_intersect <- intersect(rownames(reg_normal), rownames(reg_cancer))
    reg_intersect <- reg_intersect[!is.na(reg_intersect)]
    tar_intersect <- intersect(rownames(tar_normal), rownames(tar_cancer))
    tar_intersect <- tar_intersect[!is.na(tar_intersect)]

    reg_normal <- reg_normal[reg_intersect, , drop = FALSE]
    tar_normal <- tar_normal[tar_intersect, , drop = FALSE]
    reg_cancer <- reg_cancer[reg_intersect, , drop = FALSE]
    tar_cancer <- tar_cancer[tar_intersect, , drop = FALSE]

    all_pairs <- expand.grid(regulator = reg_intersect, target = tar_intersect, stringsAsFactors = FALSE)

    PairedSSN_net <- lapply(seq_len(ncol(reg_cancer)), function(k) {
        cat("Processing sample", k, "/", ncol(reg_cancer), "\n")

        curr_reg_cancer <- reg_cancer[, k, drop = FALSE]
        curr_tar_cancer <- tar_cancer[, k, drop = FALSE]
        curr_reg_normal <- reg_normal[, k, drop = FALSE]
        curr_tar_normal <- tar_normal[, k, drop = FALSE]

        ssn_cancer <- SSN_paired(reg_normal, tar_normal, curr_reg_cancer, curr_tar_cancer,
            cormethod, nn = ncol(reg_normal))
        ssn_normal <- SSN_paired(reg_normal, tar_normal, curr_reg_normal, curr_tar_normal,
            cormethod, nn = ncol(reg_normal))

        process_p <- function(pmat) {
            pmat[is.na(pmat)] <- 1
            pmat[pmat >= p.value.cutoff] <- 0
            pmat[pmat != 0] <- 1
            pmat
        }

        diff_matrix <- abs(process_p(ssn_cancer$p) - process_p(ssn_normal$p))

        sig_indices <- which(diff_matrix > 0, arr.ind = TRUE)
        data.frame(regulator = rownames(diff_matrix)[sig_indices[, 1]], target = colnames(diff_matrix)[sig_indices[,
            2]], stringsAsFactors = FALSE)
    })

    PairedSSN_net_graph <- lapply(seq(PairedSSN_net), function(i) make_graph(c(t(PairedSSN_net[[i]])), directed = FALSE))
    names(PairedSSN_net_graph) <- colnames(reg_cancer)

    return(PairedSSN_net_graph)
}

#' Sample-specific network inferred with the LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples) method
#'
#' @title LIONESS
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @import stats
#' @import utils
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#' @references
#' Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019;14:226-240.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_LIONESS_res <- LIONESS(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_LIONESS_res <- LIONESS(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#'
#' @export
LIONESS <- function(reg_cancer, tar_cancer, cormethod = "pearson") {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    N <- ncol(reg_cancer)
    PCC <- cor(x = t(reg_cancer), y = t(tar_cancer), method = cormethod)
    LIONESS_net <- list()

    for (k in seq(N)) {
        cat("Processing sample", k, "/", N, "\n")

        PCC1 <- cor(x = t(reg_cancer[, -k]), y = t(tar_cancer[, -k]), method = cormethod)
        x <- N * (PCC - PCC1) + PCC1
        x[is.na(x)] <- 0

        xx <- x
        xx[lower.tri(xx, diag = FALSE)] <- 0
        a_x <- abs(xx[xx != 0])
        threshold <- mean(a_x, na.rm = T) + 2 * sd(a_x, na.rm = T)

        rm(xx)
        gc()

        x[abs(x) < threshold] <- 0
        x[x != 0] <- 1
        index <- which(x != 0, arr.ind = TRUE)
        LIONESS_net[[k]] <- data.frame(regulator = rownames(x)[index[, 1]], target = colnames(x)[index[,
            2]])
    }

    LIONESS_net_graph <- lapply(seq(LIONESS_net), function(i) make_graph(c(t(LIONESS_net[[i]])), directed = FALSE))
    names(LIONESS_net_graph) <- colnames(reg_cancer)

    return(LIONESS_net_graph)
}

#' Sample-specific network inferred with SCS (single-sample controller strategy) method
#'
#' @title SCS
#' @param mut_reg A matrix or data frame object, containing a list of mutation regulators.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param number  The number of RWR with random matrix.
#' @import stats
#' @import igraph
#' @import dplyr
#' @import Matrix
#' @import methods
#' @import utils
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#' @references
#' Guo WF, Zhang SW, Liu LL, Liu F, Shi QQ, Zhang L, Tang Y, Zeng T, Chen L. Discovering personalized driver mutation profiles of single samples in cancer by network control strategy. Bioinformatics. 2018;34(11):1893-1903.
#'
#' @examples
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' data(mut_ncR)
#' scRNA_SCS_res <- SCS(mut_ncR, scRNA_ncR_cancer[, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = 'pearson', number = 10)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' data(mut_ncR)
#' ST_SCS_res <- SCS(mut_ncR, ST_ncR_cancer[, 1:5], ST_tar_cancer[1:500, 1:5], cormethod = 'pearson', number = 10)
#'
#' @export
#'
SCS <- function(mut_reg, reg_cancer, tar_cancer, cormethod = "pearson", number = 50) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    edge0 <- expand.grid(regulator = mut_reg, target = rownames(tar_cancer))
    colnames(edge0) <- c("regulator", "target")
    node0 <- unique(c(edge0$regulator, edge0$target))

    SNP_name <- mut_reg
    reg_EXPR_data <- reg_cancer
    reg_EXPR_name <- rownames(reg_cancer)
    EXPR_data <- tar_cancer
    EXPR_name <- rownames(tar_cancer)
    original_sample_names <- colnames(reg_cancer)
    rm("mut_reg", "reg_cancer", "tar_cancer")
    gc()

    z1 <- match(edge0[, 1], node0)
    z2 <- match(edge0[, 2], node0)
    z <- as.data.frame(cbind(z1, z2))
    z_genes <- as.data.frame(cbind(node0[z1], node0[z2]))
    z <- dplyr::distinct(z)
    z_genes <- dplyr::distinct(z_genes)

    omiga <- 1
    Targets_genes_sampels <- list()
    Mutation_genes_sampels <- list()
    for (i in 1:ncol(EXPR_data)) {
        target_genes <- EXPR_name[which((abs(EXPR_data[, i]) > omiga) %in% TRUE)]
        Targets_genes_sampels[[i]] <- target_genes

        SNP_genes_tmp <- c()
        if (length(abs(reg_EXPR_data[, i]) != 0) != 0) {
            SNP_genes_tmp <- reg_EXPR_name[which(abs(reg_EXPR_data[, i]) != 0)]
        }
        Mutation_gene <- intersect(SNP_name, SNP_genes_tmp)
        Mutation_genes_sampels[[i]] <- Mutation_gene
    }

    index_delete_sample <- which(sapply(Targets_genes_sampels, length) * sapply(Mutation_genes_sampels,
        length) == 0)
    if (length(index_delete_sample) > 0) {
      kept_sample_indices <- setdiff(1:length(Targets_genes_sampels), index_delete_sample)
      kept_sample_names <- original_sample_names[kept_sample_indices]

      Targets_genes_sampels <- Targets_genes_sampels[-index_delete_sample]
      Mutation_genes_sampels <- Mutation_genes_sampels[-index_delete_sample]
    } else {
      kept_sample_indices <- 1:length(Targets_genes_sampels)
      kept_sample_names <- original_sample_names
    }

    Targets_C_genes <- lapply(Targets_genes_sampels, function(x) intersect(x, node0))
    Constrained_B_genes <- lapply(Mutation_genes_sampels, function(x) intersect(x,
        node0))

    SCS_net <- list()
    for (i in seq_along(Constrained_B_genes)) {
        cat("Processing sample", i, "/", length(Constrained_B_genes), "\n")

        B_fda_genes <- Constrained_B_genes[[i]]
        C_genes <- Targets_C_genes[[i]]

        edge0_B_C <- edge0[edge0$regulator %in% B_fda_genes & edge0$target %in% C_genes,
            ]
        if (nrow(edge0_B_C) == 0)
	      next
        edge0_network <- igraph::graph_from_data_frame(edge0_B_C, directed = TRUE)
        edge0_MultiplexObject <- create.multiplex(list(edge0_network), directed = TRUE)
        AdjMatrix_edge0 <- compute.adjacency.matrix(edge0_MultiplexObject)

        RWR_results_edge0 <- RWR_MN(AdjMatrix_edge0, edge0_MultiplexObject, B_fda_genes,
            scores = 1e-04)
        colnames(RWR_results_edge0)[3] <- "gdtruth_score"

        RWR_results_dir_srand <- lapply(1:number, function(j) {
            dir_srand <- dir_generate_srand(as.matrix(AdjMatrix_edge0))
            dir_srand_AdjMatrix <- as(as.matrix(dir_srand), "dgCMatrix")
            RWR_MN(dir_srand_AdjMatrix, edge0_MultiplexObject, B_fda_genes, scores = 1e-04)
        })

        RWR_results_dir_srand_merge <- Reduce(function(x, y) merge(x, y, by = c("SeedGenes",
            "NodeNames")), RWR_results_dir_srand)
        RWR_results_edge0_dir_srand <- merge(RWR_results_edge0, RWR_results_dir_srand_merge)

        p_value <- sapply(1:nrow(RWR_results_edge0_dir_srand), function(k) {
            rand_scores <- as.numeric(RWR_results_edge0_dir_srand[k, -(1:3)])
            if (all(rand_scores <= RWR_results_edge0_dir_srand$gdtruth_score[k])) {
                return(RWR_results_edge0_dir_srand$gdtruth_score[k])
            } else {
                return(0)
            }
        })

        SCS_net[[i]] <- RWR_results_edge0_dir_srand[p_value > 0, 1:2]
    }

    SCS_net_graph <- lapply(seq(SCS_net), function(i) make_graph(c(t(SCS_net[[i]])), directed = FALSE))
    names(SCS_net_graph) <- kept_sample_names

    return(SCS_net_graph)
}


#' Sample-specific network inferred with the CSN (cell-specific network) method, adapted from the CSmiR method
#'
#' @title CSN
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param boxsize Size of neighborhood (0.1 in default).
#' @param bootstrap A logic value, TRUE for bootstrap. If the number of samples is less than 100, boostrap is recommended.
#' @param bootstrap_betw_point The number of interpolation points between each cell (5 in default), bootstrap_betw_point = 0 is used to compute.
#' @param bootstrap_num The number of bootstrapping for interpolating pseudo-cells.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network.
#' @param num.cores The number of cores to run in parallel.
#' @import igraph
#' @import stats
#' @import utils
#' @import pracma
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import SummarizedExperiment
#' @importFrom WGCNA pmedian
#' @return An igraph object: a list of regulator-target interactions for each sample.
#' @references
#' Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.2021;22(1):578.
#' Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019;47(11):e62.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_CSN_res <- CSN(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], boxsize = 0.1, p.value.cutoff = 0.05)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_CSN_res <- CSN(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], boxsize = 0.1, p.value.cutoff = 0.05)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_CSN_res <- CSN(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], boxsize = 0.1, p.value.cutoff = 0.05)
#'
#' @export
CSN <- function(reg_cancer, tar_cancer, boxsize = 0.1, bootstrap = FALSE, bootstrap_betw_point = 5,
    bootstrap_num = 100, p.value.cutoff = 0.05, num.cores = 2) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    reg_cancer <- data.frame(reg_cancer)
    tar_cancer <- data.frame(tar_cancer)

    reg_data_cancer <- t(reg_cancer)
    tar_data_cancer <- t(tar_cancer)
    rm("reg_cancer", "tar_cancer")
    gc()

    if (bootstrap) {
        CSN_net <- CSN_net_bootstrap(reg_data_cancer, tar_data_cancer, boxsize = boxsize,
            bootstrap_betw_point = bootstrap_betw_point, bootstrap_num = bootstrap_num,
            p.value.cutoff = p.value.cutoff)
    } else {
        CSN_net <- CSN_net(reg_data_cancer, tar_data_cancer, boxsize = boxsize, p.value.cutoff = p.value.cutoff,
            num.cores = num.cores)
    }

    CSN_net_graph <- lapply(seq(CSN_net), function(i) make_graph(c(t(CSN_net[[i]])), directed = FALSE))
    names(CSN_net_graph) <- rownames(reg_data_cancer)

    return(CSN_net_graph)
}

#' Identifying sample-specific ncRNA regulation using Scan (Sample-speCific miRNA regulAtioN) and statistical perturbation strategy
#'
#' @title Scan.perturb
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param method Methods for calculating correlations, one of 'pearson', 'spearman', 'kendall'.
#' @param p.value.cutoff Significance level for the identified regulator-target interactions.
#' @param num.cores Number of CPU cores when parallel computation.
#' @importFrom igraph as_data_frame
#' @importFrom igraph graph_from_biadjacency_matrix
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom WGCNA corAndPvalue
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Zhang J, Liu L, Wei X, Zhao C, Luo Y, Li J, Le TD. Scanning sample-specific miRNA regulation from bulk and single-cell RNA-sequencing data. BMC Biol. 2024;22(1):218.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_Scan.perturb_res <- Scan.perturb(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_Scan.perturb_res <- Scan.perturb(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_Scan.perturb_res <- Scan.perturb(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' @export
Scan.perturb <- function(reg_cancer, tar_cancer, method = c("pearson", "spearman",
    "kendall"), p.value.cutoff = 0.05, num.cores = 2) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    if (method == "pearson") {
        res.all <- Pearson_adj(reg_cancer, tar_cancer, p.value.cutoff)

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Pearson_adj")) %dopar% {
            Pearson_adj(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
        }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "spearman") {
        res.all <- Spearman_adj(reg_cancer, tar_cancer, p.value.cutoff)

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Spearman_adj")) %dopar% {
            Spearman_adj(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
        }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "kendall") {
        res.all <- Kendall_adj(reg_cancer, tar_cancer, p.value.cutoff)

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Kendall_adj")) %dopar% {
            Kendall_adj(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
        }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()
    }

    Scan_perturb_incident <- lapply(seq(res.single), function(i) abs(res.single[[i]] -
        res.all))

    Scan_perturb_graph <- lapply(seq(res.single), function(i) graph_from_biadjacency_matrix(t(Scan_perturb_incident[[i]])))
    names(Scan_perturb_graph) <- colnames(reg_cancer)

    return(Scan_perturb_graph)
}

#' Identifying sample-specific ncRNA regulation using Scan (Sample-speCific miRNA regulAtioN) and linear interpolation strategy
#'
#' @title Scan.interp
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param method Methods for calculating correlations, one of 'pearson', 'spearman', 'kendall'.
#' @param p.value.cutoff Significance level for the identified regulator-target interactions.
#' @param num.cores Number of CPU cores when parallel computation.
#' @importFrom igraph as_data_frame
#' @importFrom igraph graph_from_biadjacency_matrix
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import SummarizedExperiment
#' @importFrom WGCNA cor
#' @importFrom WGCNA corAndPvalue#'
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Zhang J, Liu L, Wei X, Zhao C, Luo Y, Li J, Le TD. Scanning sample-specific miRNA regulation from bulk and single-cell RNA-sequencing data. BMC Biol. 2024;22(1):218.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_Scan.interp_res <- Scan.interp(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_Scan.interp_res <- Scan.interp(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_Scan.interp_res <- Scan.interp(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], method = 'pearson', p.value.cutoff = 0.05)
#'
#' @export
Scan.interp <- function(reg_cancer, tar_cancer, method = c("pearson", "spearman",
    "kendall"), p.value.cutoff = 0.05, num.cores = 2) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    reg_cancer <- t(reg_cancer)
    tar_cancer <- t(tar_cancer)
    nsamples <- nrow(reg_cancer)

    if (method == "pearson") {

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(nrow(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Pearson_Scan", "Pearson", "matrixzscore")) %dopar%
            {
                Pearson_Scan(reg_cancer, tar_cancer, i)
            }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "spearman") {

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(nrow(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Spearman_Scan", "Spearman", "matrixzscore")) %dopar%
            {
                Spearman_Scan(reg_cancer, tar_cancer, i)
            }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    } else if (method == "kendall") {

        # get number of cores to run
        cl <- makeCluster(num.cores)
        registerDoParallel(cl)

        res.single <- foreach(i = seq(nrow(reg_cancer)), .packages = c("igraph",
            "WGCNA"), .export = c("Kendall_Scan", "Kendall", "matrixzscore")) %dopar%
            {
                Kendall_Scan(reg_cancer, tar_cancer, i)
            }

        # shut down the workers
        stopCluster(cl)
        stopImplicitCluster()

    }

    Scan_interp_incident <- lapply(seq(res.single), function(i) ifelse(res.single[[i]] <
        p.value.cutoff, 1, 0))

    Scan_interp_graph <- lapply(seq(res.single), function(i) graph_from_biadjacency_matrix(t(Scan_interp_incident[[i]])))
    names(Scan_interp_graph) <- rownames(reg_cancer)

    return(Scan_interp_graph)
}

#' Sample-specific network inferred with the SWEET (Sample-specific weighted correlation network) method
#'
#' @title SWEET
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param dismethod A character string indicating which type of distance is to be computed. One of 'euclidean' (default), or 'correlation' can be used.
#' @param threshold The cutoff of correlation strength.
#' @import stats
#' @import SummarizedExperiment
#' @importFrom igraph graph_from_biadjacency_matrix
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Chen HH, Hsueh CW, Lee CH, Hao TY, Tu TY, Chang LY, Lee JC, Lin CY. SWEET: a single-sample network inference method for deciphering individual features in disease. Brief Bioinform. 2023, 24(2):bbad032.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SWEET_res <- SWEET(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], dismethod = 'euclidean')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SWEET_res <- SWEET(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], dismethod = 'euclidean')
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_SWEET_res <- SWEET(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], dismethod = 'euclidean')
#'
#' @export
SWEET <- function(reg_cancer, tar_cancer, dismethod = "euclidean", threshold = 0.7) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

   # calculate distances between samples
   expr <- rbind(reg_cancer, tar_cancer)
   n_genes <- nrow(expr)
   n_samples <- ncol(expr)
   gene_names <- rownames(expr)
   sample_names <- colnames(expr)
   weight_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
   rownames(weight_matrix) <- colnames(expr)
   colnames(weight_matrix) <- colnames(expr)

   if (dismethod == "euclidean") {
    dist_matrix <- as.matrix(dist(t(expr), method = "euclidean"))
   } else if (dismethod == "correlation") {
    dist_matrix <- 1 - cor(expr, use = "pairwise.complete.obs")
   }

   sigma <- median(dist_matrix, na.rm = TRUE)
   weight_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))

   diag(weight_matrix) <- 0

   # calculate sample-specific correlation
   sample_specific_cors <- array(0, dim = c(n_genes, n_genes, n_samples))
   dimnames(sample_specific_cors) <- list(gene_names, gene_names, sample_names)

   for (s in 1:n_samples) {
     weights <- weight_matrix[s, ]
     weighted_expr <- t(expr) * weights
     cov_matrix <- cov.wt(t(expr), wt = weights, cor = TRUE)$cor
     sample_specific_cors[, , s] <- cov_matrix
   }

   # extract sample specific network
   SWEET_net <- list()
   for (s in 1:n_samples) {
      cor_matrix <- sample_specific_cors[, , s]
      adj_matrix <- ifelse(abs(cor_matrix) > threshold, 1, 0)
      SWEET_net[[sample_names[s]]] <- adj_matrix[1:nrow(reg_cancer), (1+nrow(reg_cancer)):(nrow(reg_cancer)+nrow(tar_cancer))]
   }

   SWEET_net_graph <- lapply(seq(SWEET_net), function(i) graph_from_biadjacency_matrix(SWEET_net[[i]]))
   names(SWEET_net_graph) <- colnames(reg_cancer)

   return(SWEET_net_graph)
}

#' Sample-specific network inferred with the SINUM (SIngle-cell Network Using Mutual information) method
#'
#' @title SINUM
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param boxsize Size of neighborhood (0.2 in default).
#' @param zscore.cutoff Zscore cutoff for identifying sample-specific regulatory network.
#' @param num.cores The number of cores to run in parallel.
#' @import utils
#' @import stats
#' @import Matrix
#' @import entropy
#' @import infotheo
#' @import foreach
#' @import doParallel
#' @import igraph
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Chang LY, Hao TY, Wang WJ, Lin CY. Inference of single-cell network using mutual information for scRNA-seq data analysis. BMC Bioinformatics. 2024;25(Suppl 2):292.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SINUM_res <- SINUM(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], boxsize = 0.2, zscore.cutoff = 0)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SINUM_res <- SINUM(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], boxsize = 0.2, zscore.cutoff = 0)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_SINUM_res <- SINUM(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], boxsize = 0.2, zscore.cutoff = 0)
#'
#' @export
SINUM <- function(reg_cancer, tar_cancer, boxsize = 0.2, zscore.cutoff = 0, num.cores = 2) {

    if(is(reg_cancer, "SummarizedExperiment") ||
       is(reg_cancer, "SingleCellExperiment") ||
       is(reg_cancer, "SpatialExperiment")){
      reg_cancer <- as.matrix(assay(reg_cancer))
      tar_cancer <- as.matrix(assay(tar_cancer))
    }

    # Combine regulatory and target cancer data
    GEM <- rbind(reg_cancer, tar_cancer)

    # Get gene and cell names
    geneLst <- rownames(GEM)
    cellLst <- colnames(GEM)

    # Convert to matrix format
    data <- as.matrix(GEM)

    # Get number of genes and cells
    geneTot <- nrow(data)
    cellTot <- ncol(data)
    RegTot <- nrow(reg_cancer)
    TarTot <- nrow(tar_cancer)

    # Calculate number of bins for discretization (using square root rule)
    binsNum <- round(sqrt(cellTot))

    # Discretize the data
    Xt <- discretize_data(data, binsNum)

    # Calculate upper and lower bounds
    ul <- get_upper_and_lower_file(data, boxsize)
    upper <- ul$upper
    lower <- ul$lower

    # Discretize upper and lower bounds using the same function
    Xupp <- discretize_data(upper, binsNum) + 1  # Add 1 to avoid 0 values
    Xlow <- discretize_data(lower, binsNum)

    # Calculate entropy matrix for each gene
    EntropyX <- get_entropyX_matrix(Xt, geneTot, cellTot, binsNum)

    # Generate all possible regulatory-target gene pair combinations
    genePairs <- expand.grid(1:RegTot, (RegTot + 1):geneTot)
    genePairnum <- nrow(genePairs)  # Number of gene pairs

    # Initialize mutual information matrix (gene pairs × cells)
    mutual_info <- matrix(0, nrow = genePairnum, ncol = cellTot)

    # Set up parallel computing
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    # Calculate mutual information for all gene pairs in parallel
    mutual_info <- foreach(s = 1:genePairnum, .combine = rbind, .packages = c("entropy", "infotheo"), .export = c("get_entropyXY_matrix")) %dopar% {
      # Get indices of current gene pair
      g1 <- genePairs[s, 1]
      g2 <- genePairs[s, 2]

      # Calculate joint entropy for the two genes
      entropies <- get_entropyXY_matrix(Xt, g1, g2, cellTot, binsNum)

      # Initialize mutual information vector for current gene pair across all cells
      mi_cell <- numeric(cellTot)

      # Calculate local mutual information for each cell
      for (cell in 1:cellTot) {
        # If either gene is not expressed in current cell, mutual information is 0
        if (data[g1, cell] * data[g2, cell] == 0) {
          mi_cell[cell] <- 0
          next
        }

        # Get upper and lower bounds for both genes in current cell
        # Add bounds checking to prevent "subscript out of bounds" error
        start1 <- max(1, min(binsNum, Xlow[g1, cell]))
        end1 <- max(1, min(binsNum, Xupp[g1, cell]))
        start2 <- max(1, min(binsNum, Xlow[g2, cell]))
        end2 <- max(1, min(binsNum, Xupp[g2, cell]))

        # Ensure start is not greater than end
        if (start1 > end1) {
          start1 <- 1
          end1 <- binsNum
        }
        if (start2 > end2) {
          start2 <- 1
          end2 <- binsNum
        }

        # Calculate local entropy and joint entropy
        HX <- -sum(EntropyX[g1, start1:end1])  # Local entropy of gene 1
        HY <- -sum(EntropyX[g2, start2:end2])  # Local entropy of gene 2
        HXY <- -sum(entropies[start1:end1, start2:end2])  # Local joint entropy of gene pair

        # Calculate mutual information: I(X;Y) = H(X) + H(Y) - H(X,Y)
        mi_cell[cell] <- round(HX + HY - HXY, 6)
      }

      # Return mutual information for current gene pair across all cells
      mi_cell
    }

    # Stop parallel computing cluster
    stopCluster(cl)

    # Convert mutual information to Z-scores
    zscore <- convert_zScore(mutual_info)

    # Initialize list for adjacency matrices
    adj_matrix <- list()

    # Build adjacency matrix for each cell
    for (cell in 1:cellTot) {
      cat(sprintf('Processing cell %d/%d...\n', cell, cellTot))

      # Set adjacency relationships based on Z-score threshold
      # Z-score > threshold indicates gene pair correlation, set to 1, otherwise 0
      interin <- ifelse(zscore[, cell] > zscore.cutoff, 1, 0)

      # Reshape to regulatory × target matrix
      adj_matrix[[cell]] <- matrix(interin, nrow = RegTot, ncol = TarTot, byrow = FALSE)
      rownames(adj_matrix[[cell]]) <- rownames(reg_cancer)
      colnames(adj_matrix[[cell]]) <- rownames(tar_cancer)
    }

    # Convert adjacency matrices to bipartite graphs
    SINUM_net_graph <- lapply(adj_matrix, function(x) {
      graph_from_biadjacency_matrix(x, directed = FALSE, weighted = NULL)
    })

    names(SINUM_net_graph) <- colnames(reg_cancer)

    return(SINUM_net_graph)
}

#' Sample-specific network inferred with the BONOBO (Bayesian Optimized Networks Obtained By assimilating Omics data) method
#'
#' @title BONOBO
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param p.value.cutoff pvalue cutoff for identifying sample-specific regulatory network.
#' @import stats
#' @import igraph
#' @import SummarizedExperiment
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Saha E, Fanfani V, Mandros P, Ben Guebila M, Fischer J, Shutta KH, DeMeo DL, Lopes-Ramos CM, Quackenbush J. Bayesian inference of sample-specific coexpression networks. Genome Res. 2024;34(9):1397-1410.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_BONOBO_res <- BONOBO(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5])
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_BONOBO_res <- BONOBO(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5])
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_BONOBO_res <- BONOBO(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5])
#'
#' @export
BONOBO <- function(reg_cancer, tar_cancer, p.value.cutoff = 0.05) {

  if(is(reg_cancer, "SummarizedExperiment") ||
     is(reg_cancer, "SingleCellExperiment") ||
     is(reg_cancer, "SpatialExperiment")){
    reg_cancer <- as.matrix(assay(reg_cancer))
    tar_cancer <- as.matrix(assay(tar_cancer))
  }

  # Number of regulators and targets
  RegTot <- nrow(reg_cancer)
  TarTot <- nrow(tar_cancer)

  # Combine regulatory and target cancer data
  expression_matrix <- as.matrix(rbind(reg_cancer, tar_cancer))
  expression_mean <- rowMeans(expression_matrix)

  # Create a matrix with expression_mean values repeated for each column
  expression_mean_matrix <- matrix(expression_mean, nrow = nrow(expression_matrix),
                                   ncol = ncol(expression_matrix))
  expression_diff <- expression_matrix - expression_mean_matrix

  mask_include <- rep(TRUE, ncol(expression_matrix))
  BONOBO_net_graph <- list()

  for (sample_idx in seq(ncol(expression_matrix))){

    # Create mask to exclude the sample of interest
    mask_include[sample_idx] <- FALSE

    # Compute covariance matrix from the remaining data (excluding the target sample)
    covariance_matrix <- stats::cov(t(expression_matrix[, mask_include]))

    # Check for NA/NaN in covariance matrix
    if (any(is.na(covariance_matrix) | any(is.nan(covariance_matrix)))) {
      warning("Covariance matrix contains NA/NaN values. Adding regularization.")
      covariance_matrix <- covariance_matrix + diag(1e-10, nrow(covariance_matrix))
      covariance_matrix[is.na(covariance_matrix) | is.nan(covariance_matrix)] <- 1e-10
    }

    # Check for constant genes (zero variance)
    diag_cov <- diag(covariance_matrix)

    zero_variance_genes <- which(diag_cov == 0)
    if (length(zero_variance_genes) > 0) {
      warning(length(zero_variance_genes), " genes have zero variance. Adding small regularization.")
      diag(covariance_matrix)[zero_variance_genes] <- 1e-10
    }

    # Compute posterior weight delta from data
    delta <- 1 / (3 + 2 * mean(sqrt(diag(covariance_matrix))) / var(diag(covariance_matrix)))

    # Compute sample-specific covariance matrix
    sscov <- delta * tcrossprod(expression_diff[, sample_idx]) +
      (1 - delta) * covariance_matrix

    # Infer sample-specific network
    g <- ncol(sscov)
    d <- g + 1 / delta
    a1 <- (d - g + 1) / ((d - g) * (d - g - 3))
    a2 <- (d - g - 1) / ((d - g) * (d - g - 3))
    v <- diag(sscov)
    # Calculate pointwise variance
    v_matrix <- a1 * (sscov * sscov) + a2 * tcrossprod(v)
    # Calculate pointwise standard deviation
    sd_matrix <- sqrt(v_matrix)

    # Check for zero or negative values in sd_matrix
    sd_matrix[sd_matrix <= 0 | is.na(sd_matrix) | is.nan(sd_matrix)] <- 1e-10

    # Calculate z-scores
    z_matrix <- sscov / sd_matrix

    # Calculate p-values
    pval <- 2 * (1 - stats::pnorm(abs(z_matrix)))
    adj_matrix <- ifelse(pval < p.value.cutoff, 1, 0)
    adj_matrix_sub <- adj_matrix[1:RegTot, (RegTot + 1):(RegTot + TarTot)]

    BONOBO_net_graph[[sample_idx]] <- graph_from_biadjacency_matrix(adj_matrix_sub, directed = FALSE, weighted = NULL)
  }

  names(BONOBO_net_graph) <- colnames(reg_cancer)

  return(BONOBO_net_graph)
}

#' Validating sample-specific network
#'
#' @title Network_validation
#' @param network_graph An igraph object: a list of sample-specific networks.
#' @param groundtruth A data frame object, experimentally validated interactions.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param num.cores Number of CPU cores when parallel computation.
#' @import igraph
#' @import foreach
#' @return A list object: a list of validated interactions for each sample.
#'
#' @export
Network_validation <- function(network_graph, groundtruth, directed = FALSE, num.cores = 2) {

    groundtruth_graph <- make_graph(c(t(groundtruth[, c(1, 2)])), directed = directed)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    network_validated <- foreach(i = seq(network_graph), .packages = c("igraph", "foreach",
        "doParallel")) %dopar% {
        as_data_frame(network_graph[[i]] %s% groundtruth_graph)
    }  # transform to a data frame

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

    return(network_validated)
}

#' Identifying sample-specific hubs
#'
#' @title Network_hub
#' @param network_graph An igraph object: a list of of sample-specific networks.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param top_percentage A number indicating the top percenatge of hub regulators, default 0.2.
#' @import igraph
#' @return A list object: a list of hubs for each sample.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' hubs <- Network_hub(bulk_LIONESS_res)
#'
#' @references
#'Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.2021; 22(1):578.
#'
#' @export
Network_hub <- function(network_graph, directed = FALSE, top_percentage = 0.2) {

    if (directed == TRUE) {
        network_degree <- lapply(seq(network_graph), function(i) igraph::degree(network_graph[[i]],
            mode = "out"))
    } else {
        network_degree <- lapply(seq(network_graph), function(i) igraph::degree(network_graph[[i]]))
    }
    network_hubs <- lapply(seq(network_graph), function(i) names(sort(network_degree[[i]][which(network_degree[[i]] !=
        0)], decreasing = TRUE))[1:ceiling(top_percentage * length(which(network_degree[[i]] !=
        0)))])

    names(network_hubs) <- names(network_graph)

    return(network_hubs)
}

#' Predicting sample-specific modules
#'
#' @title Network_module
#' @param network_graph An igraph object, a list of sample-specific networks.
#' @param method Possible methods include FN, MCL, LINKCOMM, MCODE, betweenness, infomap, prop, eigen, louvain, walktrap.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param modulesize The size of modules, e.g. 3.
#' @import igraph
#' @importFrom MCL mcl
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.flush
#' @importFrom grDevices dev.hold
#' @importFrom graphics abline
#' @importFrom graphics barplot
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics points
#' @importFrom graphics polygon
#' @importFrom graphics text
#' @return A list object: a list of sample-specific modules.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' Network_module_res <- Network_module(bulk_LIONESS_res)
#'
#' @export
Network_module <- function(network_graph, method = "MCL", directed = FALSE, modulesize = 3) {

    network_Cluster <- list()
    edgelist <- list()
    network_Cluster_result <- list()
    size <- list()

    for (k in seq(network_graph)) {

        if (method == "FN" | method == "MCL" | method == "MCODE") {
            network_Cluster[[k]] <- cluster(network_graph[[k]],
                method = method, directed = directed)
        } else if (method == "LINKCOMM") {
            edgelist[[k]] <- get.edgelist(network_graph[[k]])
            network_Cluster[[k]] <- getLinkCommunities(edgelist[[k]], directed = directed)$nodeclusters
        }

        if (method == "FN" | method == "MCL") {
            network_Cluster_result[[k]] <- lapply(seq_len(max(network_Cluster[[k]])),
                function(i) rownames(as.matrix(network_Cluster[[k]]))[which(network_Cluster[[k]] ==
                  i)])
            size[[k]] <- unlist(lapply(seq_len(max(network_Cluster[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "LINKCOMM") {
            network_Cluster_result[[k]] <- lapply(seq_len(max(c(network_Cluster[[k]]$cluster))),
                function(i) as.character(network_Cluster[[k]]$node[which(c(network_Cluster[[k]]$cluster) ==
                  i)]))
            size[[k]] <- unlist(lapply(seq_len(max(c(network_Cluster[[k]]$cluster))),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "MCODE") {
            network_Cluster[[k]] <- network_Cluster[[k]] + 1
            network_Cluster_result[[k]] <- lapply(seq_len(max(network_Cluster[[k]])),
                function(i) rownames(as.matrix(network_Cluster[[k]]))[which(network_Cluster[[k]] ==
                  i)])
            size[[k]] <- unlist(lapply(seq_len(max(network_Cluster[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "betweenness") {
            network_Cluster_result[[k]] <- cluster_edge_betweenness(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "infomap") {
            network_Cluster_result[[k]] <- cluster_infomap(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "prop") {
            network_Cluster_result[[k]] <- cluster_label_prop(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "eigen") {
            network_Cluster_result[[k]] <- cluster_leading_eigen(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "louvain") {
            network_Cluster_result[[k]] <- cluster_louvain(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        } else if (method == "walktrap") {
            network_Cluster_result[[k]] <- cluster_walktrap(network_graph[[k]])
            size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])),
                function(i) length(network_Cluster_result[[k]][[i]])))
            network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize),
                function(i) network_Cluster_result[[k]][[i]])
        }
    }

    names(network_Cluster_result) <- names(network_graph)

    return(network_Cluster_result)

}

#' Calculating sample-sample similarity matrix between two lists of sample-specific networks
#'
#' @title Sim.network
#' @param net1 The first list of sample-specific networks, an igraph object.
#' @param net2 The second list of sample-specific networks, an igraph object.
#' @param method Methods for calculating similatiry between two two lists of sample-specific networks, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @import igraph
#' @return A matrix object, a sample-sample similarity matrix in terms of sample-specific networks.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' res <- Sim.network(bulk_LIONESS_res,  bulk_LIONESS_res)
#'
#' @export
Sim.network <- function(net1, net2, method = "Simpson") {

    if (class(net1) != "list" | class(net2) != "list") {
        stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)

    if (method == "Simpson") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- ecount(net1[[i]] %s% net2[[j]])
                Sim[i, j] <- overlap_interin/min(ecount(net1[[i]]), ecount(net2[[j]]))
            }
        }
    } else if (method == "Jaccard") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- ecount(net1[[i]] %s% net2[[j]])
                Sim[i, j] <- overlap_interin/ecount(net1[[i]] %u% net2[[j]])
            }
        }
    } else if (method == "Lin") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- ecount(net1[[i]] %s% net2[[j]])
                Sim[i, j] <- 2 * overlap_interin/(ecount(net1[[i]]) + ecount(net2[[j]]))
            }
        }
    }

    rownames(Sim) <- names(net1)
    colnames(Sim) <- names(net2)

    return(Sim)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific networks
#'
#' @title Network_sim_hc
#' @param network_graph An igraph object, a list of sample-specific networks.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param method Methods for calculating similatiry between two two lists of sample-specific networks, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @import igraph
#' @return Hierarchical cluster analysis result of samples.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' res <- Network_sim_hc(bulk_LIONESS_res)
#' plot(res)
#'
#' @export
Network_sim_hc <- function(network_graph, directed = FALSE, method = "Simpson") {

    network_Sim <- Sim.network(network_graph, network_graph, method = method)
    hclust_res_network <- hclust(as.dist(1 - network_Sim), "complete")

    return(hclust_res_network)
}


#' Calculating sample-sample similarity matrix between two lists of sample-specific hubs
#'
#' @title Sim.hub
#' @param hub1 The first list of sample-specific hubs.
#' @param hub2 The second list of sample-specific hubs.
#' @param method Methods for calculating similarity between two two lists of sample-specific hubs, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return A matrix object, a sample-sample similarity matrix in terms of sample-specific hubs.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' hubs <- Network_hub(bulk_LIONESS_res)
#' res <- Sim.hub(hubs, hubs)
#'
#' @export
Sim.hub <- function(hub1, hub2, method = "Simpson") {

    if (class(hub1) != "list" | class(hub2) != "list") {
        stop("Please check your input hub! The input hub should be list object! \n")
    }
    m <- length(hub1)
    n <- length(hub2)
    Sim <- matrix(NA, m, n)

    if (method == "Simpson") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
                Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
            }
        }
    } else if (method == "Jaccard") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
                Sim[i, j] <- overlap_interin/length(union(hub1[[i]], hub2[[j]]))
            }
        }
    } else if (method == "Lin") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
                Sim[i, j] <- 2 * overlap_interin/(length(hub1[[i]]) + length(hub2[[i]]))
            }
        }
    }

    rownames(Sim) <- names(hub1)
    colnames(Sim) <- names(hub2)

    return(Sim)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific hubs
#'
#' @title Hub_sim_hc
#' @param hub A list object, a list of sample-specific hubs.
#' @param method Methods for calculating similatiry between two two lists of sample-specific hubs, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Hierarchical cluster analysis result of samples.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' hubs <- Network_hub(bulk_LIONESS_res)
#' res <- Hub_sim_hc(hubs)
#' plot(res)
#'
#' @export
Hub_sim_hc <- function(hub, method = "Simpson") {

    hub_Sim <- Sim.hub(hub, hub, method = method)
    hclust_res_hub <- hclust(as.dist(1 - hub_Sim), "complete")

    return(hclust_res_hub)
}

#' Calculating sample-sample similarity matrix between two lists of sample-specific modules
#'
#' @title Sim.module
#' @param Module1 A list object, the first list of modules.
#' @param Module2 A list object, the second list of modules.
#' @param method Methods for calculating similatiry between two modules, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return A matrix object, a sample-sample similarity matrix in terms of sample-specific modules.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' Network_module_res <- Network_module(bulk_LIONESS_res)
#' res <- Sim.module(Network_module_res[[1]], Network_module_res[[2]])
#'
#' @export
Sim.module <- function(Module1, Module2, method = "Simpson") {

    if (class(Module1) != "list" | class(Module2) != "list") {
        stop("Please check your input modules! The input modules should be list object! \n")
    }

    m <- length(Module1)
    n <- length(Module2)
    Sim <- matrix(NA, m, n)

    if (method == "Simpson") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
                min_vertex <- min(length(Module1[[i]]), length(Module2[[j]]))
                Sim[i, j] <- overlap_vertex/min_vertex
            }
        }
    } else if (method == "Jaccard") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
                union_vertex <- length(union(Module1[[i]], Module2[[j]]))
                Sim[i, j] <- overlap_vertex/union_vertex
            }
        }
    } else if (method == "Lin") {
        for (i in seq(m)) {
            for (j in seq(n)) {
                overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
                sum_vertex <- length(Module1[[i]]) + length(Module2[[j]])
                Sim[i, j] <- 2 * overlap_vertex/sum_vertex
            }
        }
    }

    if (m < n) {
        GS <- mean(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]]))) *
            m/n
    } else if (m == n) {
        GS <- mean(c(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])),
            unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i]))))
    } else if (m > n) {
        GS <- mean(unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i]))) *
            n/m
    }

    return(GS)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific modules
#'
#' @title Module_sim_hc
#' @param module A list object, a list of sample-specific module groups. Each module group contain a list of modules.
#' @param method Methods for calculating similarity between two lists of sample-specific modules, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Hierarchical cluster analysis result of samples.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' Network_module_res <- Network_module(bulk_LIONESS_res)
#' res <- Module_sim_hc(Network_module_res)
#' plot(res)
#'
#' @export
Module_sim_hc <- function(module, method = "Simpson") {

    n <- length(module)
    module_Sim <- matrix(NA, n, n)
    for (i in seq(n)) {
        for (j in seq(n)) {
            module_Sim[i, j] <- Sim.module(module[[i]], module[[j]], method = method)
        }
    }
    rownames(module_Sim) <- colnames(module_Sim) <- names(module)
    hclust_res_module <- hclust(as.dist(1 - module_Sim), "complete")

    return(hclust_res_module)
}
