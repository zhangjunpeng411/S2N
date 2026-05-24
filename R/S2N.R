#' @useDynLib S2N, .registration = TRUE
NULL

#' Sample-specific network inferred with the SSN (single-sample network) method
#'
#' @title SSN
#' @param reg_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed.
#'   One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import utils
#' @import stats
#' @import SummarizedExperiment
#' @import parallel
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SSN_res <- SSN(bulk_ncR_normal[1:100, 1:5], bulk_tar_normal[1:500, 1:5],
#'                     bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5],
#'                     cormethod = 'pearson')
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SSN_res <- SSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5],
#'                      scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5],
#'                      cormethod = 'pearson')
#'
#' @references
#' Liu X, Wang Y, Ji H, Aihara K, Chen L. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016;44(22):e164.
#'
#' @export
SSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer,
                cormethod = "pearson", p.value.cutoff = 0.05,
                num.cores = 2) {

    if (is(reg_normal, "SummarizedExperiment") ||
        is(reg_normal, "SingleCellExperiment") ||
        is(reg_normal, "SpatialExperiment")) {
        reg_normal <- as.matrix(assay(reg_normal))
        tar_normal <- as.matrix(assay(tar_normal))
        reg_cancer <- as.matrix(assay(reg_cancer))
        tar_cancer <- as.matrix(assay(tar_cancer))
    }

    # Common regulators and targets
    valid_reg <- intersect(rownames(reg_normal), rownames(reg_cancer))
    valid_reg <- valid_reg[!is.na(valid_reg)]
    valid_tar <- intersect(rownames(tar_normal), rownames(tar_cancer))
    valid_tar <- valid_tar[!is.na(valid_tar)]
    all_pairs <- expand.grid(regulator = valid_reg, target = valid_tar, stringsAsFactors = FALSE)

    # Pre-compute normal correlations and expression vectors for each pair
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

    # Function to process one cancer sample
    process_sample <- function(k) {
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
        # Filter by p-value and return only regulator-target pairs (columns 1:2)
        res[which(as.numeric(res[, 7]) < p.value.cutoff), seq(2), drop = FALSE]
    }

    # Execute serially or in parallel
    if (num.cores == 1) {
        SSN_net <- lapply(seq_len(samples_num), function(k) {
            cat("Processing sample", k, "/", samples_num, "\n")
            process_sample(k)
        })
    } else {
        if (.Platform$OS.type == "windows") {
            cl <- parallel::makeCluster(num.cores)
            on.exit(parallel::stopCluster(cl))
            parallel::clusterExport(cl, varlist = c("process_sample", "pre_r1_info",
                                                    "reg_cancer", "tar_cancer",
                                                    "cormethod", "nn", "p.value.cutoff",
                                                    "ssn_zscore", "pnorm"),
                                    envir = environment())
            SSN_net <- parallel::parLapply(cl, seq_len(samples_num), process_sample)
        } else {
            SSN_net <- parallel::mclapply(seq_len(samples_num), process_sample,
                                          mc.cores = num.cores, mc.preschedule = TRUE)
        }
    }

    # Build igraph objects
    SSN_net_graph <- lapply(SSN_net, function(edges) {
        if (nrow(edges) > 0) {
            igraph::make_graph(c(t(edges)), directed = TRUE)
        } else {
            NULL
        }
    })
    names(SSN_net_graph) <- colnames(reg_cancer)

    return(SSN_net_graph)
}

#' Sample-specific network inferred with the Paired-SSN (paired single sample network) method
#'
#' @title PairedSSN
#' @param reg_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed.
#'   One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import utils
#' @import stats
#' @import SummarizedExperiment
#' @import parallel
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_PaiedSSN_res <- PairedSSN(bulk_ncR_normal[1:100, 1:5], bulk_tar_normal[1:500, 1:5],
#'                                bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5],
#'                                cormethod = 'pearson')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_PaiedSSN_res <- PairedSSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5],
#'                                 scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5],
#'                                 cormethod = 'pearson')
#'
#' @references
#' Guo WF, Zhang SW, Zeng T, Li Y, Gao J, Chen L. A novel network control model for identifying personalized driver genes in cancer. PLoS Comput Biol. 2019;15(11):e1007520.
#'
#' @export
PairedSSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer,
                      cormethod = "pearson", p.value.cutoff = 0.05,
                      num.cores = 2) {

    if (is(reg_normal, "SummarizedExperiment") ||
        is(reg_normal, "SingleCellExperiment") ||
        is(reg_normal, "SpatialExperiment")) {
        reg_normal <- as.matrix(assay(reg_normal))
        tar_normal <- as.matrix(assay(tar_normal))
        reg_cancer <- as.matrix(assay(reg_cancer))
        tar_cancer <- as.matrix(assay(tar_cancer))
    }

    # Intersect features
    reg_intersect <- intersect(rownames(reg_normal), rownames(reg_cancer))
    reg_intersect <- reg_intersect[!is.na(reg_intersect)]
    tar_intersect <- intersect(rownames(tar_normal), rownames(tar_cancer))
    tar_intersect <- tar_intersect[!is.na(tar_intersect)]

    reg_normal <- reg_normal[reg_intersect, , drop = FALSE]
    tar_normal <- tar_normal[tar_intersect, , drop = FALSE]
    reg_cancer <- reg_cancer[reg_intersect, , drop = FALSE]
    tar_cancer <- tar_cancer[tar_intersect, , drop = FALSE]

    all_pairs <- expand.grid(regulator = reg_intersect, target = tar_intersect,
                             stringsAsFactors = FALSE)
    samples_num <- ncol(reg_cancer)
    nn <- ncol(reg_normal)

    # Internal function to process one sample (cancer and corresponding normal)
    process_sample <- function(k) {
        curr_reg_cancer <- reg_cancer[, k, drop = FALSE]
        curr_tar_cancer <- tar_cancer[, k, drop = FALSE]
        curr_reg_normal <- reg_normal[, k, drop = FALSE]
        curr_tar_normal <- tar_normal[, k, drop = FALSE]

        ssn_cancer <- SSN_paired(reg_normal, tar_normal, curr_reg_cancer, curr_tar_cancer,
                                 cormethod, nn = nn)
        ssn_normal <- SSN_paired(reg_normal, tar_normal, curr_reg_normal, curr_tar_normal,
                                 cormethod, nn = nn)

        # Convert p-value matrix to binary (0/1) based on cutoff
        process_p <- function(pmat) {
            pmat[is.na(pmat)] <- 1
            pmat[pmat >= p.value.cutoff] <- 0
            pmat[pmat != 0] <- 1
            pmat
        }

        diff_matrix <- abs(process_p(ssn_cancer$p) - process_p(ssn_normal$p))
        sig_indices <- which(diff_matrix > 0, arr.ind = TRUE)
        data.frame(regulator = rownames(diff_matrix)[sig_indices[, 1]],
                   target = colnames(diff_matrix)[sig_indices[, 2]],
                   stringsAsFactors = FALSE)
    }

    # Execute serially or in parallel
    if (num.cores == 1) {
        PairedSSN_net <- lapply(seq_len(samples_num), function(k) {
            cat("Processing sample", k, "/", samples_num, "\n")
            process_sample(k)
        })
    } else {
        if (.Platform$OS.type == "windows") {
            cl <- parallel::makeCluster(num.cores)
            on.exit(parallel::stopCluster(cl))
            parallel::clusterExport(cl, varlist = c("process_sample", "reg_normal",
                                                    "tar_normal", "reg_cancer",
                                                    "tar_cancer", "cormethod",
                                                    "p.value.cutoff", "nn",
                                                    "SSN_paired"),
                                    envir = environment())
            PairedSSN_net <- parallel::parLapply(cl, seq_len(samples_num), function(k) {
                process_sample(k)
            })
        } else {
            PairedSSN_net <- parallel::mclapply(seq_len(samples_num), function(k) {
                process_sample(k)
            }, mc.cores = num.cores, mc.preschedule = TRUE)
        }
    }

    # Build igraph objects
    PairedSSN_net_graph <- lapply(seq_along(PairedSSN_net), function(i) {
        edges <- PairedSSN_net[[i]]
        if (nrow(edges) > 0) {
            igraph::make_graph(c(t(edges)), directed = TRUE)
        } else {
            NULL
        }
    })
    names(PairedSSN_net_graph) <- colnames(reg_cancer)

    return(PairedSSN_net_graph)
}

#' Sample-specific network inferred with the LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples) method
#'
#' @title LIONESS
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed.
#'   One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import stats
#' @import utils
#' @import SummarizedExperiment
#' @import parallel
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019;14:226-240.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5],
#'                             cormethod = 'pearson')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_LIONESS_res <- LIONESS(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5],
#'                              cormethod = 'pearson')
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_LIONESS_res <- LIONESS(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5],
#'                           cormethod = 'pearson')
#'
#' @export
LIONESS <- function(reg_cancer, tar_cancer, cormethod = "pearson", num.cores = 2) {

    if (is(reg_cancer, "SummarizedExperiment") ||
        is(reg_cancer, "SingleCellExperiment") ||
        is(reg_cancer, "SpatialExperiment")) {
        reg_cancer <- as.matrix(assay(reg_cancer))
        tar_cancer <- as.matrix(assay(tar_cancer))
    }

    N <- ncol(reg_cancer)
    PCC <- cor(x = t(reg_cancer), y = t(tar_cancer), method = cormethod)
    sample_names <- colnames(reg_cancer)

    # Function to process one sample
    process_sample <- function(k) {
        # Compute correlation without sample k
        PCC1 <- cor(x = t(reg_cancer[, -k, drop = FALSE]),
                    y = t(tar_cancer[, -k, drop = FALSE]),
                    method = cormethod)
        x <- N * (PCC - PCC1) + PCC1
        x[is.na(x)] <- 0

        # Thresholding
        xx <- x
        xx[lower.tri(xx, diag = FALSE)] <- 0
        a_x <- abs(xx[xx != 0])
        threshold <- mean(a_x, na.rm = TRUE) + 2 * sd(a_x, na.rm = TRUE)

        x[abs(x) < threshold] <- 0
        x[x != 0] <- 1
        index <- which(x != 0, arr.ind = TRUE)

        if (nrow(index) == 0) {
            return(data.frame(regulator = character(0), target = character(0),
                              stringsAsFactors = FALSE))
        }
        data.frame(regulator = rownames(x)[index[, 1]],
                   target = colnames(x)[index[, 2]],
                   stringsAsFactors = FALSE)
    }

    # Execute serially or in parallel
    if (num.cores == 1) {
        LIONESS_net <- lapply(seq_len(N), function(k) {
            cat("Processing sample", k, "/", N, "\n")
            process_sample(k)
        })
    } else {
        if (.Platform$OS.type == "windows") {
            cl <- parallel::makeCluster(num.cores)
            on.exit(parallel::stopCluster(cl))
            parallel::clusterExport(cl, varlist = c("process_sample", "reg_cancer",
                                                    "tar_cancer", "cormethod",
                                                    "N", "PCC"),
                                    envir = environment())
            LIONESS_net <- parallel::parLapply(cl, seq_len(N), function(k) {
                process_sample(k)
            })
        } else {
            LIONESS_net <- parallel::mclapply(seq_len(N), function(k) {
                process_sample(k)
            }, mc.cores = num.cores, mc.preschedule = TRUE)
        }
    }

    # Build igraph objects
    LIONESS_net_graph <- lapply(seq_along(LIONESS_net), function(i) {
        edges <- LIONESS_net[[i]]
        if (nrow(edges) > 0) {
            igraph::make_graph(c(t(edges)), directed = TRUE)
        } else {
            NULL
        }
    })
    names(LIONESS_net_graph) <- sample_names

    return(LIONESS_net_graph)
}

#' Sample-specific network inferred with SCS (single-sample controller strategy) method
#'
#' @title SCS
#' @param mut_reg A matrix or data frame object, containing a list of mutation regulators.
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed.
#'   One of 'pearson' (default), 'kendall', or 'spearman' can be used.
#' @param number The number of RWR with random matrix.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import stats
#' @import igraph
#' @import dplyr
#' @import Matrix
#' @import methods
#' @import utils
#' @import SummarizedExperiment
#' @import parallel
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Guo WF, Zhang SW, Liu LL, Liu F, Shi QQ, Zhang L, Tang Y, Zeng T, Chen L. Discovering personalized driver mutation profiles of single samples in cancer by network control strategy. Bioinformatics. 2018;34(11):1893-1903.
#'
#' @examples
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' data(mut_ncR)
#' scRNA_SCS_res <- SCS(mut_ncR, scRNA_ncR_cancer[, 1:5], scRNA_tar_cancer[1:500, 1:5],
#'                      cormethod = 'pearson', number = 10)
#' # ST data
#' data(ST_ncR_tar_cancer)
#' data(mut_ncR)
#' ST_SCS_res <- SCS(mut_ncR, ST_ncR_cancer[, 1:5], ST_tar_cancer[1:500, 1:5],
#'                   cormethod = 'pearson', number = 10)
#'
#' @export
SCS <- function(mut_reg, reg_cancer, tar_cancer, cormethod = "pearson", number = 50, num.cores = 2) {

  if (is(reg_cancer, "SummarizedExperiment") ||
      is(reg_cancer, "SingleCellExperiment") ||
      is(reg_cancer, "SpatialExperiment")) {
    reg_cancer <- as.matrix(assay(reg_cancer))
    tar_cancer <- as.matrix(assay(tar_cancer))
  }

  # Initialize edge list and node set
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

  # Deduplicate edges
  z1 <- match(edge0[, 1], node0)
  z2 <- match(edge0[, 2], node0)
  z <- as.data.frame(cbind(z1, z2))
  z_genes <- as.data.frame(cbind(node0[z1], node0[z2]))
  z <- dplyr::distinct(z)
  z_genes <- dplyr::distinct(z_genes)

  # Identify active target genes and mutated regulators per sample
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

  # Remove samples with no target or no mutated regulators
  index_delete_sample <- which(sapply(Targets_genes_sampels, length) * sapply(Mutation_genes_sampels, length) == 0)
  if (length(index_delete_sample) > 0) {
    kept_sample_indices <- setdiff(1:length(Targets_genes_sampels), index_delete_sample)
    kept_sample_names <- original_sample_names[kept_sample_indices]
    Targets_genes_sampels <- Targets_genes_sampels[-index_delete_sample]
    Mutation_genes_sampels <- Mutation_genes_sampels[-index_delete_sample]
  } else {
    kept_sample_indices <- 1:length(Targets_genes_sampels)
    kept_sample_names <- original_sample_names
  }

  # Intersect with network nodes
  Targets_C_genes <- lapply(Targets_genes_sampels, function(x) intersect(x, node0))
  Constrained_B_genes <- lapply(Mutation_genes_sampels, function(x) intersect(x, node0))

  # Function to process a single sample
  process_sample <- function(i, B_list, C_list, edge0, node0, number) {
    B_fda_genes <- B_list[[i]]
    C_genes <- C_list[[i]]
    edge0_B_C <- edge0[edge0$regulator %in% B_fda_genes & edge0$target %in% C_genes, ]
    if (nrow(edge0_B_C) == 0) return(NULL)

    edge0_network <- igraph::graph_from_data_frame(edge0_B_C, directed = TRUE)
    edge0_MultiplexObject <- create.multiplex(list(edge0_network), directed = TRUE)
    AdjMatrix_edge0 <- compute.adjacency.matrix(edge0_MultiplexObject)

    # True RWR scores
    RWR_results_edge0 <- RWR_MN(AdjMatrix_edge0, edge0_MultiplexObject, B_fda_genes, scores = 1e-04)
    colnames(RWR_results_edge0)[3] <- "gdtruth_score"

    # Random matrix RWR for significance
    RWR_results_dir_srand <- lapply(1:number, function(j) {
      dir_srand <- dir_generate_srand(as.matrix(AdjMatrix_edge0))
      dir_srand_AdjMatrix <- as(as.matrix(dir_srand), "dgCMatrix")
      RWR_MN(dir_srand_AdjMatrix, edge0_MultiplexObject, B_fda_genes, scores = 1e-04)
    })

    RWR_results_dir_srand_merge <- Reduce(function(x, y) merge(x, y, by = c("SeedGenes", "NodeNames")),
                                          RWR_results_dir_srand)
    RWR_results_edge0_dir_srand <- merge(RWR_results_edge0, RWR_results_dir_srand_merge)

    # Compute p-values
    p_value <- sapply(1:nrow(RWR_results_edge0_dir_srand), function(k) {
      rand_scores <- as.numeric(RWR_results_edge0_dir_srand[k, -(1:3)])
      if (all(rand_scores <= RWR_results_edge0_dir_srand$gdtruth_score[k])) {
        return(RWR_results_edge0_dir_srand$gdtruth_score[k])
      } else {
        return(0)
      }
    })

    edges <- RWR_results_edge0_dir_srand[p_value > 0, 1:2]
    if (nrow(edges) == 0) return(NULL)
    return(list(index = i, edges = edges))
  }

  # Execute either serial or parallel
  if (num.cores == 1) {
    res_list <- list()
    for (i in seq_along(Constrained_B_genes)) {
      cat("Processing sample", i, "/", length(Constrained_B_genes), "\n")
      res <- process_sample(i, Constrained_B_genes, Targets_C_genes, edge0, node0, number)
      if (!is.null(res)) res_list[[length(res_list) + 1]] <- res
    }
  } else {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(num.cores)
      on.exit(parallel::stopCluster(cl))
      # Load required packages on each worker
      parallel::clusterEvalQ(cl, {
        library(igraph)
        library(dplyr)
        library(Matrix)
      })
      # Export only variables and functions (not package names)
      parallel::clusterExport(cl, varlist = c("process_sample", "Constrained_B_genes",
                                              "Targets_C_genes", "edge0", "node0", "number",
                                              "RWR_MN", "dir_generate_srand", "create.multiplex",
                                              "compute.adjacency.matrix"),
                              envir = environment())
      res_list <- parallel::parLapply(cl, seq_along(Constrained_B_genes), function(i) {
        process_sample(i, Constrained_B_genes, Targets_C_genes, edge0, node0, number)
      })
    } else {
      res_list <- parallel::mclapply(seq_along(Constrained_B_genes), function(i) {
        process_sample(i, Constrained_B_genes, Targets_C_genes, edge0, node0, number)
      }, mc.cores = num.cores, mc.preschedule = TRUE)
    }
    res_list <- res_list[!sapply(res_list, is.null)]
  }

  # If no results, return empty list
  if (length(res_list) == 0) return(list())

  # Order results by original sample index
  indices <- sapply(res_list, function(x) x$index)
  order_idx <- order(indices)
  res_list <- res_list[order_idx]

  # Build igraph objects and assign sample names
  SCS_net_graph <- list()
  for (k in seq_along(res_list)) {
    edges <- res_list[[k]]$edges
    if (nrow(edges) > 0) {
      SCS_net_graph[[k]] <- igraph::make_graph(c(t(edges)), directed = TRUE)
    }
  }
  names(SCS_net_graph) <- kept_sample_names[indices[order_idx]]

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
    bootstrap_num = 10, p.value.cutoff = 0.05, num.cores = 2) {

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

    # Build igraph objects
    CSN_net_graph <- lapply(seq(CSN_net), function(i) make_graph(c(t(CSN_net[[i]])), directed = TRUE))
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

    Scan_perturb_incident <- lapply(Scan_perturb_incident, function(mat) {
        mat[is.na(mat)] <- 0
        return(mat)
    })

    # Build igraph objects
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
#' @importFrom WGCNA corAndPvalue
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

    Scan_interp_incident <- lapply(Scan_interp_incident, function(mat) {
        mat[is.na(mat)] <- 0
        return(mat)
    })

    # Build igraph objects
    Scan_interp_graph <- lapply(seq(res.single), function(i) graph_from_biadjacency_matrix(t(Scan_interp_incident[[i]])))
    names(Scan_interp_graph) <- rownames(reg_cancer)

    return(Scan_interp_graph)
}

#' Sample-specific network inferred with the SWEET (Sample-specific weighted correlation network) method
#'
#' @title SWEET
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param dismethod A character string indicating which type of distance is to be computed.
#'   One of 'euclidean' (default), or 'correlation' can be used.
#' @param threshold The cutoff of correlation strength.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import stats
#' @import SummarizedExperiment
#' @import parallel
#' @importFrom igraph graph_from_biadjacency_matrix
#' @return An igraph object: a list of regulator-target interactions for each sample.
#'
#' @references
#' Chen HH, Hsueh CW, Lee CH, Hao TY, Tu TY, Chang LY, Lee JC, Lin CY. SWEET: a single-sample network inference method for deciphering individual features in disease. Brief Bioinform. 2023, 24(2):bbad032.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SWEET_res <- SWEET(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5],
#'                         dismethod = 'euclidean')
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SWEET_res <- SWEET(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5],
#'                          dismethod = 'euclidean')
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_SWEET_res <- SWEET(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5],
#'                       dismethod = 'euclidean')
#'
#' @export
SWEET <- function(reg_cancer, tar_cancer, dismethod = "euclidean",
                  threshold = 0.7, num.cores = 2) {

    if (is(reg_cancer, "SummarizedExperiment") ||
        is(reg_cancer, "SingleCellExperiment") ||
        is(reg_cancer, "SpatialExperiment")) {
        reg_cancer <- as.matrix(assay(reg_cancer))
        tar_cancer <- as.matrix(assay(tar_cancer))
    }

    # Combine regulators and targets
    expr <- rbind(reg_cancer, tar_cancer)
    n_genes <- nrow(expr)
    n_samples <- ncol(expr)
    gene_names <- rownames(expr)
    sample_names <- colnames(expr)

    # Compute sample-wise distance matrix
    if (dismethod == "euclidean") {
        dist_matrix <- as.matrix(dist(t(expr), method = "euclidean"))
    } else if (dismethod == "correlation") {
        dist_matrix <- 1 - cor(expr, use = "pairwise.complete.obs")
    }

    # Convert distance to weight matrix (Gaussian kernel)
    sigma <- median(dist_matrix, na.rm = TRUE)
    weight_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))
    diag(weight_matrix) <- 0

    # Determine regulator and target indices
    n_reg <- nrow(reg_cancer)
    n_tar <- nrow(tar_cancer)
    reg_indices <- 1:n_reg
    tar_indices <- (n_reg + 1):(n_reg + n_tar)

    # Internal function to process one sample
    process_sample <- function(s) {
        weights <- weight_matrix[s, ]
        # Weighted correlation (cov.wt returns correlation matrix)
        cor_mat <- cov.wt(t(expr), wt = weights, cor = TRUE)$cor
        # Extract submatrix and apply threshold
        submat <- cor_mat[reg_indices, tar_indices, drop = FALSE]
        adj <- ifelse(abs(submat) > threshold, 1, 0)
        adj[is.na(adj)] <- 0
        return(adj)
    }

    # Execute serially or in parallel
    if (num.cores == 1) {
        SWEET_net <- lapply(seq_len(n_samples), function(s) {
            cat("Processing sample", s, "/", n_samples, "\n")
            process_sample(s)
        })
    } else {
        if (.Platform$OS.type == "windows") {
            cl <- parallel::makeCluster(num.cores)
            on.exit(parallel::stopCluster(cl))
            parallel::clusterExport(cl, varlist = c("process_sample", "expr",
                                                    "weight_matrix", "reg_indices",
                                                    "tar_indices", "threshold"),
                                    envir = environment())
            SWEET_net <- parallel::parLapply(cl, seq_len(n_samples), function(s) {
                process_sample(s)
            })
        } else {
            SWEET_net <- parallel::mclapply(seq_len(n_samples), function(s) {
                process_sample(s)
            }, mc.cores = num.cores, mc.preschedule = TRUE)
        }
    }

    # Build igraph objects from biadjacency matrices
    SWEET_net_graph <- lapply(SWEET_net, function(mat) {
        if (all(dim(mat) > 0) && any(mat > 0)) {
            igraph::graph_from_biadjacency_matrix(mat)
        } else {
            NULL  # Return NULL for empty graphs to avoid errors
        }
    })
    names(SWEET_net_graph) <- sample_names

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
      graph_from_biadjacency_matrix(x, directed = TRUE, weighted = NULL)
    })

    names(SINUM_net_graph) <- colnames(reg_cancer)

    return(SINUM_net_graph)
}

#' Sample-specific network inferred with the BONOBO (Bayesian Optimized Networks Obtained By assimilating Omics data) method
#'
#' @title BONOBO
#' @param reg_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A dataframe or matrix or SummarizedExperiment or SingleCellExperiment or SpatialExperiment object,
#'   the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param p.value.cutoff pvalue cutoff for identifying sample-specific regulatory network.
#' @param num.cores Integer, number of cores to use for parallel processing. Default is 2.
#'
#' @import stats
#' @import igraph
#' @import SummarizedExperiment
#' @import parallel
#' @return An igraph object: a list of regulator-target interactions for each sample.
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
BONOBO <- function(reg_cancer, tar_cancer, p.value.cutoff = 0.05, num.cores = 2) {

    if (is(reg_cancer, "SummarizedExperiment") ||
        is(reg_cancer, "SingleCellExperiment") ||
        is(reg_cancer, "SpatialExperiment")) {
        reg_cancer <- as.matrix(assay(reg_cancer))
        tar_cancer <- as.matrix(assay(tar_cancer))
    }

    RegTot <- nrow(reg_cancer)
    TarTot <- nrow(tar_cancer)
    expression_matrix <- as.matrix(rbind(reg_cancer, tar_cancer))
    expression_mean <- rowMeans(expression_matrix)
    expression_mean_matrix <- matrix(expression_mean,
                                     nrow = nrow(expression_matrix),
                                     ncol = ncol(expression_matrix))
    expression_diff <- expression_matrix - expression_mean_matrix

    # Pre-compute indices for all samples
    sample_indices <- seq_len(ncol(expression_matrix))

    # Internal function to process one sample
    process_sample <- function(sample_idx) {
        # Exclude current sample from covariance estimation
        mask <- setdiff(sample_indices, sample_idx)
        cov_mat <- stats::cov(t(expression_matrix[, mask, drop = FALSE]))

        # Regularization if necessary
        if (any(is.na(cov_mat) | is.nan(cov_mat))) {
            cov_mat <- cov_mat + diag(1e-10, nrow(cov_mat))
            cov_mat[is.na(cov_mat) | is.nan(cov_mat)] <- 1e-10
        }
        diag_cov <- diag(cov_mat)
        zero_var_genes <- which(diag_cov == 0)
        if (length(zero_var_genes) > 0) {
            diag(cov_mat)[zero_var_genes] <- 1e-10
        }

        # Compute delta
        delta <- 1 / (3 + 2 * mean(sqrt(diag_cov)) / var(diag_cov))

        # Sample-specific covariance
        sscov <- delta * tcrossprod(expression_diff[, sample_idx]) +
            (1 - delta) * cov_mat

        # Compute z-scores and p-values
        g <- ncol(sscov)
        d <- g + 1 / delta
        a1 <- (d - g + 1) / ((d - g) * (d - g - 3))
        a2 <- (d - g - 1) / ((d - g) * (d - g - 3))
        v <- diag(sscov)
        v_matrix <- a1 * (sscov * sscov) + a2 * tcrossprod(v)
        sd_matrix <- sqrt(v_matrix)
        sd_matrix[sd_matrix <= 0 | is.na(sd_matrix) | is.nan(sd_matrix)] <- 1e-10

        z_matrix <- sscov / sd_matrix
        pval <- 2 * (1 - stats::pnorm(abs(z_matrix)))
        adj_matrix <- ifelse(pval < p.value.cutoff, 1, 0)
        adj_matrix_sub <- adj_matrix[1:RegTot, (RegTot + 1):(RegTot + TarTot), drop = FALSE]

        # Return adjacency submatrix
        return(adj_matrix_sub)
    }

    # Execute serially or in parallel
    if (num.cores == 1) {
        BONOBO_adj <- lapply(sample_indices, function(idx) {
            cat("Processing sample", idx, "/", length(sample_indices), "\n")
            process_sample(idx)
        })
    } else {
        if (.Platform$OS.type == "windows") {
            cl <- parallel::makeCluster(num.cores)
            on.exit(parallel::stopCluster(cl))
            parallel::clusterExport(cl, varlist = c("process_sample", "expression_matrix",
                                                    "expression_diff", "RegTot", "TarTot",
                                                    "p.value.cutoff", "sample_indices"),
                                    envir = environment())
            BONOBO_adj <- parallel::parLapply(cl, sample_indices, function(idx) {
                process_sample(idx)
            })
        } else {
            BONOBO_adj <- parallel::mclapply(sample_indices, function(idx) {
                process_sample(idx)
            }, mc.cores = num.cores, mc.preschedule = TRUE)
        }
    }

    # Build igraph objects from biadjacency matrices
    BONOBO_net_graph <- lapply(BONOBO_adj, function(adj) {
        if (any(adj > 0)) {
            igraph::graph_from_biadjacency_matrix(adj, directed = TRUE, weighted = NULL)
        } else {
            # Return an empty graph or NULL? igraph can handle empty adjacency, but we'll return NULL
            NULL
        }
    })
    names(BONOBO_net_graph) <- colnames(reg_cancer)

    return(BONOBO_net_graph)
}

#' Validating sample-specific network
#'
#' @param network_graph List of igraph objects.
#' @param groundtruth Data frame of validated edges (two columns).
#' @param num.cores CPU cores.
#' @param calculate_metrics If TRUE, also returns Accuracy, F1, AUPRC, Validated_Percentage.
#'
#' @import igraph
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @return If calculate_metrics=FALSE, list of igraph intersection graphs.
#'   If TRUE, list with `edges` (igraph objects) and `metrics`.
#'
#' @references
#' Cui S, Yu S, Huang H-Y, Lin Y-C-D, Huang Y, Zhang B, Xiao J, Zuo H, Wang J, Li Z, et al. miRTarBase 2025: updates to the collection of experimentally validated microRNA-target interactions. Nucleic Acids Research. 2025;53(D1):D147–D156.
#' Skoufos G, Kakoulidis P, Tastsoglou S, Zacharopoulou E, Kotsira V, Miliotis M, Mavromati G, Grigoriadis D, Zioga M, Velli A, et al. TarBase-v9.0 extends experimentally supported miRNA-gene interactions to cell-types and virally encoded miRNAs. Nucleic Acids Research. 2024;52(D1):D304–D310.
#' Müller-Dott S, Tsirvouli E, Vazquez M, Ramirez Flores RO, Badia-I-Mompel P, Fallegger R, Türei D, Lægreid A, Saez-Rodriguez J. Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities. Nucleic Acids Research. 2023;51(20):10934–10949.
#' Li B, Wang C, Wang Y, Li P, Liu Z-P. RegNetwork 2025: an integrative data repository for gene regulatory networks in human and mouse. Nucleic Acids Research. 2025 Aug 13:gkaf779.
#'
#' @export
Network_validation <- function(network_graph, groundtruth, num.cores = 2,
                               calculate_metrics = FALSE) {

  cl <- parallel::makeCluster(num.cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  res <- foreach::foreach(g = network_graph, .packages = "igraph") %dopar% {
    # If the network is NULL, return placeholders with NA metrics
    if (is.null(g)) {
      if (!calculate_metrics) {
        return(NULL)  # For edges-only mode, return NULL to keep position
      } else {
        return(list(edges = NULL,
                    metrics = data.frame(Accuracy = NA_real_,
                                         F1 = NA_real_,
                                         AUPRC = NA_real_,
                                         Validated_Percentage = NA_real_)))
      }
    }

    # g is a valid igraph object
    dir_g <- igraph::is_directed(g)
    gt_graph <- igraph::make_graph(c(t(groundtruth[, 1:2])), directed = dir_g)
    overlap_graph <- g %s% gt_graph

    if (!calculate_metrics) {
      return(overlap_graph)
    }

    # Compute metrics
    nodes <- igraph::V(g)$name
    if (is.null(nodes)) nodes <- seq_len(igraph::vcount(g))
    gt <- groundtruth[groundtruth[,1] %in% nodes & groundtruth[,2] %in% nodes, , drop = FALSE]

    if (dir_g) {
      pairs <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
      pairs <- pairs[pairs$from != pairs$to, ]
    } else {
      pairs <- as.data.frame(t(combn(nodes, 2)), stringsAsFactors = FALSE)
      colnames(pairs) <- c("from", "to")
    }
    key <- paste(pairs[,1], pairs[,2], sep = "|")
    true <- key %in% paste(gt[,1], gt[,2], sep = "|")
    pred <- key %in% paste(igraph::as_data_frame(g)[,1], igraph::as_data_frame(g)[,2], sep = "|")

    TP <- sum(true & pred)
    TN <- sum(!true & !pred)
    FP <- sum(!true & pred)
    FN <- sum(true & !pred)

    acc  <- (TP + TN) / length(true)
    prec <- if (TP + FP > 0) TP / (TP + FP) else NA
    rec  <- if (TP + FN > 0) TP / (TP + FN) else NA
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA

    pi <- sum(true) / length(true)
    auprc <- if (!is.na(prec) && !is.na(rec) && prec >= pi) {
      prec * rec + (1 - rec) * (prec + pi) / 2
    } else NA

    n_g <- igraph::ecount(g)
    n_ov <- igraph::ecount(overlap_graph)
    val_pct <- if (n_g > 0) n_ov / n_g else NA

    list(edges = overlap_graph,
         metrics = data.frame(Accuracy = acc,
                              F1 = f1,
                              AUPRC = auprc,
                              Validated_Percentage = val_pct,
                              stringsAsFactors = FALSE))
  }

  if (!calculate_metrics) {
    # Filter out NULLs? Keep them as NULL in list
    names(res) <- names(network_graph)
    return(res)
  }

  # Extract edges and metrics
  edges <- lapply(res, `[[`, "edges")
  names(edges) <- names(network_graph)

  metrics <- do.call(rbind, lapply(res, function(x) {
    if (is.null(x$metrics)) {
      # Should not happen for calculate_metrics=TRUE, but just in case
      data.frame(Accuracy = NA, F1 = NA, AUPRC = NA, Validated_Percentage = NA)
    } else {
      x$metrics
    }
  }))
  metrics$Sample <- names(edges)
  rownames(metrics) <- NULL

  list(edges = edges, metrics = metrics)
}

#' Identifying sample-specific hubs with multiple centrality measures
#'
#' @title Network_hub
#' @param network_graph A list of igraph objects, each representing a sample-specific network.
#' @param directed Logical, whether the network is directed (TRUE) or undirected (FALSE).
#' @param top_percentage Numeric, top percentage of nodes to consider as hubs (default 0.2).
#' @param centrality_methods Character vector of centrality measures to use.
#'   Supported: `"degree"`, `"strength"`, `"betweenness"`, `"closeness"`,
#'   `"eigenvector"`, `"pagerank"`, `"hub_score"`, `"authority_score"`,
#'   `"subgraph"`, `"harmonic"`, `"kcore"`. Default `"degree"` ensures backward compatibility.
#' @param combine_method Character, how to combine hubs from multiple centrality methods.
#'   Options:
#'   \describe{
#'     \item{`"intersection"`}{Node must be a hub in **all** selected centrality methods.}
#'     \item{`"union"`}{Node must be a hub in **at least one** centrality method. }
#'     \item{`"majority"`}{Node must be a hub in **more than half** of the selected methods. }
#'     \item{`"rank_mean"`}{Averages the percentile rank of each node across all centralities, then selects the top `top_percentage` nodes. }
#'   }
#'   Ignored when only one centrality method is used (default for multiple methods: `"intersection"`).
#' @param mode Character, mode for directed graphs in degree, strength, closeness, harmonic, kcore.
#'   If `NULL` (default), `"out"` is used for directed, `"all"` for undirected networks.
#' @import igraph
#' @return A list object: a list of hubs for each sample.
#'
#' @references
#' Barabási A-L, Oltvai ZN. Network biology: understanding the cell’s functional organization. Nature Reviews. Genetics. 2004;5(2):101–113.
#' Wang M, Wang H, Zheng H. A mini review of node centrality metrics in biological networks. International Journal of Network Dynamics and Intelligence. 2022;1(1):99–110.
#'
#' @examples
#' \donttest{
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5],
#'                             bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' # Classic degree‑based hubs
#' hubs_degree <- Network_hub(bulk_LIONESS_res)
#' # Robust hubs (intersection of degree, betweenness, eigenvector)
#' hubs_robust <- Network_hub(bulk_LIONESS_res,
#'                            centrality_methods = c("degree", "betweenness", "eigenvector"))
#' }
#' @export
Network_hub <- function(network_graph, directed = TRUE, top_percentage = 0.2,
                        centrality_methods = "degree",
                        combine_method = c("intersection", "union", "majority", "rank_mean"),
                        mode = NULL) {
  combine_method <- match.arg(combine_method)
  if (is.null(mode)) mode <- if (directed) "out" else "all"
  centrality_methods <- tolower(centrality_methods)

  # Helper to compute a single centrality
  get_centrality <- function(g, method) {
    switch(method,
      degree          = igraph::degree(g, mode = mode),
      strength        = igraph::strength(g, mode = mode),
      betweenness     = igraph::betweenness(g, directed = directed),
      closeness       = igraph::closeness(g, mode = mode),
      eigenvector     = {
        if (directed) {
          warning("eigenvector centrality not suitable for directed graphs; using pagerank instead")
          igraph::page_rank(g, directed = TRUE)$vector
        } else igraph::eigen_centrality(g, directed = TRUE)$vector
      },
      pagerank        = igraph::page_rank(g, directed = directed)$vector,
      hub_score       = igraph::hub_score(g)$vector,
      authority_score = igraph::authority_score(g)$vector,
      subgraph        = igraph::subgraph_centrality(g),
      harmonic        = igraph::harmonic_centrality(g, mode = mode),
      kcore           = igraph::coreness(g, mode = mode)
    )
  }

  network_hubs <- lapply(seq_along(network_graph), function(i) {
    g <- network_graph[[i]]
    nms <- igraph::V(g)$name
    if (is.null(nms)) nms <- as.character(seq_len(igraph::vcount(g)))

    # Per‑method hubs and percentile ranks
    hub_lists <- list()
    rank_list <- list()
    for (met in centrality_methods) {
      val <- tryCatch(get_centrality(g, met), error = function(e) rep(0, length(nms)))
      names(val) <- nms
      val[is.na(val) | is.nan(val)] <- 0
      non0 <- val != 0
      if (sum(non0) > 0) {
        k <- ceiling(top_percentage * sum(non0))
        hub_lists[[met]] <- names(sort(val[non0], decreasing = TRUE))[1:min(k, sum(non0))]
      } else {
        hub_lists[[met]] <- character(0)
      }
      rank_list[[met]] <- rank(val, ties.method = "average") / length(val)
    }

    # Only one method: return its hubs directly
    if (length(centrality_methods) == 1) return(hub_lists[[1]])

    # Multiple methods: combine according to chosen strategy
    if (combine_method == "rank_mean") {
      avg_rank <- rowMeans(do.call(cbind, rank_list), na.rm = TRUE)
      k <- ceiling(top_percentage * length(nms))
      return(names(sort(avg_rank, decreasing = TRUE))[1:min(k, length(nms))])
    }

    hub_mat <- sapply(hub_lists, function(h) nms %in% h)
    selected <- switch(combine_method,
      intersection = rowSums(hub_mat) == ncol(hub_mat),
      union        = rowSums(hub_mat) > 0,
      majority     = rowMeans(hub_mat) > 0.5
    )
    nms[selected]
  })

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
#' @references
#' Csárdi G, Nepusz T (2006). “The igraph software package for complex network research.” InterJournal, Complex Systems, 1695.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' Network_module_res <- Network_module(bulk_LIONESS_res)
#'
#' @export
Network_module <- function(network_graph, method = "MCL", directed = TRUE, modulesize = 3) {

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
#' @references
#' Simpson EH. Measurement of diversity. Nature. 1949;163(4148):688–688.
#' Jaccard P. The distribution of the flora in the alpine zone. The New Phytologist. 1912;11(2):37–50.
#' Lin D. An information-theoretic definition of similarity. In: Proceedings of the Fifteenth International Conference on Machine Learning. San Francisco, CA, USA: Morgan Kaufmann Publishers Inc.; 1998. p. 296–304. (ICML ’98).
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


#' Performing clustering analysis of samples using the calculated sample-sample similarity matrix
#'
#' @param sim_matrix A square similarity matrix (values between 0 and 1, 1 = identical).
#' @param method Clustering method to use. Options are:
#'   \describe{
#'     \item{`"hclust"`}{Hierarchical clustering (complete linkage). Returns an `hclust` object and does **not** require the `k` parameter.}
#'     \item{`"pam"`}{Partitioning Around Medoids. Robust to outliers and works directly on the distance matrix. Requires `k`.}
#'     \item{`"spectral"`}{Spectral clustering using a normalized Laplacian. Can capture non‑convex cluster shapes. Requires `k`.}
#'     \item{`"kmeans"`}{Standard K‑means applied after multidimensional scaling of the distance matrix. Requires `k`.}
#'   }
#' @param k Number of clusters (required for all methods except `"hclust"`).
#'
#' @import stats
#' @import cluster
#'
#' @return For `"hclust"`, an `hclust` object. For other methods, a list with
#'   `cluster` (assignment vector) and the model object.
#'
#' @references
#' Murtagh F, Contreras P. Algorithms for hierarchical clustering: an overview. WIREs Data Mining and Knowledge Discovery. 2012;2(1):86–97.
#' Kaufman L, Rousseeuw PJ. Finding groups in data: an introduction to cluster analysis. John Wiley & Sons; 2009.
#' Ng A, Jordan M, Weiss Y. On spectral clustering: analysis and an algorithm. In: Dietterich T, Becker S, Ghahramani Z, editors. Advances in Neural Information Processing Systems. Vol. 14. MIT Press; 2001.
#' Steinley D. K-means clustering: a half-century synthesis. The British Journal of Mathematical and Statistical Psychology. 2006;59(Pt 1):1–34.
#'
#' @examples
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = 'pearson')
#' network_sim <- Sim.network(bulk_LIONESS_res,  bulk_LIONESS_res)
#' res <- sample_cluster(network_sim)
#' plot(res)
#'
#' @export
sample_cluster <- function(sim_matrix,
                           method = "hclust",
                           k = NULL) {
  method <- match.arg(method)
  dist_mat <- stats::as.dist(1 - sim_matrix)

  if (method == "hclust") {
    return(stats::hclust(dist_mat, method = "complete"))
  }

  if (is.null(k)) stop("'k' must be specified for method '", method, "'")

  if (method == "pam") {
    if (!requireNamespace("cluster", quietly = TRUE))
      stop("Package 'cluster' is required for PAM. Please install it.")
    pam_fit <- cluster::pam(dist_mat, k = k, diss = TRUE)
    return(list(cluster = pam_fit$clustering, pam_obj = pam_fit))
  }

  if (method == "kmeans") {
    mat <- as.matrix(dist_mat)
    mds <- stats::cmdscale(mat, k = min(nrow(mat) - 1, 10))
    km <- stats::kmeans(mds, centers = k)
    return(list(cluster = km$cluster, kmeans_obj = km))
  }

  if (method == "spectral") {
    S <- as.matrix(sim_matrix)
    diag(S) <- 0
    D <- diag(rowSums(S))
    L <- D - S
    D_inv_sqrt <- diag(1 / sqrt(diag(D)))
    L_sym <- D_inv_sqrt %*% L %*% D_inv_sqrt
    eig <- eigen(L_sym, symmetric = TRUE)
    vectors <- eig$vectors[, (ncol(eig$vectors) - k + 1):ncol(eig$vectors), drop = FALSE]
    vectors <- vectors / sqrt(rowSums(vectors^2))
    set.seed(123)
    km <- stats::kmeans(vectors, centers = k)
    return(list(cluster = km$cluster, spectral_obj = list(eigenvalues = eig$values, vectors = vectors)))
  }
}
