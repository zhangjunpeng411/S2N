#' Sample-specific network inferred with the SSN (single-sample network) method
#'
#' @title SSN
#' @param reg_normal A data frame object, the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A data frame object, the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A data frame object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A data frame object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network
#' @import utils
#' @import stats
#' @return A list object containing regulator-target interactions for each sample.
#' @examples
#'
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_SSN_res <- SSN(bulk_ncR_normal[1:100, ], bulk_tar_normal[1:500, ], bulk_ncR_cancer[1:100, ], bulk_tar_cancer[1:500, ], cormethod = "pearson")
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_SSN_res <- SSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5], scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = "pearson")
#'
#' @references
#' Liu X, Wang Y, Ji H, Aihara K, Chen L. Personalized characterization of diseases using sample-specific networks. Nucleic Acids Res. 2016;44(22):e164.
#'
#' @export
SSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer, cormethod = "pearson", p.value.cutoff = 0.05) {

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

    list(
      regulator = reg,
      target = tar,
      r1 = r1,
      x_norm = x,
      y_norm = y
    )

  })

  samples_num <- ncol(reg_cancer)
  nn <- ncol(reg_normal)

  ssn_net <- lapply(seq_len(samples_num), function(k) {
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

  ssn_net <- lapply(seq(ssn_net), function(i) ssn_net[[i]][which(as.numeric(ssn_net[[i]][, 7]) < p.value.cutoff), seq(2)])

  return(ssn_net)
}

#' Sample-specific network inferred with the Paired-SSN (paired single sample network) method
#'
#' @title PairedSSN
#' @param reg_normal A data frame object, the input regulator expression data in normal samples, columns are samples and rows are regulators.
#' @param tar_normal A data frame object, the input target gene expression data in normal samples, columns are samples and rows are target genes.
#' @param reg_cancer A data frame object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A data frame object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be used.
#' @param p.value.cutoff Significance p-value for identifying sample-specific regulatory network
#' @import utils
#' @import stats
#' @return A list object containing regulator-target interactions for each sample.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_PaiedSSN_res <- PairedSSN(bulk_ncR_normal[1:100, ], bulk_tar_normal[1:500, ], bulk_ncR_cancer[1:100, ], bulk_tar_cancer[1:500, ], cormethod = "pearson")
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_PaiedSSN_res <- PairedSSN(scRNA_ncR_normal[1:100, 1:5], scRNA_tar_normal[1:500, 1:5], scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = "pearson")
#'
#' @references
#' Guo WF, Zhang SW, Zeng T, Li Y, Gao J, Chen L. A novel network control model for identifying personalized driver genes in cancer. PLoS Comput Biol. 2019;15(11):e1007520.
#'
#' @export
PairedSSN <- function(reg_normal, tar_normal, reg_cancer, tar_cancer, cormethod = "pearson", p.value.cutoff = 0.05) {

  reg_intersect <- intersect(rownames(reg_normal), rownames(reg_cancer))
  reg_intersect <- reg_intersect[!is.na(reg_intersect)]
  tar_intersect <- intersect(rownames(tar_normal), rownames(tar_cancer))
  tar_intersect <- tar_intersect[!is.na(tar_intersect)]

  reg_normal <- reg_normal[reg_intersect, , drop = FALSE]
  tar_normal <- tar_normal[tar_intersect, , drop = FALSE]
  reg_cancer <- reg_cancer[reg_intersect, , drop = FALSE]
  tar_cancer <- tar_cancer[tar_intersect, , drop = FALSE]

  all_pairs <- expand.grid(
    regulator = reg_intersect,
    target = tar_intersect,
    stringsAsFactors = FALSE
  )

  pairedssn_net <- lapply(seq_len(ncol(reg_cancer)), function(k) {
    cat("Processing sample", k, "/", ncol(reg_cancer), "\n")

    curr_reg_cancer <- reg_cancer[, k, drop = FALSE]
    curr_tar_cancer <- tar_cancer[, k, drop = FALSE]
    curr_reg_normal <- reg_normal[, k, drop = FALSE]
    curr_tar_normal <- tar_normal[, k, drop = FALSE]

    ssn_cancer <- SSN_paired(reg_normal, tar_normal, curr_reg_cancer, curr_tar_cancer, cormethod, nn = ncol(reg_normal))
    ssn_normal <- SSN_paired(reg_normal, tar_normal, curr_reg_normal, curr_tar_normal, cormethod, nn = ncol(reg_normal))

    process_p <- function(pmat) {
      pmat[is.na(pmat)] <- 1
      pmat[pmat >= p.value.cutoff] <- 0
      pmat[pmat != 0] <- 1
      pmat
    }

    diff_matrix <- abs(process_p(ssn_cancer$p) - process_p(ssn_normal$p))

    sig_indices <- which(diff_matrix > 0, arr.ind = TRUE)
    data.frame(
      regulator = rownames(diff_matrix)[sig_indices[, 1]],
      target = colnames(diff_matrix)[sig_indices[, 2]],
      stringsAsFactors = FALSE
    )
  })

  return(pairedssn_net)
}

#' Sample-specific network inferred with the LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples) method
#'
#' @title LIONESS
#' @param reg_cancer A data frame object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A data frame object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be used.
#' @import stats
#' @import utils
#' @return A list object containing regulator-target interactions for each sample.
#' @references
#' Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating Sample-Specific Regulatory Networks. iScience. 2019;14:226-240.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_LIONESS_res <- LIONESS(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = "pearson")
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_LIONESS_res <- LIONESS(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = "pearson")
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_LIONESS_res <- LIONESS(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], cormethod = "pearson")
#'
#' @export
LIONESS <- function(reg_cancer, tar_cancer, cormethod = "pearson") {

  N <- ncol(reg_cancer)
  PCC <- cor(x = t(reg_cancer), y = t(tar_cancer), method = cormethod)
  lioness_net <- list()

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
    lioness_net[[k]] <- data.frame("regulator" = rownames(x)[index[, 1]], "target" = colnames(x)[index[, 2]])
  }

  return(lioness_net)
}

#' Sample-specific network inferred with SCS (single-sample controller strategy) method
#'
#' @title SCS
#' @param mut_reg A matrix object, containing a list of mutation regulators.
#' @param reg_cancer A data frame object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A data frame object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
#' @param cormethod A character string indicating which type of correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman" can be used.
#' @param number  The number of RWR with random matrix.
#' @import stats
#' @import igraph
#' @import dplyr
#' @import Matrix
#' @import methods
#' @import utils
#' @return A list object containing regulator-target interactions for each sample.
#' @references
#' Guo WF, Zhang SW, Liu LL, Liu F, Shi QQ, Zhang L, Tang Y, Zeng T, Chen L. Discovering personalized driver mutation profiles of single samples in cancer by network control strategy. Bioinformatics. 2018;34(11):1893-1903.
#'
#' @examples
#' #bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' data(mut_ncR)
#' bulk_scs_res <- SCS(mut_ncR, bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], cormethod = "pearson", number = 10)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' data(mut_ncR)
#' scRNA_scs_res <- SCS(mut_ncR, scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], cormethod = "pearson", number = 10)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' data(mut_ncR)
#' ST_scs_res <- SCS(mut_ncR, ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], cormethod = "pearson", number = 10)
#'
#' @export
#'
SCS <- function(mut_reg, reg_cancer, tar_cancer, cormethod = "pearson", number = 50) {

  edge0 <- expand.grid(regulator = mut_reg, target = rownames(tar_cancer))
  colnames(edge0) <- c("regulator", "target")
  node0 <- unique(c(edge0$regulator, edge0$target))

  SNP_name <- mut_reg
  reg_EXPR_data <- reg_cancer
  reg_EXPR_name <- rownames(reg_cancer)
  EXPR_data <- tar_cancer
  EXPR_name <- rownames(tar_cancer)
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

  index_delete_sample <- which(sapply(Targets_genes_sampels, length) *
                                sapply(Mutation_genes_sampels, length) == 0)
  if (length(index_delete_sample) > 0) {
    Targets_genes_sampels <- Targets_genes_sampels[-index_delete_sample]
    Mutation_genes_sampels <- Mutation_genes_sampels[-index_delete_sample]
  }

  Targets_C_genes <- lapply(Targets_genes_sampels, function(x) intersect(x, node0))
  Constrained_B_genes <- lapply(Mutation_genes_sampels, function(x) intersect(x, node0))

  SCS_net <- list()
  for (i in seq_along(Constrained_B_genes)) {
    cat("Processing sample", i, "/", length(Constrained_B_genes), "\n")

    B_fda_genes <- Constrained_B_genes[[i]]
    C_genes <- Targets_C_genes[[i]]

    edge0_B_C <- edge0[
      edge0$regulator %in% B_fda_genes &
      edge0$target %in% C_genes,
    ]
    if (nrow(edge0_B_C) == 0) next

    edge0_network <- igraph::graph_from_data_frame(edge0_B_C, directed = TRUE)
    edge0_MultiplexObject <- create.multiplex(list(edge0_network), directed = TRUE)
    AdjMatrix_edge0 <- compute.adjacency.matrix(edge0_MultiplexObject)

    RWR_results_edge0 <- RWR_MN(AdjMatrix_edge0, edge0_MultiplexObject, B_fda_genes, scores = 0.0001)
    colnames(RWR_results_edge0)[3] <- "gdtruth_score"

    RWR_results_dir_srand <- lapply(1:number, function(j) {
      dir_srand <- dir_generate_srand(as.matrix(AdjMatrix_edge0))
      dir_srand_AdjMatrix <- as(as.matrix(dir_srand), "dgCMatrix")
      RWR_MN(dir_srand_AdjMatrix, edge0_MultiplexObject, B_fda_genes, scores = 0.0001)
    })

    RWR_results_dir_srand_merge <- Reduce(
      function(x, y) merge(x, y, by = c("SeedGenes", "NodeNames")),
      RWR_results_dir_srand
    )
    RWR_results_edge0_dir_srand <- merge(RWR_results_edge0, RWR_results_dir_srand_merge)

    p_value <- sapply(1:nrow(RWR_results_edge0_dir_srand), function(k) {
      rand_scores <- as.numeric(RWR_results_edge0_dir_srand[k, -(1:3)])
      if (all(rand_scores <= RWR_results_edge0_dir_srand$gdtruth_score[k])) {
        return(RWR_results_edge0_dir_srand$gdtruth_score[k])
      } else {
        return(0)
      }
    })

    SCS_net[[i]] <- RWR_results_edge0_dir_srand[p_value > 0, 1:3]
  }

  return(SCS_net)
}


#' Sample-specific network inferred with the CSN (cell-specific network) method, adapted from the CSmiR method
#'
#' @title CSN
#' @param reg_cancer A data frame object, the input regulator expression data in tumor samples, columns are samples and rows are regulators.
#' @param tar_cancer A data frame object, the input target gene expression data in tumor samples, columns are samples and rows are target genes.
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
#' @importFrom WGCNA pmedian
#' @return A list object containing regulator-target interactions for each sample.
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
CSN <- function(reg_cancer, tar_cancer, boxsize = 0.1, bootstrap = FALSE, bootstrap_betw_point = 5, bootstrap_num = 100, p.value.cutoff = 0.05, num.cores = 2) {

  reg_cancer <- data.frame(reg_cancer)
  tar_cancer <- data.frame(tar_cancer)

  reg_data_cancer <- t(reg_cancer)
  tar_data_cancer <- t(tar_cancer)
  rm("reg_cancer", "tar_cancer")
  gc()

  if(bootstrap){
    CSN_net <- CSN_net_bootstrap(reg_data_cancer,
                                 tar_data_cancer,
                                 boxsize = boxsize,
                                 bootstrap_betw_point = bootstrap_betw_point,
                                 bootstrap_num = bootstrap_num,
                                 p.value.cutoff = p.value.cutoff
    )
  } else {
    CSN_net <- CSN_net(reg_data_cancer,
                       tar_data_cancer,
                       boxsize = boxsize,
                       p.value.cutoff = p.value.cutoff,
                       num.cores = num.cores
    )
  }

  return(CSN_net)
}

#' Identifying sample-specific miRNA regulation using Scan (Sample-speCific miRNA regulAtioN) and statistical perturbation strategy
#'
#' @title Scan.perturb
#' @param reg_cancer miRNA expression data with columns are samples and rows are regulators.
#' @param tar_cancer mRNA expression data with columns are samples and rows are mRNAs.
#' @param method Methods for calculating correlations, one of "Pearson", "Spearman", "Kendall".
#' @param p.value.cutoff Significance level for the identified regulator-target interactions.
#' @param num.cores Number of CPU cores when parallel computation.
#' @importFrom igraph as_data_frame
#' @importFrom igraph graph_from_biadjacency_matrix
#' @import parallel
#' @import foreach
#' @import doParallel
#' @importFrom WGCNA cor
#' @importFrom WGCNA corAndPvalue
#'
#' @return A list object containing RNA-RNA interactions for each sample.
#' @references
#' Zhang J, Liu L, Wei X, Zhao C, Luo Y, Li J, Le TD. Scanning sample-specific miRNA regulation from bulk and single-cell RNA-sequencing data. BMC Biol. 2024;22(1):218.
#'
#' @examples
#' # bulk data
#' data(bulk_ncR_tar_normal_cancer)
#' bulk_Scan.perturb_res <- Scan.perturb(bulk_ncR_cancer[1:100, 1:5], bulk_tar_cancer[1:500, 1:5], method = "Pearson", p.value.cutoff = 0.05)
#'
#' # scRNA data
#' data(scRNA_ncR_tar_normal_cancer)
#' scRNA_Scan.perturb_res <- Scan.perturb(scRNA_ncR_cancer[1:100, 1:5], scRNA_tar_cancer[1:500, 1:5], method = "Pearson", p.value.cutoff = 0.05)
#'
#' # ST data
#' data(ST_ncR_tar_cancer)
#' ST_Scan.perturb_res <- Scan.perturb(ST_ncR_cancer[1:100, 1:5], ST_tar_cancer[1:500, 1:5], method = "Pearson", p.value.cutoff = 0.05)
#'
#' @export
Scan.perturb <- function(reg_cancer, tar_cancer,
                         method = c("Pearson", "Spearman", "Kendall"),
                         p.value.cutoff = 0.05,
                         num.cores = 2){

  if (method == "Pearson"){
    res.all <- Pearson(reg_cancer, tar_cancer, p.value.cutoff)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph", "WGCNA"), .export = c("Pearson")) %dopar% {
      Pearson(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

  } else if (method == "Spearman"){
    res.all <- Spearman(reg_cancer, tar_cancer, p.value.cutoff)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph", "WGCNA"), .export = c("Spearman")) %dopar% {
      Spearman(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

  } else if (method == "Kendall"){
    res.all <- Kendall(reg_cancer, tar_cancer, p.value.cutoff)

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)

    res.single <- foreach(i = seq(ncol(reg_cancer)), .packages = c("igraph", "WGCNA"), .export = c("Kendall")) %dopar% {
      Kendall(reg_cancer[, -i], tar_cancer[, -i], p.value.cutoff)
    }

    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()
  }

  Scan.perturb.incident <- lapply(seq(res.single), function(i)
    abs(res.single[[i]] - res.all))

  Scan.perturb.graph <- lapply(seq(res.single), function(i)
    graph_from_biadjacency_matrix(t(Scan.perturb.incident[[i]])))

  Scan.perturb_net <- lapply(seq(res.single), function(i)
    as_data_frame(Scan.perturb.graph[[i]]))

  return(Scan.perturb_net)
}


#' Validating sample-specific network
#'
#' @title Network_validation
#' @param network List object: a list of sample-specific networks.
#' @param groundtruth A data frame object, experimentally validated RNA-RNA interactions
#' @import igraph
#' @import foreach
#' @return List object: a list of validated RNA-RNA interactions for each sample.
#'
#' @export
Network_validation <- function(network, groundtruth) {

  network_graph <- lapply(seq(network), function(i) igraph::make_graph(c(t(network[[i]][, 1:2])), directed = TRUE))
  groundtruth_graph <- make_graph(c(t(groundtruth[, c(1, 2)])), directed = T)
  network_validated <- foreach::foreach(i = seq(network),.packages = c("igraph", 'foreach', 'doParallel')) %dopar% {as_data_frame(network_graph[[i]] %s% groundtruth_graph)} # transform to a data frame

  return(network_validated)
}

#' Identifying sample-specific hub regulators
#'
#' @title Network_hub
#' @param network List object: a list of of sample-specific networks.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param top_percentage A number indicating the top percenatge of hub regulators, default 0.2.
#' @import igraph
#' @return List object: a list of hub regulators for each sample.
#' @export
#'
#' @examples
#' data(network)
#' hubs <- Network_hub(network)
#' @references
#'Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD. Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data. BMC Bioinformatics.2021; 22(1):578.
Network_hub <- function(network, directed = TRUE, top_percentage = 0.2) {

  network_graph <- lapply(seq(network), function(i) igraph::make_graph(c(t(na.omit(network[[i]][, 1:2]))), directed = directed))
  if(directed == TRUE){
    network_degree <- lapply(seq(network), function(i) igraph::degree(network_graph[[i]], mode = "out"))
  } else {
    network_degree <- lapply(seq(network), function(i) igraph::degree(network_graph[[i]]))
  }
  network_hubs <- lapply(seq(network), function(i) names(sort(network_degree[[i]][which(network_degree[[i]] != 0)], decreasing = TRUE))[1:ceiling(top_percentage * length(which(network_degree[[i]] != 0)))])

  return(network_hubs)
}

#' Predicting sample-specific modules
#'
#' @title Network_module
#' @param network List object, a list of sample-specific networks.
#' @param method Possible methods include FN, MCL, LINKCOMM, MCODE, betweenness, infomap, prop, eigen, louvain, walktrap.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param modulesize The size of modules, e.g. 3.
#' @import igraph
#' @importFrom MCL mcl
#' @importFrom linkcomm getLinkCommunities
#' @return List object: a list of sample-specific modules.
#' @export
#'
#' @examples
#' data(network)
#' Network_module_res <- Network_module(network)
Network_module <- function(network, method = "MCL", directed = TRUE, modulesize = 3) {

  network_Cluster <- list()
  edgelist <- list()
  network_Cluster_result <- list()
  size <- list()

  for (k in seq(network)) {

    if (method == "FN" | method == "MCL" | method == "MCODE") {
      network_Cluster[[k]] <- cluster(graph_from_data_frame(network[[k]], directed = directed),
                                     method = method, directed = directed)
    } else if (method == "LINKCOMM") {
      edgelist[[k]] <- get.edgelist(graph_from_data_frame(network[[k]], directed = directed))
      network_Cluster[[k]] <- getLinkCommunities(edgelist[[k]], directed = directed)$nodeclusters
    }

    if (method == "FN" | method == "MCL") {
      network_Cluster_result[[k]] <- lapply(seq_len(max(network_Cluster[[k]])),
                                           function(i) rownames(as.matrix(network_Cluster[[k]]))[which(network_Cluster[[k]] == i)])
      size[[k]] <- unlist(lapply(seq_len(max(network_Cluster[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "LINKCOMM") {
      network_Cluster_result[[k]] <- lapply(seq_len(max(c(network_Cluster[[k]]$cluster))),
                                           function(i) as.character(network_Cluster[[k]]$node[which(c(network_Cluster[[k]]$cluster) == i)]))
      size[[k]] <- unlist(lapply(seq_len(max(c(network_Cluster[[k]]$cluster))), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "MCODE") {
      network_Cluster[[k]] <- network_Cluster[[k]] + 1
      network_Cluster_result[[k]] <- lapply(seq_len(max(network_Cluster[[k]])),
                                           function(i) rownames(as.matrix(network_Cluster[[k]]))[which(network_Cluster[[k]] == i)])
      size[[k]] <- unlist(lapply(seq_len(max(network_Cluster[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "betweenness") {
      network_Cluster_result[[k]] <- cluster_edge_betweenness(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "infomap") {
      network_Cluster_result[[k]] <- cluster_infomap(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "prop") {
      network_Cluster_result[[k]] <- cluster_label_prop(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "eigen") {
      network_Cluster_result[[k]] <- cluster_leading_eigen(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "louvain") {
      network_Cluster_result[[k]] <- cluster_louvain(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    } else if (method == "walktrap") {
      network_Cluster_result[[k]] <- cluster_walktrap(graph_from_data_frame(network[[k]], directed = directed))
      size[[k]] <- unlist(lapply(seq_len(length(network_Cluster_result[[k]])), function(i) length(network_Cluster_result[[k]][[i]])))
      network_Cluster_result[[k]] <- lapply(which(size[[k]] >= modulesize), function(i) network_Cluster_result[[k]][[i]])
    }
  }

  return(network_Cluster_result)

}

#' Calculating sample-sample similarity matrix between two lists of sample-specific networks
#'
#' @title Sim.network
#' @param net1 The first list of sample-specific networks.
#' @param net2 The second list of sample-specific networks.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param method Methods for calculating similatiry between two two lists of sample-specific networks, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @import igraph
#' @return Matrix object, a sample-sample similarity matrix in terms of sample-specific networks.
#' @export
#'
#' @examples
#' data(network)
#' res <- Sim.network(network, network)
Sim.network <- function(net1, net2, directed = TRUE, method = "Simpson") {

  if (class(net1) != "list" | class(net2) != "list") {
    stop("Please check your input network! The input network should be list object! \n")
  }

  m <- length(net1)
  n <- length(net2)
  Sim <- matrix(NA, m, n)

  if (method == "Simpson") {
    for (i in seq(m)) {
      for (j in seq(n)) {
        net1_graph_interin <- igraph::make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
        net2_graph_interin <- igraph::make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
        overlap_interin <- nrow(igraph::as_data_frame(net1_graph_interin %s% net2_graph_interin))
        Sim[i, j] <- overlap_interin / min(nrow(net1[[i]]), nrow(net2[[j]]))
      }
    }
  } else if (method == "Jaccard") {
    for (i in seq(m)) {
      for (j in seq(n)) {
        net1_graph_interin <- igraph::make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
        net2_graph_interin <- igraph::make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
        overlap_interin <- nrow(igraph::as_data_frame(net1_graph_interin %s% net2_graph_interin))
        Sim[i, j] <- overlap_interin / nrow(unique(rbind(net1[[i]], net2[[j]])))
      }
    }
  } else if (method == "Lin") {
    for (i in seq(m)) {
      for (j in seq(n)) {
        net1_graph_interin <- igraph::make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
        net2_graph_interin <- igraph::make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
        overlap_interin <- nrow(igraph::as_data_frame(net1_graph_interin %s% net2_graph_interin))
        Sim[i, j] <- 2 * overlap_interin / (nrow(net1[[i]]) + nrow(net2[[j]]))
      }
    }
  }

  return(Sim)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific networks
#'
#' @title Network_sim_hc
#' @param network List object, a list of sample-specific networks.
#' @param directed A logical value indicating a network directed (TRUE) or undirected (FALSE).
#' @param method Methods for calculating similatiry between two two lists of sample-specific networks, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @import igraph
#' @return Hierarchical cluster analysis result of samples.
#' @export
#'
#' @examples
#' data(network)
#' res <- Network_sim_hc(network)
Network_sim_hc <- function(network, directed = TRUE, method = "Simpson") {

  network_Sim <- Sim.network(network, network, directed = directed, method = method)
  hclust_res_network <- hclust(as.dist(1 - network_Sim), "complete")

  return(hclust_res_network)
}


#' Calculating sample-sample similarity matrix between two lists of sample-specific hubs
#'
#' @title Sim.hub
#' @param hub1 The first list of sample-specific hubs.
#' @param hub2 The second list of sample-specific hubs.
#' @param method Methods for calculating similatiry between two two lists of sample-specific hubs, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Matrix object, a sample-sample similarity matrix in terms of sample-specific hubs.
#' @export
#'
#' @examples
#' data(hub_ncRNAs)
#' res <- Sim.hub(hub_ncRNAs, hub_ncRNAs)
Sim.hub <- function(hub1, hub2, method = "Simpson"){

  if(class(hub1)!="list" | class(hub2)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
  }
  m <- length(hub1)
  n <- length(hub2)
  Sim <- matrix(NA, m, n)

  if (method == "Simpson") {
    for (i in seq(m)){
      for (j in seq(n)){
        overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
        Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
      }
    }
  } else if (method == "Jaccard") {
    for (i in seq(m)){
      for (j in seq(n)){
        overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
        Sim[i, j] <- overlap_interin/length(union(hub1[[i]], hub2[[j]]))
      }
    }
  } else if (method == "Lin") {
    for (i in seq(m)){
      for (j in seq(n)){
        overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
        Sim[i, j] <- 2 * overlap_interin/(length(hub1[[i]]) + length(hub2[[i]]))
      }
    }
  }

  return(Sim)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific hubs
#'
#' @title Hub_sim_hc
#' @param hub List object, a list of sample-specific hubs.
#' @param method Methods for calculating similatiry between two two lists of sample-specific hubs, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Hierarchical cluster analysis result of samples.
#' @export
#'
#' @examples
#' data(hub_ncRNAs)
#' res <- Hub_sim_hc(hub_ncRNAs)
Hub_sim_hc <- function(hub, method = "Simpson"){

  hub_Sim <- Sim.hub(hub, hub, method = method)
  hclust_res_hub <- hclust(as.dist(1 - hub_Sim), "complete")

  return(hclust_res_hub)
}

#' Calculating sample-sample similarity matrix between two lists of sample-specific modules
#'
#' @title Sim.module
#' @param Module1 List object, the first list of modules.
#' @param Module2 List object, the second list of modules.
#' @param method Methods for calculating similatiry between two modules, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Matrix object, a sample-sample similarity matrix in terms of sample-specific modules.
#' @export
#'
Sim.module <- function(Module1, Module2, method = "Simpson"){

    if(class(Module1)!="list" | class(Module2)!="list") {
    stop("Please check your input modules! The input modules should be list object! \n")
    }

    m <- length(Module1)
    n <- length(Module2)
    Sim <- matrix(NA, m, n)

    if (method == "Simpson") {
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
	    min_vertex <- min(length(Module1[[i]]), length(Module2[[j]]))
	    Sim[i, j] <- overlap_vertex/min_vertex
	}
    }
    } else if (method == "Jaccard") {
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
	    union_vertex <- length(union(Module1[[i]], Module2[[j]]))
	    Sim[i, j] <- overlap_vertex/union_vertex
	}
    }
    } else if (method == "Lin") {
    for (i in seq(m)){
        for (j in seq(n)){
	    overlap_vertex <- length(intersect(Module1[[i]], Module2[[j]]))
            sum_vertex <- length(Module1[[i]]) + length(Module2[[j]])
            Sim[i, j] <- 2 * overlap_vertex/sum_vertex
	}
    }
    }

    if (m < n) {
	GS <- mean(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])))*m/n
    } else if (m == n) {
        GS <- mean(c(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])),
	           unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i]))))
    } else if (m > n) {
        GS <- mean(unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i])))*n/m
    }

    return(GS)
}

#' Performing hierarchical clustering analysis of samples using the calculated sample-sample similarity matrix in terms of sample-specific modules
#'
#' @title Module_sim_hc
#' @param module List object, a list of sample-specific module groups. Each module group contain a list of modules.
#' @param method Methods for calculating similatiry between two lists of sample-specific modules, select one of three mrthods (Simpson, Jaccard and Lin). Default method is Simpson.
#' @import stats
#' @return Hierarchical cluster analysis result of samples.
#' @export
#'
Module_sim_hc <- function(module, method = "Simpson"){

  n <- length(module)
  module_Sim <- matrix(NA, n, n)
  for (i in seq(n)) {
    for (j in seq(n)) {
        module_Sim[i, j] <- Sim.module(module[[i]], module[[j]], method = method)
    }
  }
  hclust_res_module <- hclust(as.dist(1 - module_Sim), "complete")

  return(hclust_res_module)
}
