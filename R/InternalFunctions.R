## Calulate z-score value for the SSN method Internal function
ssn_zscore <- function(deta, pcc, nn) {

    if (any(is.na(c(deta, pcc, nn))))
        return(NA)
    if (nn <= 1)
        stop("nn should be larger than 1")

    if (abs(pcc - 1) < 1e-10) {
        pcc <- 0.9999999999
    } else if (abs(pcc + 1) < 1e-10) {
        pcc <- -0.9999999999
    }

    denominator <- (1 - pcc^2)/(nn - 1)
    if (abs(denominator) < 1e-10) {
        warning("The denominator is close to zero, and the results are not stable.")
        return(NA)
    }

    z <- deta/denominator
    return(z)
}

## Calculate ssn score using paired normal and tumor data for the PairedSSN
## method Internal function
SSN_paired <- function(reg_normal_mat, tar_normal_mat, reg_sample, tar_sample, cormethod = "pearson",
    nn) {

    R0 <- cor(x = t(reg_normal_mat), y = t(tar_normal_mat), method = cormethod)

    inter1 <- intersect(rownames(reg_normal_mat), rownames(reg_sample))
    reg_normal_mat <- reg_normal_mat[inter1, ]
    reg_sample <- reg_sample[inter1, ]

    inter2 <- intersect(rownames(tar_normal_mat), rownames(tar_sample))
    tar_normal_mat <- tar_normal_mat[inter2, ]
    tar_sample <- tar_sample[inter2, ]

    R1 <- cor(x = t(cbind(reg_normal_mat, reg_sample)), y = t(cbind(tar_normal_mat,
        tar_sample)), method = cormethod)
    detal <- R1 - R0
    Z <- detal/((1 - R0^2)/(nn - 1))
    p <- 1 - pnorm(abs(Z))
    ssn_result <- list(detal = detal, p = p)

    return(ssn_result)
}

## Functions to perform Random Walk with Restart on Multiplex Networks
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
simplify.layers <- function(Input_Layer) {

    ## Undirected Graphs
    Layer <- as_undirected(Input_Layer, mode = c("collapse"), edge.attr.comb = igraph_opt("edge.attr.comb"))

    ## Unweighted or Weigthed Graphs
    if (is_weighted(Layer)) {
        b <- 1
        weigths_layer <- E(Layer)$weight
        if (min(weigths_layer) != max(weigths_layer)) {
            a <- min(weigths_layer)/max(weigths_layer)
            range01 <- (b - a) * (weigths_layer - min(weigths_layer))/(max(weigths_layer) -
                min(weigths_layer)) + a
            E(Layer)$weight <- range01
        } else {
            E(Layer)$weight <- rep(1, length(weigths_layer))
        }
    } else {
        E(Layer)$weight <- rep(1, ecount(Layer))
    }

    ## Simple Graphs
    Layer <- igraph::simplify(Layer, remove.multiple = TRUE, remove.loops = TRUE,
        edge.attr.comb = mean)

    return(Layer)
}

## Add missing nodes in some of the layers
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
add.missing.nodes <- function(Layers, Nr_Layers, NodeNames) {

    add_vertices(Layers, length(NodeNames[which(!NodeNames %in% V(Layers)$name)]),
        name = NodeNames[which(!NodeNames %in% V(Layers)$name)])
}

## Is this R object a Multiplex object?
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
isMultiplex <- function(x) {
    is(x, "Multiplex")
}

## Scores for the seeds of the multiplex network
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
get.seed.scoresMultiplex <- function(Seeds, Number_Layers, tau) {

    Nr_Seeds <- length(Seeds)

    Seeds_Seeds_Scores <- rep(tau/Nr_Seeds, Nr_Seeds)
    Seed_Seeds_Layer_Labeled <- paste0(rep(Seeds, Number_Layers), sep = "_", rep(seq(Number_Layers),
        length.out = Nr_Seeds * Number_Layers, each = Nr_Seeds))

    Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled, Score = Seeds_Seeds_Scores,
        stringsAsFactors = FALSE)

    return(Seeds_Score)
}

## Geometric mean: R computation
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
geometric.mean <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)) {
        FinalScore[i] <- prod(Scores[seq(from = i, to = N * L, by = N)])^(1/L)
    }

    return(FinalScore)
}

## Create Multiplex and Multiplex-Heterogeneous objects
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
create.multiplex <- function(LayersList, ...) {

    if (!class(LayersList) == "list") {
        stop("The input object should be a list of graphs.")
    }


    Number_of_Layers <- length(LayersList)
    SeqLayers <- seq(Number_of_Layers)
    Layers_Name <- names(LayersList)

    if (!all(sapply(SeqLayers, function(x) is_igraph(LayersList[[x]])))) {
        stop("Not igraph objects")
    }

    Layer_List <- lapply(SeqLayers, function(x) {
        if (is.null(V(LayersList[[x]])$name)) {
            LayersList[[x]] <- set_vertex_attr(LayersList[[x]], "name", value = seq(1,
                vcount(LayersList[[x]]), by = 1))
        } else {
            LayersList[[x]]
        }
    })

    ## We simplify the layers
    Layer_List <- lapply(SeqLayers, function(x) simplify.layers(Layer_List[[x]]))

    ## We set the names of the layers.

    if (is.null(Layers_Name)) {
        names(Layer_List) <- paste0("Layer_", SeqLayers)
    } else {
        names(Layer_List) <- Layers_Name
    }

    ## We get a pool of nodes (Nodes in any of the layers.)
    Pool_of_Nodes <- sort(unique(unlist(lapply(SeqLayers, function(x) V(Layer_List[[x]])$name))))

    Number_of_Nodes <- length(Pool_of_Nodes)

    Layer_List <- lapply(Layer_List, add.missing.nodes, Number_of_Layers, Pool_of_Nodes)

    # We set the attributes of the layer
    counter <- 0
    Layer_List <- lapply(Layer_List, function(x) {
        counter <<- counter + 1
        set_edge_attr(x, "type", E(x), value = names(Layer_List)[counter])
    })


    MultiplexObject <- c(Layer_List, list(Pool_of_Nodes = Pool_of_Nodes, Number_of_Nodes_Multiplex = Number_of_Nodes,
        Number_of_Layers = Number_of_Layers))

    class(MultiplexObject) <- "Multiplex"

    return(MultiplexObject)
}

## Compute the matrices and perform the Random Walks
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
#' @export
#' @method compute.adjacency.matrix
compute.adjacency.matrix <- function(x, delta = 0.5, ...) {
    if (!isMultiplex(x)) {
        stop("Not a Multiplex object")
    }
    if (delta > 1 || delta <= 0) {
        stop("Delta should be between 0 and 1")
    }

    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers

    ## We impose delta=0 in the monoplex case.
    if (L == 1) {
        delta = 0
    }

    Layers_Names <- names(x)[seq(L)]

    ## IDEM_MATRIX.
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)

    counter <- 0
    Layers_List <- lapply(x[Layers_Names], function(x) {

        counter <<- counter + 1
        if (is_weighted(x)) {
            Adjacency_Layer <- as_adjacency_matrix(x, sparse = TRUE, attr = "weight")
        } else {
            Adjacency_Layer <- as_adjacency_matrix(x, sparse = TRUE)
        }

        Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)), order(colnames(Adjacency_Layer))]
        colnames(Adjacency_Layer) <- paste0(colnames(Adjacency_Layer), "_", counter)
        rownames(Adjacency_Layer) <- paste0(rownames(Adjacency_Layer), "_", counter)
        Adjacency_Layer
    })

    MyColNames <- unlist(lapply(Layers_List, function(x) unlist(colnames(x))))
    MyRowNames <- unlist(lapply(Layers_List, function(x) unlist(rownames(x))))
    names(MyColNames) <- c()
    names(MyRowNames) <- c()
    SupraAdjacencyMatrix <- (1 - delta) * (bdiag(unlist(Layers_List)))
    colnames(SupraAdjacencyMatrix) <- MyColNames
    rownames(SupraAdjacencyMatrix) <- MyRowNames

    offdiag <- (delta/(L - 1)) * Idem_Matrix

    i <- seq_len(L)
    Position_ini_row <- 1 + (i - 1) * N
    Position_end_row <- N + (i - 1) * N
    j <- seq_len(L)
    Position_ini_col <- 1 + (j - 1) * N
    Position_end_col <- N + (j - 1) * N

    for (i in seq_len(L)) {
        for (j in seq_len(L)) {
            if (j != i) {
                SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]), (Position_ini_col[j]:Position_end_col[j])] <- offdiag
            }
        }
    }

    SupraAdjacencyMatrix <- methods::as(SupraAdjacencyMatrix, "dgCMatrix")
    return(SupraAdjacencyMatrix)
}

## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
regular.mean <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)) {
        FinalScore[i] <- mean(Scores[seq(from = i, to = N * L, by = N)])
    }

    return(FinalScore)
}

## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
sumValues <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)) {
        FinalScore[i] <- sum(Scores[seq(from = i, to = N * L, by = N)])
    }

    return(FinalScore)
}

## Computes column normalization of an adjacency matrix
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
normalize.multiplex.adjacency <- function(x) {
    if (!is(x, "dgCMatrix")) {
        stop("Not a dgCMatrix object of Matrix package")
    }

    Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1, sparseResult = FALSE)))

    return(Adj_Matrix_Norm)
}

## Performs Random Walk with Restart on a Multiplex Network
## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/) Internal
## function
Random.Walk.Restart.Multiplex <- function(x, MultiplexObject, Seeds, r = 0.7, tau,
    MeanType = "Geometric", DispResults = "TopScores", ...) {

    ### We control the different values.
    if (!is(x, "dgCMatrix")) {
        stop("Not a dgCMatrix object of Matrix package")
    }

    if (!isMultiplex(MultiplexObject)) {
        stop("Not a Multiplex object")
    }

    L <- MultiplexObject$Number_of_Layers
    N <- MultiplexObject$Number_of_Nodes

    Seeds <- as.character(Seeds)
    if (length(Seeds) < 1 | length(Seeds) >= N) {
        stop("The length of the vector containing the seed nodes is not
         correct")
    } else {
        if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)) {
            stop("Some of the seeds are not nodes of the network")

        }
    }

    if (r >= 1 || r <= 0) {
        stop("Restart partameter should be between 0 and 1")
    }

    if (missing(tau)) {
        tau <- rep(1, L)/L
    } else {
        tau <- as.numeric(tau)
        if (sum(tau)/L != 1) {
            stop("The sum of the components of tau divided by the number of
           layers should be 1")
        }
    }

    if (!(MeanType %in% c("Geometric", "Arithmetic", "Sum"))) {
        stop("The type mean should be Geometric, Arithmetic or Sum")
    }

    if (!(DispResults %in% c("TopScores", "Alphabetic"))) {
        stop("The way to display RWRM results should be TopScores or
         Alphabetic")
    }

    ## We define the threshold and the number maximum of iterations for the
    ## random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)

    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1

    ## We compute the scores for the different seeds.
    Seeds_Score <- get.seed.scoresMultiplex(Seeds, L, tau)

    ## We define the prox_vector(The vector we will move after the first RWR
    ## iteration. We start from The seed. We have to take in account that the
    ## walker with restart in some of the Seed nodes, depending on the score we
    ## gave in that file).
    prox_vector <- matrix(0, nrow = NetworkSize, ncol = 1)

    prox_vector[which(colnames(x) %in% Seeds_Score[, 1])] <- (Seeds_Score[, 2])

    prox_vector <- prox_vector/sum(prox_vector)
    restart_vector <- prox_vector

    while (residue >= Threeshold) {

        old_prox_vector <- prox_vector
        prox_vector <- (1 - r) * (x %*% prox_vector) + r * restart_vector
        residue <- sqrt(sum((prox_vector - old_prox_vector)^2))
        iter <- iter + 1
    }

    NodeNames <- character(length = N)
    Score = numeric(length = N)

    rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
    rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])

    if (MeanType == "Geometric") {
        rank_global$Score <- geometric.mean(as.vector(prox_vector[, 1]), L, N)
    } else {
        if (MeanType == "Arithmetic") {
            rank_global$Score <- regular.mean(as.vector(prox_vector[, 1]), L, N)
        } else {
            rank_global$Score <- sumValues(as.vector(prox_vector[, 1]), L, N)
        }
    }

    if (DispResults == "TopScores") {
        ## We sort the nodes according to their score.
        Global_results <- rank_global[with(rank_global, order(-Score, NodeNames)),
            ]

        ### We remove the seed nodes from the Ranking and we write the results.
        Global_results <- Global_results[which(!Global_results$NodeNames %in% Seeds),
            ]
    } else {
        Global_results <- rank_global
    }

    rownames(Global_results) <- c()

    RWRM_ranking <- list(RWRM_Results = Global_results, Seed_Nodes = Seeds)

    class(RWRM_ranking) <- "RWRM_Results"
    return(RWRM_ranking)
}

## Perform Random Walk with Restart on a multiplex network Internal function
RWR_MN <- function(AdjMatrix_edge0 = AdjMatrix_edge0, edge0_MultiplexObject = edge0_MultiplexObject,
    SNP_name = SNP_name, scores = 1e-04) {

    AdjMatrixNorm_edge0 <- normalize.multiplex.adjacency(AdjMatrix_edge0)
    SeedGenes <- edge0_MultiplexObject[["Pool_of_Nodes"]][na.omit(match(SNP_name,
        edge0_MultiplexObject[["Pool_of_Nodes"]]))]
    RWR_edge0_Results_top <- list()
    for (i in 1:length(SeedGenes)) {
        RWR_edge0_Results_tmp <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_edge0,
            edge0_MultiplexObject, SeedGenes[i], r = 0.6)
        RWR_edge0_Results_tmp1 <- data.frame(SeedGenes = RWR_edge0_Results_tmp[["Seed_Nodes"]],
            RWR_edge0_Results_tmp[[1]])
        RWR_edge0_Results_top[[i]] <- RWR_edge0_Results_tmp1[RWR_edge0_Results_tmp1[,
            3] > scores, ]
    }
    RWR_edge0_Results_top_all <- Reduce(rbind, RWR_edge0_Results_top)

    return(RWR_edge0_Results_top_all)
}


## Construct random network with the same degree as the original one Internal
## function
dir_generate_srand <- function(s1, ntry) {

    # s1 - the adjacency matrix of a directed network ntry - (optional) the
    # number of rewiring steps. If none is given ntry=4*(# of edges in the
    # network) Output: dir_srand - the adjacency matrix of a randomized network
    # with the same set of in- and out-degrees as the original one
    dir_srand <- s1
    ij_srand <- which(dir_srand > 0, arr.ind = TRUE)
    i_srand <- ij_srand[, 1]
    j_srand <- ij_srand[, 2]
    Ne <- nrow(ij_srand)
    nargin <- 1
    if (nargin < 2) {
        ntry <- 4 * Ne
    }
    nrew <- 0
    for (l in 1:ntry) {
        e1 <- 1 + floor(Ne * runif(1, 0, 1))
        e2 <- 1 + floor(Ne * runif(1, 0, 1))
        v1 <- i_srand[e1]
        v2 <- j_srand[e1]

        v3 <- i_srand[e2]
        v4 <- j_srand[e2]
        if ((v1 != v3) & (v1 != v4) & (v2 != v4) & (v2 != v3)) {
            if ((dir_srand[v1, v4] == 0) & (dir_srand[v3, v2] == 0)) {
                dir_srand[v1, v4] <- dir_srand[v1, v2]
                dir_srand[v3, v2] <- dir_srand[v3, v4]
                dir_srand[v1, v2] <- 0
                dir_srand[v3, v4] <- 0
                nrew <- nrew + 1
                i_srand[e1] <- v1
                j_srand[e1] <- v4
                i_srand[e2] <- v3
                j_srand[e2] <- v2
            }
        }
    }

    return(dir_srand)
}


## Internal function csn_edge from CSmiR methods
## (https://github.com/zhangjunpeng411/CSmiR) with GPL-2 license.
csn_edge <- function(gx, gy, boxsize = 0.1) {

    # Define the neighborhood of each plot
    n <- length(gx)
    upper <- pracma::zeros(1, n)
    lower <- pracma::zeros(1, n)
    a <- pracma::zeros(2, n)
    B <- list()

    for (i in seq_len(2)) {
        g <- gx * (i == 1) + gy * (i == 2)
        s1 <- sort(g, index.return = TRUE)[[1]]
        s2 <- sort(g, index.return = TRUE)[[2]]
        n0 <- n - sum(sign(s1))
        h <- round(boxsize/2 * sum(sign(s1)) + pracma::eps(1))
        k <- 1
        while (k <= n) {
            s <- 0
            while ((n >= k + s + 1) && (s1[k + s + 1] == s1[k])) {
                s <- s + 1
            }
            if (s >= h) {
                upper[s2[k:(k + s)]] <- g[s2[k]]
                lower[s2[k:(k + s)]] <- g[s2[k]]
            } else {
                upper[s2[k:(k + s)]] <- g[s2[min(n, k + s + h)]]
                lower[s2[k:(k + s)]] <- g[s2[max(n0 * (n0 > h) + 1, k - h)]]
            }
            k <- k + s + 1
        }

        B[[i]] <- (do.call(cbind, lapply(seq_len(n), function(i) g <= upper[i]))) &
            (do.call(cbind, lapply(seq_len(n), function(i) g >= lower[i])))
        a[i, ] <- colSums(B[[i]])
    }

    # Calculate the normalized statistic of edge gx-gy
    res <- (colSums(B[[1]] & B[[2]]) * n - a[1, ] * a[2, ])/sqrt(a[1, ] * a[2, ] *
        (n - a[1, ]) * (n - a[2, ])/(n - 1) + pracma::eps(1))

    return(res)
}

## Internal function cluster from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster <- function(graph, method = "MCL", expansion = 2, inflation = 2, hcmethod = "average",
    directed = FALSE, outfile = NULL, ...) {
    # method<-match.arg(method)
    if (method == "FN") {
        graph <- simplify(graph)
        fc <- fastgreedy.community(graph, merges = TRUE, modularity = TRUE)
        membership <- membership(fc)
        if (!is.null(V(graph)$name)) {
            names(membership) <- V(graph)$name
        }
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), outfile = outfile)
        } else {
            return(membership)
        }
    } else if (method == "LINKCOMM") {
        edgelist <- get.edgelist(graph)
        if (!is.null(E(graph)$weight)) {
            edgelist <- cbind(edgelist, E(graph)$weight)
        }
        lc <- getLinkCommunities(edgelist, plot = FALSE, directed = directed, hcmethod = hcmethod)
        if (!is.null(outfile)) {
            cluster.save(lc$nodeclusters, outfile = outfile)
        } else {
            return(lc$nodeclusters)
        }
    } else if (method == "MCL") {
        adj <- matrix(rep(0, length(V(graph))^2), nrow = length(V(graph)), ncol = length(V(graph)))
        for (i in seq_along(V(graph))) {
            neighbors <- neighbors(graph, v = V(graph)$name[i], mode = "all")
            j <- match(neighbors$name, V(graph)$name, nomatch = 0)
            adj[i, j] = 1
        }
        lc <- mcl(adj, addLoops = TRUE, expansion = expansion, inflation = inflation,
            allow1 = TRUE, max.iter = 100, ESM = FALSE)
        lc$name <- V(graph)$name
        lc$Cluster <- lc$Cluster

        if (!is.null(outfile)) {
            cluster.save(cbind(lc$name, lc$Cluster), outfile = outfile)
        } else {
            result <- lc$Cluster
            names(result) <- V(graph)$name
            return(result)
        }
    } else if (method == "MCODE") {
        compx <- mcode(graph, vwp = 0.9, haircut = T, fluff = T, fdt = 0.1)
        index <- which(!is.na(compx$score))
        membership <- rep(0, vcount(graph))
        for (i in seq_along(index)) {
            membership[compx$COMPLEX[[index[i]]]] <- i
        }
        if (!is.null(V(graph)$name))
            names(membership) <- V(graph)$name
        if (!is.null(outfile)) {
            cluster.save(cbind(names(membership), membership), outfile = outfile)
            invisible(NULL)
        } else {
            return(membership)
        }
    }
}

## Internal function cluster.save from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster.save <- function(membership, outfile) {
    wd <- dirname(outfile)
    wd <- ifelse(wd == ".", paste(wd, "/", sep = ""), wd)
    filename <- basename(outfile)
    if ((filename == "") || (grepl(":", filename))) {
        filename <- "membership.txt"
    } else if (grepl("\\.", filename)) {
        filename <- sub("\\.(?:.*)", ".txt", filename)
    }
    write.table(membership, file = paste(wd, filename, sep = "/"), row.names = FALSE,
        col.names = c("node", "cluster"), quote = FALSE)
}

## Internal function mcode.vertex.weighting from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.vertex.weighting <- function(graph, neighbors) {
    stopifnot(is_igraph(graph))
    weight <- lapply(seq_len(vcount(graph)), function(i) {
        subg <- induced.subgraph(graph, neighbors[[i]])
        core <- graph.coreness(subg)
        k <- max(core)
        ### k-coreness
        kcore <- induced.subgraph(subg, which(core == k))
        if (vcount(kcore) > 1) {
            if (any(is.loop(kcore))) {
                k * ecount(kcore)/choose(vcount(kcore) + 1, 2)
            } else {
                k * ecount(kcore)/choose(vcount(kcore), 2)
            }
        } else {
            0
        }
    })

    return(unlist(weight))
}

## Internal function mcode.find.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complex <- function(neighbors, neighbors.indx, vertex.weight, vwp, seed.vertex,
    seen) {

    res <- .C("complex", as.integer(neighbors), as.integer(neighbors.indx), as.single(vertex.weight),
        as.single(vwp), as.integer(seed.vertex), seen = as.integer(seen), COMPLEX = as.integer(rep(0,
            length(seen))), PACKAGE = "S2N")

    return(list(seen = res$seen, COMPLEX = which(res$COMPLEX != 0)))
}

## Internal function mcode.find.complexex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complexex <- function(graph, neighbors, vertex.weight, vwp) {
    seen <- rep(0, vcount(graph))

    neighbors <- lapply(neighbors, function(item) {
        item[-1]
    })
    neighbors.indx <- cumsum(unlist(lapply(neighbors, length)))

    neighbors.indx <- c(0, neighbors.indx)
    neighbors <- unlist(neighbors) - 1

    COMPLEX <- list()
    n <- 1
    w.order <- order(vertex.weight, decreasing = TRUE)
    for (i in w.order) {
        if (!(seen[i])) {
            res <- mcode.find.complex(neighbors, neighbors.indx, vertex.weight, vwp,
                i - 1, seen)
            if (length(res$COMPLEX) > 1) {
                COMPLEX[[n]] <- res$COMPLEX
                seen <- res$seen
                n <- n + 1
            }
        }
    }
    rm(neighbors)
    return(list(COMPLEX = COMPLEX, seen = seen))
}

## Internal function mcode.fluff.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.fluff.complex <- function(graph, vertex.weight, fdt = 0.8, complex.g, seen) {
    seq_complex.g <- seq_along(complex.g)
    for (i in seq_complex.g) {
        node.neighbor <- unlist(neighborhood(graph, 1, complex.g[i]))
        if (length(node.neighbor) > 1) {
            subg <- induced.subgraph(graph, node.neighbor)
            if (graph.density(subg, loops = FALSE) > fdt) {
                complex.g <- c(complex.g, node.neighbor)
            }
        }
    }

    return(unique(complex.g))
}

## Internal function mcode.post.process from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.post.process <- function(graph, vertex.weight, haircut, fluff, fdt = 0.8, set.complex.g,
    seen) {
    indx <- unlist(lapply(set.complex.g, function(complex.g) {
        if (length(complex.g) <= 2)
            0 else 1
    }))
    set.complex.g <- set.complex.g[indx != 0]
    set.complex.g <- lapply(set.complex.g, function(complex.g) {
        coreness <- graph.coreness(induced.subgraph(graph, complex.g))
        if (fluff) {
            complex.g <- mcode.fluff.complex(graph, vertex.weight, fdt, complex.g,
                seen)
            if (haircut) {
                ## coreness needs to be recalculated
                coreness <- graph.coreness(induced.subgraph(graph, complex.g))
                complex.g <- complex.g[coreness > 1]
            }
        } else if (haircut) {
            complex.g <- complex.g[coreness > 1]
        }
        return(complex.g)
    })
    set.complex.g <- set.complex.g[lapply(set.complex.g, length) > 2]
    return(set.complex.g)
}

## Internal function mcode from ProNet package (https://github.com/cran/ProNet)
## with GPL-2 license.
mcode <- function(graph, vwp = 0.5, haircut = FALSE, fluff = FALSE, fdt = 0.8, loops = TRUE) {
    stopifnot(is_igraph(graph))
    if (vwp > 1 | vwp < 0) {
        stop("vwp must be between 0 and 1")
    }
    if (!loops) {
        graph <- simplify(graph, remove.multiple = FALSE, remove.loops = TRUE)
    }
    neighbors <- neighborhood(graph, 1)
    W <- mcode.vertex.weighting(graph, neighbors)
    res <- mcode.find.complexex(graph, neighbors = neighbors, vertex.weight = W,
        vwp = vwp)
    COMPLEX <- mcode.post.process(graph, vertex.weight = W, haircut = haircut, fluff = fluff,
        fdt = fdt, res$COMPLEX, res$seen)
    score <- unlist(lapply(COMPLEX, function(complex.g) {
        complex.g <- induced.subgraph(graph, complex.g)
        if (any(is.loop(complex.g)))
            score <- ecount(complex.g)/choose(vcount(complex.g) + 1, 2) * vcount(complex.g) else score <- ecount(complex.g)/choose(vcount(complex.g), 2) * vcount(complex.g)
        return(score)
    }))
    order_score <- order(score, decreasing = TRUE)
    return(list(COMPLEX = COMPLEX[order_score], score = score[order_score]))
}

## Internal function integer.edgelist from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
integer.edgelist <- function(network)
  # Returns an edge list with integer numbers replacing elements of the network.
{
  if(!is.character(network)){
    cn <- cbind(as.character(network[,1]),as.character(network[,2]))
  }else{
    cn <- network
  }
  nodes <- unique(as.character(t(cn)))
  ids <- seq(nodes)
  names(ids) <- nodes
  g <- matrix(ids[t(cn)],nrow(cn),ncol(cn),byrow=TRUE)
  ret <- list()
  ret$edges <- g
  ret$nodes <- ids
  return(ret)
}

## Internal function edge.duplicates from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
edge.duplicates <- function(network, verbose = TRUE)
  # Finds and removes loops, duplicate edges, and bi-directional edges.
{
  xx <- cbind(as.character(network[,1]),as.character(network[,2]))
  edges <- integer.edgelist(network)$edges
  ne <- nrow(edges)
  loops <- rep(0,ne)
  dups <- rep(0,ne)

  out <- .C("edgeDuplicates",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(ne), loops = as.integer(loops), dups = as.integer(dups), as.logical(verbose), PACKAGE = "S2N")

  if(verbose){cat("\n")}

  loops <- which(out$loops == 1)
  dups <- which(out$dups == 1)
  inds <- unique(c(loops,dups))
  ret <- list()
  ret$inds <- inds
  if(length(inds)>0){
    ret$edges <- xx[-inds,]
  }else{
    ret$edges <- xx
  }
  if(verbose){
    if(length(loops)>0){
      cat("   Found and removed ",length(loops)," loop(s)\n",sep="")
    }
    if(length(dups)>0){
      cat("   Found and removed ",length(dups)," duplicate edge(s)\n",sep="")
    }
  }
  return(ret)
}

## Internal function getLinkCommunities from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
getLinkCommunities <- function(network, hcmethod = "average", use.all.edges = FALSE, edglim = 10^4, directed = FALSE, dirweight = 0.5, bipartite = FALSE, dist = NULL, plot = TRUE, check.duplicates = TRUE, removetrivial = TRUE, verbose = TRUE)
  # network is an edge list. Nodes can be ASCII names or integers, but are always treated as character names in R.
  # If plot is true (default), a dendrogram and partition density score as a function of dendrogram height are plotted side-by-side.
  # When there are more than "edglim" edges, hierarchical clustering is carried out via temporary files written to disk using compiled C++ code.
{

  if(is.character(network) && !is.matrix(network)){
    if(file.access(network) == -1){
      stop(cat("\nfile not found: \"",network,"\"\n",sep=""))
    }else{
      network <- read.table(file = network, header = FALSE)
    }
  }
  x <- network
  rm(network)

  if(ncol(x)==3){
    wt <- as.numeric(as.character(x[,3]))
    if(length(which(is.na(wt)==TRUE))>0){
      stop("\nedge weights must be numerical values\n")
    }
    x <- cbind(as.character(x[,1]),as.character(x[,2]))
  }else if(ncol(x)==2){
    x <- cbind(as.character(x[,1]),as.character(x[,2]))
    wt <- NULL
  }else{
    stop("\ninput data must be an edge list with 2 or 3 columns\n")
  }

  if(check.duplicates){
    dret <- edge.duplicates(x, verbose = verbose)
    x <- dret$edges
    if(!is.null(wt)){
      if(length(dret$inds) > 0){
        wt <- wt[-dret$inds]
      }
    }
    rm(dret)
  }

  el <- x # Modified edge list returned to user.
  rm(x)
  len <- nrow(el) # Number of edges.
  nnodes <- length(unique(c(as.character(el[,1]),as.character(el[,2])))) # Number of nodes.

  intel <- integer.edgelist(el) # Edges with numerical node IDs.
  edges <- intel$edges
  node.names <- names(intel$nodes)
  numnodes <- length(node.names)

  if(bipartite){
    # Check that network is bipartite.
    big <- graph_from_edgelist(as.matrix(el), directed = directed)
    bip.test <- bipartite.mapping(big)
    if(!bip.test$res){
      stop("\nnetwork is not bi-partite; change bipartite argument to FALSE\n")
    }
    bip <- rep(1,length(bip.test$type))
    bip[which(bip.test$type==FALSE)] <- 0
    names(bip) <- V(big)$name
    bip <- bip[match(node.names, names(bip))]
    rm(big, bip.test)
  }else{
    bip <- 0
  }

  rm(intel)

  # Switch depending on size of network.
  if(len <= edglim){
    disk <- FALSE
    if(is.null(dist)){
      emptyvec <- rep(1,(len*(len-1))/2)
      if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
      if(!use.all.edges){
        dissvec <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose), PACKAGE = "S2N")$dissvec
      }else{
        dissvec <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = as.double(emptyvec), as.logical(bipartite), as.logical(verbose), PACKAGE = "S2N")$dissvec
      }
      distmatrix <- matrix(1,len,len)
      distmatrix[lower.tri(distmatrix)] <- dissvec
      colnames(distmatrix) <- 1:len
      rownames(distmatrix) <- 1:len
      distobj <- as.dist(distmatrix) # Convert into 'dist' object for hclust.
      rm(distmatrix)
    }else{
      # Did the user provide an adequate distance matrix?
      if(!inherits(dist,"dist")){
        stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
      }else if(attr(dist,which="Size") != len){
        stop("\ndistance matrix size must equal the number of edges in the input network\n")
      }else if(length(dist) != (len*(len-1))/2){
        stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
      }
      distobj <- dist
    }
    if(verbose){
      cat("\n   Hierarchical clustering of edges...")
    }
    #if(hcmethod=="energy"){
    #	hcedges <- energy.hclust(distobj)
    #}else{
    #	hcedges <- hclust(distobj, method = hcmethod)
    #	}
    hcedges <- hclust(distobj, method = hcmethod)
    hcedges$order <- rev(hcedges$order)
    rm(distobj)
    if(verbose){cat("\n")}
  }else{
    disk <- TRUE
    if(!is.null(wt)){ weighted <- TRUE}else{ wt <- 0; weighted <- FALSE}
    if(!use.all.edges){
      rowlen <- .C("getEdgeSimilarities",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),rowlen=integer(len-1),weights=as.double(wt),as.logical(directed),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose), PACKAGE = "S2N")$rowlen
    }else{
      rowlen <- .C("getEdgeSimilarities_all",as.integer(edges[,1]),as.integer(edges[,2]),as.integer(len),as.integer(numnodes),rowlen=integer(len-1),weights=as.double(wt),as.logical(FALSE),as.double(dirweight),as.logical(weighted),as.logical(disk), dissvec = double(1), as.logical(bipartite), as.logical(verbose), PACKAGE = "S2N")$rowlen
    }
    if(verbose){cat("\n")}
    hcobj <- .C("hclustLinkComm",as.integer(len),as.integer(rowlen),heights = single(len-1),hca = integer(len-1),hcb = integer(len-1), as.logical(verbose), PACKAGE = "S2N")
    if(verbose){cat("\n")}
    hcedges<-list()
    hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
    hcedges$height <- hcobj$heights

    hcedges$order <- .C("hclustPlotOrder",as.integer(len),as.integer(hcobj$hca),as.integer(hcobj$hcb),order=integer(len), PACKAGE = "S2N")$order
    hcedges$order <- rev(hcedges$order)
    hcedges$method <- "single"
    class(hcedges) <- "hclust"

  }

  # Calculate link densities, cut the tree, and extract optimal clusters.

  hh <- unique(round(hcedges$height, digits = 5)) # Round to 5 digits to prevent numerical instability affecting community formation.
  countClusters <- function(x,ht){return(length(which(ht==x)))}
  clusnums <- sapply(hh, countClusters, ht = round(hcedges$height, digits = 5)) # Number of clusters at each height.
  numcl <- length(clusnums)

  ldlist <- .C("getLinkDensities",as.integer(hcedges$merge[,1]), as.integer(hcedges$merge[,2]), as.integer(edges[,1]), as.integer(edges[,2]), as.integer(len), as.integer(clusnums), as.integer(numcl), pdens = double(length(hh)), heights = as.double(hh), pdmax = double(1), csize = integer(1), as.logical(removetrivial), as.logical(bipartite), as.integer(bip), as.logical(verbose), PACKAGE = "S2N")

  pdens <- c(0,ldlist$pdens)
  heights <- c(0,hh)
  pdmax <- ldlist$pdmax
  csize <- ldlist$csize

  if(csize == 0){
    stop("\nno clusters were found in this network; maybe try a larger network\n")
  }

  if(verbose){
    cat("\n   Maximum partition density = ",max(pdens),"\n")
  }

  # Read in optimal clusters from a file.
  clus <- list()
  for(i in 1:csize){
    if(verbose){
      mes<-paste(c("   Finishing up...1/4... ",floor((i/csize)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    clus[[i]] <- scan(file = "linkcomm_clusters.txt", nlines = 1, skip = i-1, quiet = TRUE)
  }

  file.remove("linkcomm_clusters.txt")

  # Extract nodes for each edge cluster.
  ecn <- data.frame()
  ee <- data.frame()
  lclus <- length(clus)
  for(i in 1:lclus){
    if(verbose){
      mes<-paste(c("   Finishing up...2/4... ",floor((i/lclus)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ee <- rbind(ee,cbind(el[clus[[i]],],i))
    nodes <- node.names[unique(c(edges[clus[[i]],]))]
    both <- cbind(nodes,rep(i,length(nodes)))
    ecn <- rbind(ecn,both)
  }
  colnames(ecn) <- c("node","cluster")
  colnames(ee) <- c("node1","node2","cluster")

  # Extract the node-size of each edge cluster and order largest to smallest.
  ss <- NULL
  unn <- unique(ecn[,2])
  lun <- length(unn)
  for(i in 1:length(unn)){
    if(verbose){
      mes<-paste(c("   Finishing up...3/4... ",floor((i/lun)*100),"%"),collapse="")
      cat(mes,"\r")
      flush.console()
    }
    ss[i] <- length(which(ecn[,2]==unn[i]))
  }
  names(ss) <- unn
  ss <- sort(ss,decreasing=T)

  # Extract the number of edge clusters that each node belongs to.
  unn <- unique(ecn[,1])

  iecn <- as.integer(as.factor(ecn[,1]))
  iunn <- unique(iecn)
  lunn <- length(iunn)
  nrows <- nrow(ecn)

  oo <- rep(0,lunn)

  oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn), counts = as.integer(oo), as.integer(lunn), as.integer(nrows), as.logical(verbose), PACKAGE = "S2N")$counts

  names(oo) <- unn

  if(verbose){cat("\n")}

  pdplot <- cbind(heights,pdens)

  # Add nodeclusters of size 0.
  missnames <- setdiff(node.names,names(oo))
  m <- rep(0,length(missnames))
  names(m) <- missnames
  oo <- append(oo,m)

  all <- list()

  all$numbers <- c(len,nnodes,length(clus)) # Number of edges, nodes, and clusters.
  all$hclust <- hcedges # Return the 'hclust' object. To plot the dendrogram: 'plot(lcobj$hclust,hang=-1)'
  all$pdmax <- pdmax # Partition density maximum height.
  all$pdens <- pdplot # Add data for plotting Partition Density as a function of dendrogram height.
  all$nodeclusters <- ecn # n*2 character matrix of node names and the cluster ID they belong to.
  all$clusters <- clus # Clusters of edge IDs arranged as a list of lists.
  all$edges <- ee # Edges and the clusters they belong to, arranged so we can easily put them into an edge attribute file for Cytoscape.
  all$numclusters <- sort(oo,decreasing=TRUE) # The number of clusters that each node belongs to (named vector where the names are node names).
  all$clustsizes <- ss # Cluster sizes sorted largest to smallest (named vector where names are cluster IDs).
  all$igraph <- graph_from_edgelist(el, directed = directed) # igraph graph.
  all$edgelist <- el # Edge list.
  all$directed <- directed # Logical indicating if graph is directed or not.
  all$bipartite <- bipartite # Logical indicating if graph is bipartite or not.

  class(all) <- "linkcomm"

  if(plot){
    if(verbose){
      cat("   Plotting...\n")
    }
    if(len < 1500){ # Will be slow to plot dendrograms for large networks.
      if(len < 500){
        all <- plot(all, type="summary", droptrivial = removetrivial, verbose = verbose)
      }else{ # Slow to reverse order of large dendrograms.
        all <- plot(all, type="summary", right = FALSE, droptrivial = removetrivial, verbose = verbose)
      }
    }else if(len <= edglim){
      oldpar <- par(no.readonly = TRUE)
      par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))
      plot(hcedges,hang=-1,labels=FALSE)
      abline(pdmax,0,col='red',lty=2)
      plot(pdens,heights,type='n',xlab='Partition Density',ylab='Height')
      lines(pdens,heights,col='blue',lwd=2)
      abline(pdmax,0,col='red',lty=2)
      par(oldpar)
    }else{
      plot(heights,pdens,type='n',xlab='Height',ylab='Partition Density')
      lines(heights,pdens,col='blue',lwd=2)
      abline(v = pdmax,col='red',lwd=2)
    }
  }

  return(all)

}

## Internal function plot.linkcomm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
#' @export
#' @method plot.linkcomm
plot.linkcomm <- function(x, type = "", ...)
  # S3 method for "plot" generic function.
  # x is a "linkcomm" object.
{
  switch(type,
         summary = plotLinkCommSumm(x, ...),
         members = plotLinkCommMembers(x, ...),
         dend = plotLinkCommDend(x, ...),
         graph = plotLinkCommGraph(x, ...),
         commsumm = plotLinkCommSummComm(x, ...)
  )
}

## Internal function plotLinkCommSumm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommSumm <- function(x, col = TRUE, pal = brewer.pal(9,"Set1"), right = TRUE, droptrivial = TRUE, verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  oldpar <- par(no.readonly = TRUE)
  # Set up and colour clusters in dendrogram.
  if(is.null(x$dendr)){
    dd <- as.dendrogram(x$hclust)
  }else{
    dd <- x$dendr
  }
  if(col && is.null(x$dendr)){
    cl <- unlist(x$clusters)
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(x$clusters))
    cols <- sample(cols,length(x$clusters),replace=FALSE)
    numnodes <- nrow(x$hclust$merge) + length(which(x$hclust$merge[,1]<0)) + length(which(x$hclust$merge[,2]<0))
    dd <- dendrapply(dd,.COL,height=x$pdmax,clusters=cl,cols=cols,labels=FALSE,numnodes=numnodes,droptrivial = droptrivial,verbose=verbose)
    if(verbose){cat("\n")}
    assign("i",0,environment(.COL))
    assign("memb",0,environment(.COL))
    assign("first",0,environment(.COL))
    assign("left",0,environment(.COL))
  }
  if(is.null(x$dendr) && right){ dd <- rev(dd)}
  grid.newpage()
  plot.new()
  # Set margin.
  margin<-unit(0.045,"npc")
  pushViewport(viewport(x=margin,y=margin,width=unit(1,"npc")-2*margin,height=unit(1,"npc")-2*margin,just=c("left","bottom")))
  pushViewport(viewport(layout=grid.layout(nrow=4,ncol=4,widths=unit(c(1,0.01,0.2,0.09),units=c("null",rep("native",3))),heights=unit(c(0.04,0.79,0.025,0.05),units=rep("npc",4)),respect=TRUE)))
  # Plot dendrogram using base plot function.
  pushViewport(viewport(layout.pos.row=1:3,layout.pos.col=1))
  #return(gridPLT())
  gpl <- c(0.009,0.71,0.13,0.902)
  par(oma=rep(0,4),mar=rep(0,4),ann=FALSE,omd=c(0,1,0,1),pty="m",mgp=rep(0,3),fig = gpl,xpd=NA,new=TRUE)
  plot(dd,axes=FALSE,leaflab="none")
  popViewport(1)
  # Title.
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=1:3))
  title <- grid.text("Link Communities Dendrogram",x = unit(0.5,"npc"),y = unit(2,"npc"),draw=FALSE,name="title")
  title <- editGrob(title,gp = gpar(fontsize=14))
  grid.draw(title)
  popViewport(1)
  # Plot link partition densities.
  numzeros <- -1*log10(max(x$pdens[,2]))
  if(numzeros <= 1){ # Prevent partition density axis from being rounded to 0.
    rr <- 1
  }else{
    rr <- trunc(numzeros)+1
  }
  if(round(max(x$pdens[,2]), digits = rr) > max(x$pdens[,2])){
    xscale_add <- round(max(x$pdens[,2]), digits = rr) + 0.05*max(x$pdens[,2]) # Add 5% to part density x-axis.
    xaxs_max <- round(max(x$pdens[,2]), digits = rr)
  }else{
    xscale_add <- max(x$pdens[,2]) + 0.075*max(x$pdens[,2]) # 7.5% of max partition density added to x-axis.
    xaxs_max <- round(max(x$pdens[,2]), digits = (rr+1))
  }
  pushViewport(viewport(layout.pos.row=2, layout.pos.col=3, xscale=c(0,xscale_add),yscale=c(0,1)))
  ph <- x$pdens[,1]/max(x$pdens[,1])
  max <- x$pdmax/max(x$pdens[,1])

  grid.lines(x$pdens[,2],ph,gp=gpar(col='blue',lwd=2),default.units="native")
  xticks <- seq(0, xaxs_max, length.out=3)
  xa <- grid.xaxis(at = xticks,draw=FALSE,name="xa")
  xa <- editGrob(xa,gPath="ticks",y1 = unit(-0.02,"npc"))
  xa <- editGrob(xa,gPath="labels",gp = gpar(fontsize=10),y = unit(-0.04,"npc"))
  grid.draw(xa)
  xl <- grid.text("Partition Density",x=unit(0.5, "npc"), y = unit(-0.08, "npc"),draw=FALSE,name="xl")
  xl <- editGrob(xl,gp = gpar(fontsize=10))
  grid.draw(xl)
  popViewport(1)
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=1:3))
  grid.lines(x=c(0,1),y=c(max,max),gp = gpar(col="red",lty=2,lwd=2))
  popViewport(1)
  pushViewport(viewport(layout.pos.row=2,layout.pos.col=4))
  yticks <- c(0,0.2,0.4,0.6,0.8,1)
  ya <- grid.yaxis(at = yticks,name="ya",draw=FALSE)
  ya <- editGrob(ya,gPath="ticks",x1 = unit(0.1,"npc"))
  ya <- editGrob(ya,gPath="labels",gp = gpar(fontsize=10),x = unit(0.2,"npc"))
  ya <- editGrob(ya,gPath="labels",just = c("left","centre"))
  if(max(x$pdens[,1]<1)){
    yu <- seq(0,round(max(x$pdens[,1]),2),length.out=6)
    roundS <- function(x){return(round(x,2))}
    yu <- sapply(yu,roundS)
    ya <- editGrob(ya,gPath="labels",label=as.character(yu))
  }
  grid.draw(ya)
  yl <- grid.text("Height",x=unit(0.8, "npc"), y = unit(0.5, "npc"),rot=90,draw=FALSE,name="yl")
  yl <- editGrob(yl,gp = gpar(fontsize=10))
  grid.draw(yl)
  popViewport(1)
  # Summary statistics.
  pushViewport(viewport(layout.pos.row=4,layout.pos.col=1))
  summ <- paste("# edges = ",x$numbers[1],",   ","# nodes = ",x$numbers[2],"\n# clusters = ",x$numbers[3],",   Largest cluster = ",x$clustsizes[1]," nodes\nHclust method: ",x$hclust$method)
  ne <- grid.text(summ,x=unit(0.5,"npc"),y=unit(0.1,"npc"),draw=FALSE,name="ne")
  ne <- editGrob(ne,gp = gpar(fontsize=11))
  grid.draw(ne)
  popViewport(1)
  # Return linkcomm object with dendrogram so we don't have to generate it again in the future.
  if(is.null(x$dendr)){
    x$dendr <- dd
    return(x)
  }
  popViewport(0)
  par(oldpar)
}

## Internal function plotLinkCommMembers from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommMembers <- function(x, nodes = head(names(x$numclusters),10), pal = brewer.pal(11,"Spectral"), shape = "rect", total=TRUE, fontsize=11, nspace = 3.5, maxclusters = 20)
  # Plots community membership matrix using a community-specific colour scheme.
  # x is a "linkcomm" object.
{
  # Construct community matrix.
  comms <- unique(x$nodeclusters[as.character(x$nodeclusters[,1])%in%nodes,2]) # Community (cluster) IDs.
  if(length(comms) > maxclusters){
    comms <- comms[1:maxclusters]
  }
  commatrix <- getCommunityMatrix(x,nodes=nodes)
  crf <- colorRampPalette(pal,bias=1)
  cols <- crf(length(comms))
  grid.newpage()
  # Set margin.
  if(total){
    C <- 2; R <- 3
    nodesums <- apply(commatrix,1,sum)
    commsums <- apply(commatrix,2,sum)
  }else{
    C <- 1; R <- 2
  }
  margin<-unit(0.1,"lines")
  pushViewport(viewport(x=1,y=1,width=unit(1,"npc")-2*margin,height=unit(1,"npc")-2*margin,just=c("right","top")))
  pushViewport(viewport(layout=grid.layout(nrow=length(nodes)+R,ncol=length(comms)+C,widths=unit(c(nspace,rep(1,length(comms)+C-1)),rep("null",length(comms)+C)),heights=unit(rep(1,length(nodes)+R),rep("null",length(nodes)+R)),respect=TRUE)))
  # Titles.
  pushViewport(viewport(layout.pos.row=1,layout.pos.col=2:length(comms)+1))
  ctitle <- grid.text("Community Membership",x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ctitle")
  ctitle <- editGrob(ctitle,gp = gpar(fontsize=14))
  grid.draw(ctitle)
  popViewport(1)
  # Draw membership coloured squares/circles/polygons.
  for(i in 1:(length(nodes)+R-2)){
    if(i != length(nodes)+1){
      pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=1))
      nname <- grid.text(as.character(nodes[i]),x = unit(0.9,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="nname")
      nname <- editGrob(nname,gp = gpar(fontsize=fontsize),just="right")
      grid.draw(nname)
      popViewport(1)
    }
    for(j in 1:(length(comms)+C-1)){
      if(total && j == length(comms)+1 && i != length(nodes)+1){
        pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
        ntot <- grid.text(nodesums[i],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ntot")
        ntot <- editGrob(ntot,gp = gpar(fontsize=12))
        grid.draw(ntot)
        popViewport(1)
        if(i == 1){
          pushViewport(viewport(layout.pos.row=2,layout.pos.col=j+1))
          rt <- grid.text(expression(Sigma),x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="rt")
          rt <- editGrob(rt,gp = gpar(fontsize=12))
          grid.draw(rt)
          popViewport(1)
        }
      }else{
        if(i == 1 && j != length(comms)+1){
          pushViewport(viewport(layout.pos.row=2,layout.pos.col=j+1))
          rtitle <- grid.text(comms[j],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="rtitle")
          rtitle <- editGrob(rtitle,gp = gpar(fontsize=12))
          grid.draw(rtitle)
          popViewport(1)
        }
        if(total && i == length(nodes)+1 && j != length(comms)+1){
          if(j==1){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=1))
            ct <- grid.text(expression(Sigma),x = unit(0.9,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ct")
            ct <- editGrob(ct,gp = gpar(fontsize=12))
            grid.draw(ct)
            popViewport(1)
          }
          pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
          ctot <- grid.text(commsums[j],x = unit(0.5,"npc"),y = unit(0.5,"npc"),draw=FALSE,name="ctot")
          ctot <- editGrob(ctot,gp = gpar(fontsize=12))
          grid.draw(ctot)
          popViewport(1)
        }else if(i != length(nodes)+1 && j != length(comms)+1){
          if(commatrix[i,j] == 1){
            fill <- cols[j]
          }else{
            fill <- "white"
          }
          if(shape=="rect"){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
            grid.rect(gp=gpar(fill=fill,col="grey"),width = unit(0.9,"npc"), height = unit(0.9,"npc"),draw=TRUE)
            popViewport(1)
          }else if(shape=="circle"){
            pushViewport(viewport(layout.pos.row=i+2,layout.pos.col=j+1))
            grid.circle(x=0.5,y=0.5,r=0.45,gp=gpar(fill=fill,col="grey"),draw=TRUE)
            popViewport(1)
          }
        }
      }
    }
  }
}

## Internal function plotLinkCommDend from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommDend <- function(x, col=TRUE, pal = brewer.pal(9,"Set1"), height=x$pdmax, right = FALSE, labels=FALSE, plotcut=TRUE, droptrivial = TRUE, leaflab = "none", verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  dd <- as.dendrogram(x$hclust)
  if(col){
    cl <- unlist(x$clusters)
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(x$clusters))
    cols <- sample(cols,length(x$clusters),replace=FALSE)
    numnodes <- nrow(x$hclust$merge) + length(which(x$hclust$merge[,1]<0)) + length(which(x$hclust$merge[,2]<0))
    dd <- dendrapply(dd, .COL, height=height, clusters=cl, cols=cols, labels=labels, numnodes = numnodes, droptrivial = droptrivial, verbose = verbose)
    cat("\n")
    assign("i",0,environment(.COL))
    assign("memb",0,environment(.COL))
    assign("first",0,environment(.COL))
    assign("left",0,environment(.COL))
  }
  if(right){
    dd <- rev(dd)
  }
  plot(dd,ylab="Height", leaflab = leaflab, ...)
  if(plotcut){
    abline(h=height,col='red',lty=2,lwd=2)
  }
  #ll <- sapply(x$clusters,length)
  #maxnodes <- length(unique(x$nodeclusters[x$nodeclusters[,2]%in%which(ll==max(ll)),1]))
  summ <- paste("# clusters = ",length(x$clusters),"\nLargest cluster = ",x$clustsizes[1]," nodes")
  mtext(summ, line = -28)
}

## Internal function from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
.COL<-local({

  memb <- 0
  first <- 0
  i <- 0
  left <- 0

  colorHclusters <<- function(x, height, clusters, cols, labels, numnodes, droptrivial, verbose)
    # Adds colours to edges that belong to clusters below "height" in the dendrogram.
    # Clusters gives leaf IDs for clusters that should be coloured.
    # x is a node in the tree.
  {
    left <<- left + 1
    if(verbose){
      out <- paste(c("   Colouring dendrogram... ",floor((left/numnodes)*100),"%"),collapse="")
      cat(out,"\r")
      flush.console()
    }

    if(round(attributes(x)$height,digits=5) > height){
      return(x)
    }else{
      if(is.leaf(x)){
        if(is.na(match(as.numeric(attributes(x)$label),clusters))){
          if(!labels){
            attributes(x)$label <- NULL
          }
          return(x)
        }
      }
      a <- attributes(x)
      if(memb == 0){
        memb <<- attributes(x)$members
        if(droptrivial == TRUE && memb == 2){
          memb <<- 0
          first <<- 1
        }else{
          i <<- i+1
          first <<- 1 # Because we don't colour the edge leading to the first node in a cluster.
        }
      }
      if(first == 0){
        attr(x,"edgePar") <- c(a$edgePar,list(col = cols[i], lwd = 2))
      }
      if(is.leaf(x)){
        if(!labels){
          attributes(x)$label <- NULL
        }
        memb <<- memb-1
      }
      first <<- 0
    }

    return(x)
  }

})


## Internal function plotLinkCommGraph from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommGraph <- function(x, clusterids = 1:length(x$clusters), nodes = NULL, layout = layout.fruchterman.reingold, pal = brewer.pal(7,"Set2"), random = TRUE, node.pies = TRUE, pie.local = TRUE, vertex.radius = 0.03, scale.vertices = 0.05, edge.color = NULL, vshape = "none", vsize = 15, ewidth = 3, margin = 0, vlabel.cex = 0.8, vlabel.color = "black", vlabel.family = "Helvetica", vertex.color = "palegoldenrod", vlabel = TRUE, col.nonclusters = "black", jitter = 0.2, circle = TRUE, printids = TRUE, cid.cex = 1, shownodesin = 0, showall = FALSE, verbose = TRUE, ...)
  # x is a "linkcomm" object.
{
  if(length(nodes) > 0){
    clusterids <- which.communities(x, nodes = nodes)
  }
  clusters <- x$clusters[clusterids]
  miss <- setdiff(x$hclust$order,unlist(clusters))
  crf <- colorRampPalette(pal,bias=1)
  cols <- crf(length(clusters))
  if(random){
    cols <- sample(cols,length(clusters),replace=FALSE)
  }
  if(showall){
    # Add single edge "clusters".
    single <- setdiff(1:x$numbers[1],unlist(clusters))
    ll <- length(clusters)
    for(i in 1:length(single)){
      clusters[[(i+ll)]] <- single[i]
    }
    cols <- append(cols, rep(col.nonclusters, length(single)))
  }
  drawcircle <- FALSE
  if(inherits(layout,"character")){
    if(layout == "spencer.circle"){
      if(length(clusters) > length(x$clusters[1:x$numbers[3]])){
        clusterids <- 1:x$numbers[3]
      }
      ord <- orderCommunities(x, clusterids = clusterids, verbose = FALSE)
      clusters <- ord$ordered
      clusterids <- ord$clusids
      layout <- layout.spencer.circle(x, clusterids = clusterids, jitter = jitter, verbose = verbose)$nodes
      drawcircle <- TRUE
    }
  }
  names(cols) <- clusterids
  if(length(unlist(clusters)) < nrow(x$edgelist) || length(miss) == 0){
    # Convert old clus ids into new ones.
    edges <- x$edgelist[unlist(clusters),]
    ig <- graph.edgelist(edges, directed=x$directed)
    clen <- sapply(clusters,length)
    j<-1
    # Colour edges according to community membership.
    for(i in 1:length(clusters)){
      newcids <- j:sum(clen[1:i])
      E(ig)[newcids]$color <- cols[i]
      j <- tail(newcids,1)+1
    }
  }else{
    ig <- x$igraph
    for(i in 1:length(clusters)){
      E(ig)[clusters[[i]]]$color <- cols[i]
    }
  }

  if(shownodesin == 0){
    vnames <- V(ig)$name
  }else{ # Show nodes that belong to more than x number of communities.
    vnames <- V(ig)$name
    inds <- NULL
    for(i in 1:length(vnames)){
      if(x$numclusters[which(names(x$numclusters)==vnames[i])] < shownodesin){
        inds <- append(inds,i)
      }
    }
    vnames[inds] <- ""
  }
  if(vlabel==FALSE){
    vnames = NA
  }

  dev.hold(); on.exit(dev.flush())
  oldpar <- par(no.readonly = TRUE)
  par(mar = c(4,4,2,2))

  if(!node.pies){
    plot(ig, layout=layout, vertex.shape=vshape, edge.width=ewidth, vertex.label=vnames, vertex.label.family=vlabel.family, vertex.label.color=vlabel.color, vertex.size=vsize, vertex.color=vertex.color, margin=margin, vertex.label.cex = vlabel.cex, ...)
  }else{
    nodes <- V(ig)$name
    # Get node community membership by edges.
    if(pie.local){
      edge.memb <- numberEdgesIn(x, clusterids = clusterids, nodes = nodes)
    }else{
      edge.memb <- numberEdgesIn(x, nodes = nodes)
    }

    cat("   Getting node layout...")
    if(inherits(layout,"function")){
      lay <- layout(ig)
    }else{
      lay <- layout
    }
    lay <- layout.norm(lay, xmin=-1, xmax=1, ymin=-1, ymax=1)
    rownames(lay) <- V(ig)$name
    cat("\n")
    node.pies <- .nodePie(edge.memb=edge.memb, layout=lay, nodes=nodes, edges=100, radius=vertex.radius, scale=scale.vertices)
    cat("\n")
    # Plot graph.
    if(is.null(edge.color)){
      plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=1, edge.width=ewidth, vertex.label.color=vlabel.color, ...)
    }else{
      plot(ig, layout=lay, vertex.shape="none", vertex.label=NA, vertex.label.dist=1, edge.width=ewidth, vertex.label.color=vlabel.color, edge.color=edge.color, ...)
    }
    labels <- list()
    # Plot node pies and node names.
    for(i in 1:length(node.pies)){
      yp <- NULL
      for(j in 1:length(node.pies[[i]])){
        seg.col <- cols[which(names(cols)==names(edge.memb[[i]])[j])]
        polygon(node.pies[[i]][[j]][,1], node.pies[[i]][[j]][,2], col = seg.col)
        yp <- append(yp, node.pies[[i]][[j]][,2])
      }
      lx <- lay[which(rownames(lay)==names(node.pies[i])),1] + 0.1
      ly <- max(yp) + 0.02 # Highest point of node pie.
      labels[[i]] <- c(lx, ly)
    }
    # Plot node names after nodes so they overlay them.
    for(i in 1:length(labels)){
      text(labels[[i]][1], labels[[i]][2], labels = vnames[which(nodes==names(node.pies[i]))], cex = vlabel.cex, col = vlabel.color)
    }
  }

  if(circle && drawcircle){
    # Add circle for Spencer layout.
    cx<-NULL; for(i in 1:100){cx[i]<-1.25*cos(i*(2*pi)/100)}
    cy<-NULL; for(i in 1:100){cy[i]<-1.25*sin(i*(2*pi)/100)}
    polygon(cx-0.08,cy-0.08, border="grey",lwd=2)
    # Add community anchor points and cluster IDs.
    for(i in 1:length(clusters)){
      px <- 1.1*cos(i*(2*pi)/length(clusters))
      py <- 1.1*sin(i*(2*pi)/length(clusters))
      points(px-0.08,py-0.08, pch = 20, col = cols[i])
      if(printids){
        tx <- 1.3*cos(i*(2*pi)/length(clusters))
        ty <- 1.3*sin(i*(2*pi)/length(clusters))
        text(tx-0.08,ty-0.08, labels = clusterids[i], col = cols[i], cex = cid.cex, font=2)
      }
    }
  }
  par(oldpar)

}


## Internal function layout.spencer.circle from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
layout.spencer.circle <- function(x, clusterids = 1:x$numbers[3], verbose = TRUE, jitter = 0.2)
  # Returns x-y node coordinates for Rob Spencer's circular layout of link communities together with x-y coordinates for the community anchors.
  # x is a "linkcomm" object.
{
  clusters <- x$clusters[clusterids]
  edges <- x$edgelist[unlist(clusters),]
  ig <- graph.edgelist(edges, directed=FALSE)
  # Put communities in dendrogram order.
  clusters <- orderCommunities(x, clusterids = clusterids, verbose = verbose)$ordered
  # Set up community anchor points in Cartesian coordinates around unit circle (communities evenly spaced).
  xy_anchors <- matrix(0,length(clusters),2)
  for(i in 1:length(clusters)){
    xy_anchors[i,] <- c(cos(i*(2*pi)/length(clusters)), sin(i*(2*pi)/length(clusters)))
  }
  # Calculate community membership percentages per node.
  nodes <- c(x$edgelist[unlist(clusters),1],x$edgelist[unlist(clusters),2])
  node_names <- unique(nodes)
  xy_nodes <- matrix(0,length(node_names),2)
  for(i in 1:length(node_names)){
    if(verbose){
      mes <- paste(c("   Calculating node co-ordinates for Spencer circle...",floor(i/(length(node_names))*100),"%"),collapse="")
      cat(mes,'\r')
      flush.console()
    }
    freqs <- NULL
    total <- length(which(nodes==node_names[i]))
    for(j in 1:length(clusters)){
      freqs <- length(which(c(x$edgelist[clusters[[j]],1],x$edgelist[clusters[[j]],2])==node_names[i]))/total
      # Update x-y coordinates for this node.
      xy_nodes[i,1] <- sum(xy_nodes[i,1], freqs*xy_anchors[j,1])
      xy_nodes[i,2] <- sum(xy_nodes[i,2], freqs*xy_anchors[j,2])
    }
    # Add random jitter if this node has identical x-y coordinates to an earlier node.
    if(duplicated(xy_nodes)[i]){
      xy_nodes[i,1] <- sum(xy_nodes[i,1], runif(1, min = -jitter, max = jitter))
      xy_nodes[i,2] <- sum(xy_nodes[i,2], runif(1, min = -jitter, max = jitter))
    }
  }

  rownames(xy_nodes) <- node_names

  xy_nodes <- xy_nodes[match(V(ig)$name,rownames(xy_nodes)),]
  xy_nodes <- xy_nodes[!is.na(xy_nodes[,1]),]

  if(verbose){cat("\n")}

  xy <- list()
  xy$nodes <- xy_nodes
  xy$anchors <- xy_anchors

  return(xy)

}

## Internal function plotLinkCommSummComm from linkcomm package
## (https://github.com/alextkalinka/linkcomm).
plotLinkCommSummComm <- function(x, clusterids = 1:x$numbers[3], summary = "conn", pie = FALSE, col = TRUE, pal = brewer.pal(11,"Spectral"), random = FALSE, verbose = TRUE, ...)
  # Plots pie or bar chart summarising sizes of communities in terms of nodes, link density, community connectedness, or community modularity.
  # x is a "linkcomm" object.
{
  if(col){
    crf <- colorRampPalette(pal,bias=1)
    cols <- crf(length(clusterids))
    if(random){
      cols <- sample(cols,length(clusterids),replace=FALSE)
    }
  }else{
    cols <- "lightblue"
  }
  # Extract number of nodes per community.
  nums <- NULL
  if(summary == "nodes"){
    for(i in 1:length(clusterids)){
      nums[i] <- length(unique(c(x$edgelist[x$clusters[[clusterids[i]]],1],x$edgelist[x$clusters[[clusterids[i]]],2])))
    }
    main <- "Node density per community"
  }else if(summary == "ld"){
    nums <- LinkDensities(x, clusterids = clusterids)
    main <- "Link density per community"
  }else{
    nums <- getCommunityConnectedness(x, clusterids = clusterids, conn = summary, verbose = verbose)
    if(summary == "conn"){
      main <- "Community Connectedness"
    }else{
      main <- "Community Modularity"
    }
  }
  names(nums) <- clusterids

  if(pie){
    pie(nums, col = cols, main = main, ...)
  }else{
    barplot(nums, xlab = "Community", ylab = main, col = cols)
    abline(h=0)
  }

}


## Internal function csn_edge from CSmiR methods
## (https://github.com/zhangjunpeng411/CSmiR) with GPL-2 license.
CSN_net_bootstrap <- function(reg_cancer, tar_cancer, boxsize = 0.1, bootstrap_betw_point = 5,
    bootstrap_num = 100, p.value.cutoff = 0.05) {

    regs_num <- ncol(reg_cancer)
    tars_num <- ncol(tar_cancer)
    cell_num <- nrow(reg_cancer)
    bootstrap_sample <- lapply(seq(bootstrap_num), function(i) sample(seq(cell_num),
        bootstrap_betw_point * (cell_num - 1), replace = TRUE))
    reg_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(reg_cancer, reg_cancer[bootstrap_sample[[i]],
        ]))
    tar_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(tar_cancer, tar_cancer[bootstrap_sample[[i]],
        ]))
    res <- matrix(NA, nrow = regs_num * tars_num, ncol = cell_num + 2)
    for (i in seq(regs_num)) {
        for (j in seq(tars_num)) {
            res[(i - 1) * tars_num + j, 1] <- colnames(reg_cancer)[i]
            res[(i - 1) * tars_num + j, 2] <- colnames(tar_cancer)[j]
            res[(i - 1) * tars_num + j, 3:(cell_num + 2)] <- do.call(WGCNA::pmedian,
                lapply(seq(bootstrap_num), function(k) {
                  csn_edge(reg_bootstrap[[k]][, i], tar_bootstrap[[k]][, j], boxsize = boxsize)[seq(cell_num)]
                }))
        }
    }

    q <- -qnorm(p.value.cutoff)
    res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i +
        2]) > q), seq(2)])

    return(res_list)
}

## Internal function csn_edge from CSmiR methods
## (https://github.com/zhangjunpeng411/CSmiR) with GPL-2 license.
CSN_net <- function(reg_cancer, tar_cancer, boxsize = 0.1, p.value.cutoff = 0.05,
    num.cores = 2) {

    regs_num <- ncol(reg_cancer)
    tars_num <- ncol(tar_cancer)
    cell_num <- nrow(reg_cancer)
    int_num <- regs_num * tars_num
    res <- matrix(NA, nrow = int_num, ncol = cell_num + 2)
    index <- matrix(NA, nrow = int_num, ncol = 2)

    for (i in seq(regs_num)) {
        for (j in seq(tars_num)) {
            index[(i - 1) * tars_num + j, 1] <- i
            index[(i - 1) * tars_num + j, 2] <- j
        }
    }

    # get number of cores to run
    cl <- makeCluster(num.cores)
    registerDoParallel(cl)
    interin <- foreach(k = seq(int_num), .packages = c("pracma"), .export = c("csn_edge")) %dopar%
        {
            csn_edge(reg_cancer[, index[k, 1]], tar_cancer[, index[k, 2]], boxsize = boxsize)
        }
    # shut down the workers
    stopCluster(cl)
    stopImplicitCluster()

    interin <- do.call(rbind, interin)

    for (i in seq(int_num)) {
        res[i, 1] <- colnames(reg_cancer)[index[i, 1]]
        res[i, 2] <- colnames(tar_cancer)[index[i, 2]]
    }

    res[, 3:(cell_num + 2)] <- interin

    q <- -qnorm(p.value.cutoff)
    res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i +
        2]) > q), seq(2)])

    return(res_list)
}

## Internal Function of network inference using Pearson
## (https://github.com/zhangjunpeng411/Scan) ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs p.value: A cutoff of p-value
## Output: adj.Matrix is a zero-one matrix between ncRNAs and mRNAs
Pearson_adj <- function(ncRExp, mRExp, p.value = 0.05) {

    cor.pvalue <- corAndPvalue(t(mRExp), t(ncRExp), method = "pearson")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

    return(adj.Matrix)
}

## Internal Function of network inference using Spearman
## (https://github.com/zhangjunpeng411/Scan) ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs p.value: A cutoff of p-value
## Output: adj.Matrix is a zero-one matrix between ncRNAs and mRNAs
Spearman_adj <- function(ncRExp, mRExp, p.value = 0.05) {

    cor.pvalue <- corAndPvalue(t(mRExp), t(ncRExp), method = "spearman")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

    return(adj.Matrix)
}

## Internal Function of network inference using Kendall
## (https://github.com/zhangjunpeng411/Scan) ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs p.value: A cutoff of p-value
## Output: adj.Matrix is a zero-one matrix between ncRNAs and mRNAs
Kendall_adj <- function(ncRExp, mRExp, p.value = 0.05) {

    cor.pvalue <- corAndPvalue(t(mRExp), t(ncRExp), method = "kendall")$p
    adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

    return(adj.Matrix)
}

## Internal Function for calculating z-score of a matrix mat: Input matrix
## Output: mat.zscore is the transformed zscore matrix
matrixzscore <- function(mat) {

    mat.mean <- mean(mat[!is.na(mat)])
    mat.sd <- sd(mat[!is.na(mat)])
    mat.zscore <- (mat - mat.mean)/mat.sd

    return(mat.zscore)
}

## Internal Function of network inference using Pearson ncRExp: Gene expression
## values of ncRNAs mRExp: Gene expression values of mRNAs Output: adj.Matrix
## is a correlation matrix between ncRNAs and mRNAs
Pearson <- function(ncRExp, mRExp) {

    adj.Matrix <- corAndPvalue(mRExp, ncRExp, method = "pearson")$cor

    return(adj.Matrix)
}

## Internal Function for inferring the significant p-values between ncRNAs and
## mRNAs in a specific sample of interest ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs sample_index: the index of
## specific sample Output: Pearson.pvalue is the significant p-values between
## ncRNAs and mRNAs in a specific sample of interest
Pearson_Scan <- function(ncRExp, mRExp, sample_index) {

    nsamples <- nrow(ncRExp)
    original_all <- Pearson(ncRExp, mRExp)
    original_single <- Pearson(ncRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Pearson.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Pearson.pvalue)
}

## Internal Function of network inference using Spearman ncRExp: Gene
## expression values of ncRNAs mRExp: Gene expression values of mRNAs Output:
## adj.Matrix is a correlation matrix between ncRNAs and mRNAs
Spearman <- function(ncRExp, mRExp) {

    adj.Matrix <- corAndPvalue(mRExp, ncRExp, method = "spearman")$cor

    return(adj.Matrix)
}

## Internal Function for inferring the significant p-values between ncRNAs and
## mRNAs in a specific sample of interest ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs sample_index: the index of
## specific sample Output: Spearman.pvalue is the significant p-values between
## ncRNAs and mRNAs in a specific sample of interest
Spearman_Scan <- function(ncRExp, mRExp, sample_index) {

    nsamples <- nrow(ncRExp)
    original_all <- Spearman(ncRExp, mRExp)
    original_single <- Spearman(ncRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Spearman.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Spearman.pvalue)
}

## Internal Function of network inference using Kendall ncRExp: Gene expression
## values of ncRNAs mRExp: Gene expression values of mRNAs Output: adj.Matrix
## is a correlation matrix between ncRNAs and mRNAs
Kendall <- function(ncRExp, mRExp) {

    adj.Matrix <- corAndPvalue(mRExp, ncRExp, method = "kendall")$cor

    return(adj.Matrix)
}

## Internal Function for inferring the significant p-values between ncRNAs and
## mRNAs in a specific sample of interest ncRExp: Gene expression values of
## ncRNAs mRExp: Gene expression values of mRNAs sample_index: the index of
## specific sample Output: Kendall.pvalue is the significant p-values between
## ncRNAs and mRNAs in a specific sample of interest
Kendall_Scan <- function(ncRExp, mRExp, sample_index) {

    nsamples <- nrow(ncRExp)
    original_all <- Kendall(ncRExp, mRExp)
    original_single <- Kendall(ncRExp[-sample_index, ], mRExp[-sample_index, ])
    original <- nsamples * original_all - (nsamples - 1) * original_single
    Kendall.pvalue <- 1 - pnorm(abs(matrixzscore(original)))

    return(Kendall.pvalue)
}

## Internal Function for geting upper and lower bounds function
## Calculate upper and lower bounds for each gene in each cell
get_upper_and_lower_file <- function(data, boxsize) {

    geneNum <- nrow(data)   # Number of genes
    cellNum <- ncol(data)   # Number of cells

    # Initialize upper and lower bound matrices
    upper <- matrix(0, nrow = geneNum, ncol = cellNum)
    lower <- matrix(0, nrow = geneNum, ncol = cellNum)

    # Process each gene
    for (i in 1:geneNum) {
      # Get sorted expression values and their indices for current gene
      s2 <- order(data[i, ])
      s1 <- sort(data[i, ])

      # Count cells with expression value > 0
      sign_c <- sum(s1 > 0)
      n0 <- cellNum - sign_c  # Count cells with expression value = 0

      # Calculate window size
      h <- round(boxsize / 2 * sign_c)
      k <- 1  # Starting index

      # Iterate through all cells
      while (k <= cellNum) {
        s <- 0  # Counter for cells with same expression value

        # Count consecutive cells with same expression value
        while (k + s + 1 <= cellNum && s1[k + s + 1] == s1[k]) {
          s <- s + 1
        }

        # Set upper and lower bounds based on count of same expression values
        if (s >= h) {
          # If count of same expression values >= window size, set both bounds to current value
          upper[i, s2[k:(k + s)]] <- data[i, s2[k]]
          lower[i, s2[k:(k + s)]] <- data[i, s2[k]]
        } else {
          # Otherwise, set upper bound to value h positions ahead, lower bound to value h positions behind
          upper[i, s2[k:(k + s)]] <- data[i, s2[min(cellNum, k + s + h)]]
          lower[i, s2[k:(k + s)]] <- data[i, s2[max(1, k - h)]]
        }

        # Move to next different expression value
        k <- k + s + 1
      }
    }

    # Return upper and lower bound matrices
    return(list(upper = upper, lower = lower))
}

## Internal Function for Z-score conversion
## Convert mutual information values to Z-scores
convert_zScore <- function(mutual_info) {
  # Calculate row means and standard deviations
    mu <- rowMeans(mutual_info)
    std <- apply(mutual_info, 1, sd)

    # Apply Z-score formula: (x - mu) / std
    mutual_info <- (mutual_info - mu) / std

    # Handle NaN values
    mutual_info[is.na(mutual_info)] <- 0

    # Round to 2 decimal places
    mutual_info <- round(mutual_info, 2)

    return(mutual_info)
}

## Internal Function for calculating gene entropy matrix
## Calculate entropy for each gene
get_entropyX_matrix <- function(Xt, geneTot, cellTot, binsNum) {
    # Initialize frequency matrix
    freq <- matrix(0, nrow = geneTot, ncol = binsNum)

    # Calculate frequency distribution of discretized values for each gene
    for (r in 1:geneTot) {
      # Use factor to ensure all bins are counted, even if some have no values
      bin_counts <- table(factor(Xt[r, ], levels = 1:binsNum))
      freq[r, ] <- bin_counts
    }

    # Calculate probabilities
    prob <- freq / cellTot

    # Calculate entropy
    EntropyX <- prob * log2(prob)

    # Handle NaN values (when p(x)=0, log2(0) produces NaN)
    EntropyX[is.na(EntropyX)] <- 0

    # Return entropy matrix rounded to 6 decimal places
    return(round(EntropyX, 6))
}

## Internal Function for calculating gene pair joint entropy matrix
## Calculate joint entropy for two genes
get_entropyXY_matrix <- function(Xt, g1_site, g2_site, cellTot, binsNum) {
    # Initialize joint count matrix
    joint_counts <- matrix(0, nrow = binsNum, ncol = binsNum)

    # Calculate joint frequency distribution of discretized values for two genes
    for (i in 1:cellTot) {
      joint_counts[Xt[g1_site, i], Xt[g2_site, i]] <- joint_counts[Xt[g1_site, i], Xt[g2_site, i]] + 1
    }

    # Calculate joint probabilities
    prob <- joint_counts / cellTot

    # Calculate joint entropy
    entropies <- prob * log2(prob)

    # Handle NaN values
    entropies[is.na(entropies)] <- 0

    # Return joint entropy matrix rounded to 6 decimal places
    return(round(entropies, 6))
}

## Internal Function for discretizing gene expression data
## Use a custom discretization function to handle unique breaks
discretize_data <- function(mat, bins) {
    result <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    for (i in 1:nrow(mat)) {
      # Handle case where all values are the same
      if (diff(range(mat[i, ])) == 0) {
        result[i, ] <- 1
      } else {
         # Create unique breaks
        unique_breaks <- seq(min(mat[i, ]), max(mat[i, ]), length.out = bins + 1)
        # Ensure breaks are unique by adding a tiny epsilon if needed
        if (any(duplicated(unique_breaks))) {
         unique_breaks <- unique_breaks + seq(0, 1e-10, length.out = length(unique_breaks))
         }
         result[i, ] <- as.numeric(cut(mat[i, ], breaks = unique_breaks, labels = FALSE, include.lowest = TRUE))
       }
    }
    return(result)
}
