## Calulate z-score value for the SSN method
## Internal function
ssn_zscore <- function(deta, pcc, nn) {

  if (any(is.na(c(deta, pcc, nn)))) return(NA)
  if (nn <= 1) stop("nn should be larger than 1")

  if (abs(pcc - 1) < 1e-10) {
    pcc <- 0.9999999999
  } else if (abs(pcc + 1) < 1e-10) {
    pcc <- -0.9999999999
  }

  denominator <- (1 - pcc^2) / (nn - 1)
  if (abs(denominator) < 1e-10) {
    warning("The denominator is close to zero, and the results are not stable.")
    return(NA)
  }

  z <- deta / denominator
  return(z)
}

## Calculate ssn score using paired normal and tumor data for the PairedSSN method
## Internal function
SSN_paired <- function(reg_normal_mat, tar_normal_mat, reg_sample, tar_sample, cormethod = "pearson", nn) {

  R0 <- cor(x = t(reg_normal_mat), y = t(tar_normal_mat), method = cormethod)

  inter1 <- intersect(rownames(reg_normal_mat), rownames(reg_sample))
  reg_normal_mat <- reg_normal_mat[inter1, ]
  reg_sample <- reg_sample[inter1, ]

  inter2 <- intersect(rownames(tar_normal_mat), rownames(tar_sample))
  tar_normal_mat <- tar_normal_mat[inter2, ]
  tar_sample <- tar_sample[inter2, ]

  R1 <- cor(x = t(cbind(reg_normal_mat, reg_sample)), y = t(cbind(tar_normal_mat, tar_sample)), method = cormethod)
  detal <- R1 - R0
  Z <- detal / ((1 - R0^2) / (nn - 1))
  p <- 1 - pnorm(abs(Z))
  ssn_result <- list(detal = detal, p = p)

  return(ssn_result)
}

## Functions to perform Random Walk with Restart on Multiplex Networks (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
simplify.layers <- function(Input_Layer){

  ## Undirected Graphs
  Layer <- as_undirected(Input_Layer, mode = c("collapse"),
                         edge.attr.comb = igraph_opt("edge.attr.comb"))

  ## Unweighted or Weigthed Graphs
  if (is_weighted(Layer)){
    b <- 1
    weigths_layer <- E(Layer)$weight
    if (min(weigths_layer) != max(weigths_layer)){
      a <- min(weigths_layer)/max(weigths_layer)
      range01 <- (b-a)*(weigths_layer-min(weigths_layer))/
        (max(weigths_layer)-min(weigths_layer)) + a
      E(Layer)$weight <- range01
    } else {
      E(Layer)$weight <- rep(1, length(weigths_layer))
    }
  } else {
    E(Layer)$weight <- rep(1, ecount(Layer))
  }

  ## Simple Graphs
  Layer <-
    igraph::simplify(Layer,remove.multiple = TRUE,remove.loops = TRUE,
                     edge.attr.comb=mean)

  return(Layer)
}

## Add missing nodes in some of the layers (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {

  add_vertices(Layers,
               length(NodeNames[which(!NodeNames %in% V(Layers)$name)]),
               name=NodeNames[which(!NodeNames %in%  V(Layers)$name)])
}

## Is this R object a Multiplex object? (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
isMultiplex <- function (x)
{
  is(x,"Multiplex")
}

## Scores for the seeds of the multiplex network (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
get.seed.scoresMultiplex <- function(Seeds,Number_Layers,tau) {

  Nr_Seeds <- length(Seeds)

  Seeds_Seeds_Scores <- rep(tau/Nr_Seeds,Nr_Seeds)
  Seed_Seeds_Layer_Labeled <-
    paste0(rep(Seeds,Number_Layers),sep="_",rep(seq(Number_Layers),
                                                length.out = Nr_Seeds*Number_Layers,each=Nr_Seeds))

  Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                            Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)

  return(Seeds_Score)
}

## Geometric mean: R computation (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
geometric.mean <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
  }

  return(FinalScore)
}

## Create Multiplex and Multiplex-Heterogeneous objects (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
create.multiplex <- function(LayersList,...){

  if (!class(LayersList) == "list"){
    stop("The input object should be a list of graphs.")
  }


  Number_of_Layers <- length(LayersList)
  SeqLayers <- seq(Number_of_Layers)
  Layers_Name <- names(LayersList)

  if (!all(sapply(SeqLayers, function(x) is_igraph(LayersList[[x]])))){
    stop("Not igraph objects")
  }

  Layer_List <- lapply(SeqLayers, function (x) {
    if (is.null(V(LayersList[[x]])$name)){
      LayersList[[x]] <-
        set_vertex_attr(LayersList[[x]],"name",
                        value=seq(1,vcount(LayersList[[x]]),by=1))
    } else {
      LayersList[[x]]
    }
  })

  ## We simplify the layers
  Layer_List <-
    lapply(SeqLayers, function(x) simplify.layers(Layer_List[[x]]))

  ## We set the names of the layers.

  if (is.null(Layers_Name)){
    names(Layer_List) <- paste0("Layer_", SeqLayers)
  } else {
    names(Layer_List) <- Layers_Name
  }

  ## We get a pool of nodes (Nodes in any of the layers.)
  Pool_of_Nodes <-
    sort(unique(unlist(lapply(SeqLayers,
                              function(x) V(Layer_List[[x]])$name))))

  Number_of_Nodes <- length(Pool_of_Nodes)

  Layer_List <-
    lapply(Layer_List, add.missing.nodes,Number_of_Layers,Pool_of_Nodes)

  # We set the attributes of the layer
  counter <- 0
  Layer_List <- lapply(Layer_List, function(x) {
    counter <<- counter + 1;
    set_edge_attr(x,"type",E(x), value = names(Layer_List)[counter])
  })


  MultiplexObject <- c(Layer_List,list(Pool_of_Nodes=Pool_of_Nodes,
                                       Number_of_Nodes_Multiplex=Number_of_Nodes,
                                       Number_of_Layers=Number_of_Layers))

  class(MultiplexObject) <- "Multiplex"

  return(MultiplexObject)
}

## Compute the matrices and perform the Random Walks (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
compute.adjacency.matrix <- function(x,delta = 0.5)
{
  if (!isMultiplex(x)) {
    stop("Not a Multiplex object")
  }
  if (delta > 1 || delta <= 0) {
    stop("Delta should be between 0 and 1")
  }

  N <- x$Number_of_Nodes_Multiplex
  L <- x$Number_of_Layers

  ## We impose delta=0 in the monoplex case.
  if (L==1){
    delta = 0
  }

  Layers_Names <- names(x)[seq(L)]

  ## IDEM_MATRIX.
  Idem_Matrix <- Matrix::Diagonal(N, x = 1)

  counter <- 0
  Layers_List <- lapply(x[Layers_Names],function(x){

    counter <<- counter + 1;
    if (is_weighted(x)){
      Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE,
                                              attr = "weight")
    } else {
      Adjacency_Layer <-  as_adjacency_matrix(x,sparse = TRUE)
    }

    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),
                                       order(colnames(Adjacency_Layer))]
    colnames(Adjacency_Layer) <-
      paste0(colnames(Adjacency_Layer),"_",counter)
    rownames(Adjacency_Layer) <-
      paste0(rownames(Adjacency_Layer),"_",counter)
    Adjacency_Layer
  })

  MyColNames <- unlist(lapply(Layers_List, function (x) unlist(colnames(x))))
  MyRowNames <- unlist(lapply(Layers_List, function (x) unlist(rownames(x))))
  names(MyColNames) <- c()
  names(MyRowNames) <- c()
  SupraAdjacencyMatrix <- (1-delta)*(bdiag(unlist(Layers_List)))
  colnames(SupraAdjacencyMatrix) <-MyColNames
  rownames(SupraAdjacencyMatrix) <-MyRowNames

  offdiag <- (delta/(L-1))*Idem_Matrix

  i <- seq_len(L)
  Position_ini_row <- 1 + (i-1)*N
  Position_end_row <- N + (i-1)*N
  j <- seq_len(L)
  Position_ini_col <- 1 + (j-1)*N
  Position_end_col <- N + (j-1)*N

  for (i in seq_len(L)){
    for (j in seq_len(L)){
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                             (Position_ini_col[j]:Position_end_col[j])] <- offdiag
      }
    }
  }

  SupraAdjacencyMatrix <- methods::as(SupraAdjacencyMatrix, "dgCMatrix")
  return(SupraAdjacencyMatrix)
}

## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
regular.mean <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- mean(Scores[seq(from = i, to = N*L, by=N)])
  }

  return(FinalScore)
}

## (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
sumValues <- function(Scores, L, N) {

  FinalScore <- numeric(length = N)

  for (i in seq_len(N)){
    FinalScore[i] <- sum(Scores[seq(from = i, to = N*L, by=N)])
  }

  return(FinalScore)
}

## Computes column normalization of an adjacency matrix (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
normalize.multiplex.adjacency <- function(x)
{
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }

  Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                             sparseResult = FALSE)))

  return(Adj_Matrix_Norm)
}

## Performs Random Walk with Restart on a Multiplex Network (https://github.com/alberto-valdeolivas/RandomWalkRestartMH/)
## Internal function
Random.Walk.Restart.Multiplex <- function(x, MultiplexObject, Seeds,
  r=0.7,tau,MeanType="Geometric", DispResults="TopScores",...){

  ### We control the different values.
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }

  if (!isMultiplex(MultiplexObject)) {
    stop("Not a Multiplex object")
  }

  L <- MultiplexObject$Number_of_Layers
  N <- MultiplexObject$Number_of_Nodes

  Seeds <- as.character(Seeds)
  if (length(Seeds) < 1 | length(Seeds) >= N){
    stop("The length of the vector containing the seed nodes is not
         correct")
  } else {
    if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
      stop("Some of the seeds are not nodes of the network")

    }
  }

  if (r >= 1 || r <= 0) {
    stop("Restart partameter should be between 0 and 1")
  }

  if(missing(tau)){
    tau <- rep(1,L)/L
  } else {
    tau <- as.numeric(tau)
    if (sum(tau)/L != 1) {
      stop("The sum of the components of tau divided by the number of
           layers should be 1")
    }
    }

  if(!(MeanType %in% c("Geometric","Arithmetic","Sum"))){
    stop("The type mean should be Geometric, Arithmetic or Sum")
  }

  if(!(DispResults %in% c("TopScores","Alphabetic"))){
    stop("The way to display RWRM results should be TopScores or
         Alphabetic")
  }

  ## We define the threshold and the number maximum of iterations for
  ## the random walker.
  Threeshold <- 1e-10
  NetworkSize <- ncol(x)

  ## We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1

  ## We compute the scores for the different seeds.
  Seeds_Score <- get.seed.scoresMultiplex(Seeds,L,tau)

  ## We define the prox_vector(The vector we will move after the first RWR
  ## iteration. We start from The seed. We have to take in account
  ## that the walker with restart in some of the Seed nodes, depending on
  ## the score we gave in that file).
  prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)

  prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])

  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector

  while(residue >= Threeshold){

    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1;
  }

  NodeNames <- character(length = N)
  Score = numeric(length = N)

  rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
  rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])

  if (MeanType=="Geometric"){
    rank_global$Score <- geometric.mean(as.vector(prox_vector[,1]),L,N)
  } else {
    if (MeanType=="Arithmetic") {
      rank_global$Score <- regular.mean(as.vector(prox_vector[,1]),L,N)
    } else {
      rank_global$Score <- sumValues(as.vector(prox_vector[,1]),L,N)
    }
  }

  if (DispResults=="TopScores"){
    ## We sort the nodes according to their score.
    Global_results <-
      rank_global[with(rank_global, order(-Score, NodeNames)), ]

    ### We remove the seed nodes from the Ranking and we write the results.
    Global_results <-
      Global_results[which(!Global_results$NodeNames %in% Seeds),]
  } else {
    Global_results <- rank_global
  }

  rownames(Global_results) <- c()

  RWRM_ranking <- list(RWRM_Results = Global_results,Seed_Nodes = Seeds)

  class(RWRM_ranking) <- "RWRM_Results"
  return(RWRM_ranking)
  }

## Perform Random Walk with Restart on a multiplex network
## Internal function
RWR_MN <- function(AdjMatrix_edge0 = AdjMatrix_edge0, edge0_MultiplexObject = edge0_MultiplexObject, SNP_name = SNP_name, scores = 0.0001) {

  AdjMatrixNorm_edge0 <- normalize.multiplex.adjacency(AdjMatrix_edge0)
  SeedGenes <- edge0_MultiplexObject[["Pool_of_Nodes"]][na.omit(match(SNP_name, edge0_MultiplexObject[["Pool_of_Nodes"]]))]
  RWR_edge0_Results_top <- list()
  for (i in 1:length(SeedGenes)) {
    RWR_edge0_Results_tmp <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_edge0, edge0_MultiplexObject, SeedGenes[i], r = 0.6)
    RWR_edge0_Results_tmp1 <- data.frame("SeedGenes" = RWR_edge0_Results_tmp[["Seed_Nodes"]], RWR_edge0_Results_tmp[[1]])
    RWR_edge0_Results_top[[i]] <- RWR_edge0_Results_tmp1[RWR_edge0_Results_tmp1[, 3] > scores, ]
  }
  RWR_edge0_Results_top_all <- Reduce(rbind, RWR_edge0_Results_top)

  return(RWR_edge0_Results_top_all)
}


## Construct random network with the same degree as the original one
## Internal function
dir_generate_srand <- function(s1, ntry) {

  # s1 - the adjacency matrix of a directed network
  # ntry - (optional) the number of rewiring steps. If none is given ntry=4*(# of edges in the network)
  # Output: dir_srand - the adjacency matrix of a randomized network with the same set of in- and out-degrees as the original one
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
    h <- round(boxsize / 2 * sum(sign(s1)) + pracma::eps(1))
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
  res <- (colSums(B[[1]] & B[[2]]) * n - a[1, ] * a[2, ]) / sqrt(a[1, ] * a[2, ] * (n - a[1, ]) * (n - a[2, ]) / (n - 1) + pracma::eps(1))

  return(res)
}

## Internal function cluster from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster <- function(graph, method="MCL", expansion = 2, inflation = 2,
                  hcmethod = "average", directed = FALSE, outfile = NULL, ...)
{
  #method<-match.arg(method)
  if(method=="FN"){
	  graph <- simplify(graph)
	  fc <- fastgreedy.community(graph, merges = TRUE, modularity = TRUE)
	  membership <- membership(fc)
	  if(!is.null(V(graph)$name)){
              names(membership) <- V(graph)$name
	  }
	  if(!is.null(outfile)){
	      cluster.save(cbind(names(membership),membership),outfile=outfile)
	  }else{
	      return(membership)
	  }
  }else if(method=="LINKCOMM"){
	  edgelist <- get.edgelist(graph)
	  if(!is.null(E(graph)$weight)){
              edgelist <- cbind(edgelist,E(graph)$weight)
	  }
	  lc <- getLinkCommunities(edgelist,plot=FALSE,directed=directed,hcmethod=hcmethod)
	  if(!is.null(outfile)){
		  cluster.save(lc$nodeclusters,outfile=outfile)
	  }else{
		  return(lc$nodeclusters)
	  }
  }else if(method=="MCL"){
        adj <- matrix(rep(0,length(V(graph))^2),nrow=length(V(graph)),ncol=length(V(graph)))
        for(i in seq_along(V(graph))){
            neighbors <- neighbors(graph,v=V(graph)$name[i],mode="all")
            j <- match(neighbors$name,V(graph)$name,nomatch=0)
            adj[i,j] = 1
        }
        lc <- mcl(adj,addLoops=TRUE,expansion=expansion,inflation=inflation,allow1=TRUE,max.iter=100,ESM=FALSE)
        lc$name <- V(graph)$name
        lc$Cluster <- lc$Cluster

        if(!is.null(outfile)){
            cluster.save(cbind(lc$name,lc$Cluster),outfile=outfile)
        }else{
        result <- lc$Cluster
        names(result) <- V(graph)$name
        return(result)
        }
  }else if(method=="MCODE"){
	  compx <- mcode(graph,vwp=0.9,haircut=T,fluff=T,fdt=0.1)
	  index <- which(!is.na(compx$score))
	  membership <- rep(0,vcount(graph))
	  for(i in seq_along(index)){
	      membership[compx$COMPLEX[[index[i]]]]<-i
	  }
	      if(!is.null(V(graph)$name)) names(membership)<-V(graph)$name
	  if(!is.null(outfile)){
		  cluster.save(cbind(names(membership),membership),outfile=outfile)
		  invisible(NULL)
	  }else{
		  return(membership)
	  }
  }
}

## Internal function cluster.save from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster.save <- function(membership, outfile){
	wd <- dirname(outfile)
	wd <- ifelse(wd==".",paste(wd,"/",sep=""),wd)
	filename <- basename(outfile)
	if((filename=="")||(grepl(":",filename))){
		filename <- "membership.txt"
	}else if(grepl("\\.",filename)){
		filename <- sub("\\.(?:.*)",".txt", filename)
	}
	write.table(membership,file=paste(wd,filename,sep="/"),
              row.names=FALSE,col.names=c("node","cluster"),quote =FALSE)
}

## Internal function mcode.vertex.weighting from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.vertex.weighting<-function(graph, neighbors){
	stopifnot(is_igraph(graph))
	weight <- lapply(seq_len(vcount(graph)),
                 function(i){
		              subg<-induced.subgraph(graph,neighbors[[i]])
		              core<-graph.coreness(subg)
		              k<-max(core)
				          ### k-coreness
				          kcore<-induced.subgraph(subg,which(core==k))
				          if(vcount(kcore)>1){
					          if(any(is.loop(kcore))){
						          k*ecount(kcore)/choose(vcount(kcore)+1,2)
					          }else{
						          k*ecount(kcore)/choose(vcount(kcore),2)
					          }
				          }else{
                                             0
				          }
				         }
			 )

	return(unlist(weight))
}

## Internal function mcode.find.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complex <- function(neighbors, neighbors.indx, vertex.weight,
                             vwp, seed.vertex, seen)
{

    res<-.C("complex",as.integer(neighbors),as.integer(neighbors.indx),
          as.single(vertex.weight),as.single(vwp),as.integer(seed.vertex),
          seen=as.integer(seen),COMPLEX=as.integer(rep(0,length(seen))), PACKAGE = "S2N"
          )

	  return(list(seen=res$seen,COMPLEX=which(res$COMPLEX!=0)))
}

## Internal function mcode.find.complexex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complexex <- function(graph, neighbors, vertex.weight, vwp)
{
	seen<-rep(0,vcount(graph))

	neighbors<-lapply(neighbors,function(item){item[-1]})
	neighbors.indx<-cumsum(unlist(lapply(neighbors,length)))

	neighbors.indx<-c(0,neighbors.indx)
	neighbors<-unlist(neighbors)-1

	COMPLEX<-list()
	n<-1
        w.order<-order(vertex.weight,decreasing=TRUE)
	for(i in w.order){
		if(!(seen[i])){
			res<-mcode.find.complex(neighbors,neighbors.indx,vertex.weight,vwp,i-1,seen)
			if(length(res$COMPLEX)>1){
				COMPLEX[[n]]<-res$COMPLEX
				seen<-res$seen
				n<-n+1
			}
		}
	}
	rm(neighbors)
	return(list(COMPLEX=COMPLEX,seen=seen))
}

## Internal function mcode.fluff.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.fluff.complex <- function(graph, vertex.weight, fdt=0.8, complex.g, seen)
{
	seq_complex.g<-seq_along(complex.g)
	for(i in seq_complex.g){
	    node.neighbor<-unlist(neighborhood(graph,1,complex.g[i]))
	    if(length(node.neighbor)>1){
                subg<-induced.subgraph(graph,node.neighbor)
            if(graph.density(subg, loops=FALSE)>fdt){
                complex.g<-c(complex.g,node.neighbor)
       }
		 }
	}

	return(unique(complex.g))
}

## Internal function mcode.post.process from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.post.process<-function(graph, vertex.weight, haircut, fluff, fdt=0.8,
                             set.complex.g, seen)
{
	indx<-unlist(lapply(set.complex.g,
                      function(complex.g){
		          if(length(complex.g)<=2)
			      0
			  else
		              1
			  }
		    ))
	set.complex.g<-set.complex.g[indx!=0]
	set.complex.g<-lapply(set.complex.g,
                        function(complex.g){
			    coreness<-graph.coreness(induced.subgraph(graph,complex.g))
			    if(fluff){
				complex.g<-mcode.fluff.complex(graph,vertex.weight,fdt,complex.g,seen)
			    if(haircut){
				## coreness needs to be recalculated
				coreness<-graph.coreness(induced.subgraph(graph,complex.g))
				complex.g<-complex.g[coreness>1]
				}
				}else if(haircut){
				complex.g<-complex.g[coreness>1]
				}
				return(complex.g)
				})
	set.complex.g<-set.complex.g[lapply(set.complex.g,length)>2]
	return(set.complex.g)
}

## Internal function mcode from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode <- function(graph, vwp=0.5, haircut=FALSE, fluff=FALSE, fdt=0.8, loops=TRUE)
{
	stopifnot(is_igraph(graph))
	if(vwp>1 | vwp <0){
            stop("vwp must be between 0 and 1")
	}
	if(!loops){
            graph<-simplify(graph,remove.multiple=FALSE,remove.loops=TRUE)
	}
	neighbors<-neighborhood(graph,1)
	W<-mcode.vertex.weighting(graph,neighbors)
	res<-mcode.find.complexex(graph,neighbors=neighbors,vertex.weight=W,vwp=vwp)
	COMPLEX<-mcode.post.process(graph,vertex.weight=W,haircut=haircut,fluff=fluff,
                              fdt=fdt,res$COMPLEX,res$seen)
	score<-unlist(lapply(COMPLEX,
                       function(complex.g){
		           complex.g<-induced.subgraph(graph,complex.g)
			   if(any(is.loop(complex.g)))
			   score<-ecount(complex.g)/choose(vcount(complex.g)+1,2)*vcount(complex.g)
			   else
			   score<-ecount(complex.g)/choose(vcount(complex.g),2)*vcount(complex.g)
			   return(score)
			}
		    ))
	order_score<-order(score,decreasing=TRUE)
	return(list(COMPLEX=COMPLEX[order_score],score=score[order_score]))
}

## Internal function csn_edge from CSmiR methods
## (https://github.com/zhangjunpeng411/CSmiR) with GPL-2 license.
CSN_net_bootstrap <- function(reg_cancer, tar_cancer, boxsize = 0.1, bootstrap_betw_point = 5, bootstrap_num = 100, p.value.cutoff = 0.05) {

  regs_num <- ncol(reg_cancer)
  tars_num <- ncol(tar_cancer)
  cell_num <- nrow(reg_cancer)
  bootstrap_sample <- lapply(seq(bootstrap_num), function(i) sample(seq(cell_num), bootstrap_betw_point * (cell_num - 1), replace = TRUE))
  reg_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(reg_cancer, reg_cancer[bootstrap_sample[[i]], ]))
  tar_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(tar_cancer, tar_cancer[bootstrap_sample[[i]], ]))
  res <- matrix(NA, nrow = regs_num * tars_num, ncol = cell_num + 2)
  for (i in seq(regs_num)) {
    for (j in seq(tars_num)) {
      res[(i - 1) * tars_num + j, 1] <- colnames(reg_cancer)[i]
      res[(i - 1) * tars_num + j, 2] <- colnames(tar_cancer)[j]
      res[(i - 1) * tars_num + j, 3:(cell_num + 2)] <- do.call(WGCNA::pmedian, lapply(
        seq(bootstrap_num),
        function(k) {
          csn_edge(reg_bootstrap[[k]][, i],
                   tar_bootstrap[[k]][, j],
                   boxsize = boxsize
          )[seq(cell_num)]
        }
      ))
    }
  }

  q <- -qnorm(p.value.cutoff)
  res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i + 2]) > q), seq(2)])

  return(res_list)
}

## Internal function csn_edge from CSmiR methods
## (https://github.com/zhangjunpeng411/CSmiR) with GPL-2 license.
CSN_net <- function(reg_cancer, tar_cancer, boxsize = 0.1, p.value.cutoff = 0.05, num.cores = 2) {

  regs_num <- ncol(reg_cancer)
  tars_num <- ncol(tar_cancer)
  cell_num <- nrow(reg_cancer)
  int_num <- regs_num * tars_num
  res <- matrix(NA, nrow = int_num, ncol = cell_num + 2)
  index <- matrix(NA, nrow = int_num, ncol = 2)

  for (i in seq(regs_num)){
    for (j in seq(tars_num)){
      index[(i-1)*tars_num+j, 1] <- i
      index[(i-1)*tars_num+j, 2] <- j
    }
  }

  # get number of cores to run
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  interin <- foreach(k = seq(int_num),
                     .packages = c("pracma"),
                     .export = c("csn_edge")) %dopar% {
                       csn_edge(reg_cancer[, index[k, 1]], tar_cancer[, index[k, 2]], boxsize = boxsize)
                     }
  # shut down the workers
  stopCluster(cl)
  stopImplicitCluster()

  interin <- do.call(rbind, interin)

  for (i in seq(int_num)){
    res[i, 1] <- colnames(reg_cancer)[index[i, 1]]
    res[i, 2] <- colnames(tar_cancer)[index[i, 2]]
  }

  res[, 3:(cell_num + 2)] <- interin

  q <- -qnorm(p.value.cutoff)
  res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i+2]) > q), seq(2)])

  return(res_list)
}

## Internal Function of network inference using Pearson (https://github.com/zhangjunpeng411/Scan)
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Pearson <- function(miRExp, mRExp, p.value = 0.05){

  cor.pvalue <- corAndPvalue(t(mRExp), t(miRExp), method = "pearson")$p
  adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

  return(adj.Matrix)
}

## Internal Function of network inference using Spearman (https://github.com/zhangjunpeng411/Scan)
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Spearman <- function(miRExp, mRExp, p.value = 0.05){

  cor.pvalue <- corAndPvalue(t(mRExp), t(miRExp), method = "spearman")$p
  adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

  return(adj.Matrix)
}

## Internal Function of network inference using Kendall (https://github.com/zhangjunpeng411/Scan)
# miRExp: Gene expression values of miRNAs
# mRExp: Gene expression values of mRNAs
# p.value: A cutoff of p-value
# Output: adj.Matrix is a zero-one matrix between miRNAs and mRNAs
Kendall <- function(miRExp, mRExp, p.value = 0.05){

  cor.pvalue <- corAndPvalue(t(mRExp), t(miRExp), method = "kendall")$p
  adj.Matrix <- ifelse(cor.pvalue < p.value, 1, 0)

  return(adj.Matrix)
}
