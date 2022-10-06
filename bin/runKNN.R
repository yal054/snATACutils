runKNN <- function(obj, eigs.dims, weight.by.lambda, k, nn.eps, save.knn,
                    filename, snn, snn.prune, method) {
   UseMethod("runKNN", obj);
 }

#' @export
runKNN.default <- function(
  obj,
  eigs.dims,
  weight.by.lambda = FALSE,
  k = 15,
  nn.eps = 0,
   save.knn = FALSE,
   filename = NULL,
   snn = FALSE,
  snn.prune = 1/15,
   method = "RANN"
 ){

     cat("Epoch: checking input parameters\n", file = stderr())
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "snap")){
            stop("obj is not a snap obj")
        }
    }
    
    if(!(isDimReductComplete(obj@smat))){
        stop("dimentionality reduction is not complete, run 'runDimReduct' first")
    }
    
    ncell = nrow(obj);
    nvar = dimReductDim(obj@smat);
    
    if(missing(eigs.dims)){
        stop("eigs.dims is missing")
    }else{
        if(is.null(eigs.dims)){
            eigs.dims=1:nvar;    
        }else{
            if(any(eigs.dims > nvar) ){
                stop("'eigs.dims' exceeds PCA dimentions number");
            }        
        }
    }
    
    if(save.knn){
        if(is.null(filename)){
            stop("save.knn is TRUE but filename is NULL")
        }else{            
            if(!file.create(filename)){
                stop("fail to create filename")
            }            
        }
    }
    
    if(!is.logical(weight.by.lambda)){
        stop("weight.by.lambda must be a logical variable")
    }
    
    data.use = weightDimReduct(obj@smat, eigs.dims, weight.by.lambda);
    
    if (ncell < k) {
      warning("k set larger than number of cells. Setting k to number of cells - 1.")
      k <- ncell - 1
    }
    
    if(is.na(as.integer(k))){
        stop("k must be an integer")
    }else{
        if(k < 10 || k > 50){
            warning("too small or too large k, recommend to set k within range [10 - 50]")
        }
    }
    
     cat("Epoch: computing nearest neighbor graph\n", file = stderr())

     # exclude self neibours
   if (method == "RANN") {
     message("Using RANN to get the KNN result.")
     nn.res <- nn2(
      data = data.use,
       k = k,
       searchtype = 'standard',
       eps = nn.eps)
     nn.ranked <- nn.res$nn.idx
   } else {
     message("Using Annoy to get the KNN result.")
     nn.res <- AnnoyNN(data = data.use, k = k)
     nn.ranked <- nn.res$nn.idx
   }

     j <- as.numeric(x = t(x = nn.ranked))
     i <- ((1:length(x = j)) - 1) %/% k + 1    
    edgeList <- data.frame(i, j, 1)

     if(snn){
         cat("Epoch: converting knn graph into snn graph\n", file = stderr())    
        g = graph_from_edgelist(as.matrix(edgeList[,c(1,2)]), directed=FALSE);
        adj = as(similarity(g), "sparseMatrix");
        i = adj@i+1;
        j = findInterval(seq(adj@x)-1,adj@p[-1])+1;
        w = adj@x;
        idx = which(w >= snn.prune);
        edgeList = data.frame(i[idx], j[idx], w[idx]);
    }
    
    if(save.knn){
        cat("Epoch: writing resulting graph into a file\n", file = stderr())
        writeEdgeListToFile(edgeList, filename);
        obj@graph = newKgraph(file=filename, k=k, snn=snn, snn.prune=snn.prune);
    }else{
        kmat = Matrix(0, ncell, ncell, sparse=TRUE);
        kmat = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3]);
        obj@graph = newKgraph(mat=kmat, k=k, snn=snn, snn.prune=snn.prune);
    }
    gc();
    return(obj);
} 

####################
 ## From Seurat
 ####################
 # Run annoy
 #
 # @param data Data to build the index with
 # @param query A set of data to be queried against data
 # @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
 # "hamming"
 # @param n.trees More trees gives higher precision when querying
 # @param k Number of neighbors
 # @param search.k During the query it will inspect up to search_k nodes which
 # gives you a run-time tradeoff between better accuracy and speed.
 # @param include.distance Include the corresponding distances
 # @param index optional index object, will be recomputed if not provided
 #
 AnnoyNN <- function(data,
                     query = data,
                     metric = "euclidean",
                     n.trees = 50,
                     k,
                     search.k = -1,
                     include.distance = TRUE,
                     index = NULL
 ) {
   idx <- index %||% AnnoyBuildIndex(
     data = data,
     metric = metric,
     n.trees = n.trees)
   nn <- AnnoySearch(
     index = idx,
     query = query,
     k = k,
     search.k = search.k,
     include.distance = include.distance)
   nn$idx <- idx
   nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
   return(nn)
 }

 # Build the annoy index
 #
 # @param data Data to build the index with
 # @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
 # "hamming"
 # @param n.trees More trees gives higher precision when querying
 #
 #' @importFrom RcppAnnoy AnnoyEuclidean AnnoyAngular AnnoyManhattan AnnoyHamming
 #
 AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
   f <- ncol(x = data)
   a <- switch(
     EXPR = metric,
     "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
     "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
     "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
     "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
     stop ("Invalid metric")
   )
   for (ii in seq(nrow(x = data))) {
     a$addItem(ii - 1, data[ii, ])
   }
   a$build(n.trees)
   return(a)
 }

 # Search an Annoy approximate nearest neighbor index
 #
 # @param Annoy index, built with AnnoyBuildIndex
 # @param query A set of data to be queried against the index
 # @param k Number of neighbors
 # @param search.k During the query it will inspect up to search_k nodes which
 # gives you a run-time tradeoff between better accuracy and speed.
 # @param include.distance Include the corresponding distances in the result
 #
 # @return A list with 'nn.idx' (for each element in 'query', the index of the
 # nearest k elements in the index) and 'nn.dists' (the distances of the nearest
 # k elements)
 #
 #' @importFrom future plan
 #' @importFrom future.apply future_lapply
 #
 AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
   n <- nrow(x = query)
   idx <- matrix(nrow = n,  ncol = k)
   dist <- matrix(nrow = n, ncol = k)
   convert <- methods::is(index, "Rcpp_AnnoyAngular")
   if (!inherits(x = future::plan(), what = "multicore")) {
     oplan <- future::plan(strategy = "sequential")
     on.exit(future::plan(oplan), add = TRUE)
   }
   res <- future.apply::future_lapply(X = 1:n, FUN = function(x) {
     res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
     # Convert from Angular to Cosine distance
     if (convert) {
       res$dist <- 0.5 * (res$dist * res$dist)
     }
     list(res$item + 1, res$distance)
   })
   for (i in 1:n) {
     idx[i, ] <- res[[i]][[1]]
     if (include.distance) {
       dist[i, ] <- res[[i]][[2]]
     }
   }
   return(list(nn.idx = idx, nn.dists = dist))
 }

 `%||%` <- function(lhs, rhs) {
   if (!is.null(x = lhs)) {
     return(lhs)
   } else {
     return(rhs)
   }
 }

isDimReductComplete <- function(obj) {
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "dim.reduct")){
            stop("obj is not a dim.reduct")
        }
    }
    
    if(!((nrow(obj@dmat) > 0) && (length(obj@sdev) > 0))){
        return(FALSE);
    }    
    
    if(ncol(obj@dmat) != length(obj@sdev)){
        return(FALSE);
    }
    
    return(TRUE);
}


weightDimReduct <- function(obj, pca.dims, weight.by.sd=TRUE){
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "dim.reduct")){
            stop("obj is not a dim.reduct")
        }
    }
    
    if(weight.by.sd){
        data.use = obj@dmat[,pca.dims] %*% diag(obj@sdev[pca.dims]) ;
    }else{
        data.use = obj@dmat[,pca.dims];
    }
    return(data.use);
}


#' @importFrom methods is
dimReductDim <- function(obj){
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "dim.reduct")){
            stop("obj is not a dim.reduct")
        }
    }
    return(length(obj@sdev));
}


newKgraph <- function (mat=NULL, file=NULL, k=NULL, snn=NULL, snn.prune=NULL) {
    if(is.null(mat)){
        mat = Matrix::Matrix(0,0,0, sparse=TRUE);
    }
    if(is.null(file)){
        file = character();
    }

    if(is.null(k)){
        k = numeric();
    }

    if(is.null(snn)){
        snn = FALSE;
    }

    if(is.null(snn.prune)){
        snn.prune = numeric();
    }
    
    res = new("kgraph", 
              mat=mat,
              file=file,
              k=k,
              snn=snn,
              snn.prune=snn.prune
              );
    return(res)    
}
