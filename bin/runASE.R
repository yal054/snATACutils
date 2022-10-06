runEigDecomp2 <- function(obj, num.eigs){
    nmat = obj@jmat@nmat;
    diag(nmat) = 0;
    norm_p1 = nmat
    d_norm1 <- Matrix::rowSums(norm_p1);
    d_rot1 <- Diagonal(x = d_norm1 ^ -.5);
    transitions <- as.matrix(d_rot1 %*% norm_p1 %*% d_rot1);
    diag(transitions) = 0
    message("-- using Adjacency spectral embedding...");
    eig_transitions = eig_decomp2(transitions, num.eigs)
    obj@smat@dmat = eig_transitions$vectors[,2:(num.eigs+1)];
    obj@smat@sdev = eig_transitions$value[2:(num.eigs+1)];
    return(obj)
}


#' Eig decomposition.
 #' @param M Matrix
 #' @param n_eighs integer
 #' @param sym bool, is M symmetric, default is isSymmetric(M)
 #' @param method character, "arpack" or "RSpectra" for decomposition, default is arpack
 #' 
 #' @export
eig_decomp2 <- function(M, n_eigs, sym = isSymmetric(M)) {
    n <- nrow(M)
    f <- function(x, A = NULL) as.matrix(A %*% x)
    #wh <- if (sym) 'LA' else 'LM'
    wh <- 'LM'
    #constraints: n >= ncv > nev
    message("Use RSpectra::eigs for eig decomposition.")
    ar <- RSpectra::eigs(A = M, k = n_eigs + 1, which = wh)
    if (!sym) {
        ar$vectors <- Re(ar$vectors)
        ar$values  <- Re(ar$values)
    }
    if (length(dim(ar$vectors)) == 0L) {
         ar$vectors <- matrix(ar$vectors, ncol = 1L)
    }
    return(ar)
}

#' Dimentionality Reduction by Adjacency spectral embedding
#'
#' This function takes a snap obj as input and runs adjacency spectral embedding
#' for dimentionality reduction. 
#' 
#' @param obj A snap obj
#' @param input.mat Input matrix c("bmat", "pmat").
#' @param num.eigs Number of eigenvectors to be computed [20].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runDiffusionMaps(
#'    obj=demo.sp, 
#'    input.mat="bmat", 
#'    num.eigs=20
#' );
#' @import Matrix
#' @export
runAdjSpectral <- function(
    obj,
    input.mat=c("bmat", "pmat"), 
    num.eigs=20
){
    nmat.outlier = 0.999
    message("Epoch: checking the inputs ...");
    if(missing(obj)){
        stop("obj is missing");
    }else{
        if(!is(obj, "snap")){
            stop("obj is not a snap obj")
        }        
    }
    
    # 2. check if input matrix exists
    input.mat = match.arg(input.mat);    
    if(input.mat == "bmat"){
        data.use = obj@bmat;
        peak.use = obj@feature;
    }else if(input.mat == "pmat"){
        data.use = obj@pmat;
        peak.use = obj@peak;
    }else if(input.mat == "gmat"){
        data.use = obj@gmat;        
        peak.use = GenomicRanges::GRanges();
    }else{
        stop("input.mat does not exist in obj")
    }
    
    if((x = nrow(data.use)) == 0L){
        stop("input matrix is empty");
    }else{
        if((x = max(data.use)) > 1L){
            stop("input.mat is not a binary matrix")
        }
    }
    
    # check if empty rows exist in bmat;
    if(any(Matrix::rowSums(data.use) == 0)){
        stop("input matrix contains empty rows, remove empty rows first")    
    }
    
    rm(data.use); # free memory
    rm(peak.use); # free memory
    
    message("Epoch: computing jaccard similarity matrix ...");
    obj = runJaccard2(obj, obj, input.mat=input.mat);

    message("Epoch: fitting regression model ...");
    obj = trainRegression(obj);

    message("Epoch: performing normalization ...");
    obj = normJaccard(obj, obj@regModel[1], obj@regModel[2], obj@regModel[3]);

    # remove the outliers
    nmat.cutoff = quantile(obj@jmat@nmat, nmat.outlier);
    obj@jmat@nmat[obj@jmat@nmat > nmat.cutoff] = nmat.cutoff;
    
    message("Epoch: computing eigen decomposition ...");
    obj = runEigDecomp2(obj, num.eigs);

    obj@smat@method = "runASEMaps";
    message("Epoch: Done");
    return(obj);
}

runJaccard2 <- function(
    obj1,
    obj2,
    input.mat=c("bmat", "pmat", "gmat")
){    
    calJaccard <- function(X_i, X_j){
        A = Matrix::tcrossprod(X_i, X_j);
        bi = Matrix::rowSums(X_i);
        bj = Matrix::rowSums(X_j);
        jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
        rm(A);
        rm(bi);
        rm(bj);
        gc();    
        return(jmat);                
    }
    
    input.mat = match.arg(input.mat);    
    if(input.mat == "bmat"){
        data.use1 = obj1@bmat;
        data.use2 = obj2@bmat;
    }else if(input.mat == "pmat"){
        data.use1 = obj1@pmat;
        data.use2 = obj2@pmat;
    }else if(input.mat == "gmat"){
        data.use1 = obj1@gmat;
        data.use2 = obj2@gmat;
    }
    
    obj1@jmat@jmat = calJaccard(data.use1, data.use2);
    obj1@jmat@p1 = Matrix::rowMeans(data.use1);
    obj1@jmat@p2 = Matrix::rowMeans(data.use2);
    obj1@jmat@norm = FALSE;
    obj1@jmat@nmat = matrix(0,0,0);
    obj1@jmat@input.mat = input.mat;
    return(obj1);
}

trainRegression <- function(obj){
    row.covs = log(Matrix::rowSums(obj@bmat)+1,10);        
    row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
    sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
    idx.ds <- sort(sample(x = seq(row.covs), size = min(1000, length(row.covs)), prob = sampling_prob));
    jmat.tr = obj@jmat@jmat[idx.ds,idx.ds];
    b1.tr = obj@jmat@p1[idx.ds];
    b2.tr = obj@jmat@p2[idx.ds];
    # calculate the expected jaccard index matrix given the read depth
    emat.tr = .normOVE(b1.tr, b2.tr);
    # estimate the global scaling factor
    data = data.frame(x=emat.tr[upper.tri(emat.tr)], y=jmat.tr[upper.tri(jmat.tr)])    
    model <- lm(y ~ x + I(x^2), data);
    beta0 = as.numeric(model$coefficients)[1]
    beta1 = as.numeric(model$coefficients)[2]
    beta2 = as.numeric(model$coefficients)[3]
    obj@regModel = c(beta0, beta1, beta2);
    rm(jmat.tr);
    rm(emat.tr);
    rm(data);
    rm(model);
    rm(row.covs);
    return(obj)
}

normJaccard <- function(obj, beta0, beta1, beta2){
    b1.te = obj@jmat@p1;
    b2.te = obj@jmat@p2;
    jmat.te = obj@jmat@jmat;
    emat.te = .normOVE(b1.te, b2.te);
    preds = beta0 + beta1 * emat.te + beta2 * (emat.te ** 2);
    nmat.te = jmat.te/preds;
    obj@jmat@nmat = nmat.te;
    obj@jmat@norm = TRUE;
    rm(jmat.te);
    rm(emat.te);
    rm(nmat.te);
    rm(preds);
    return(obj);
}

.normOVE <- function(p1, p2){
    pp = tcrossprod(p1, p2);
    ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
    ee = pp/(ss - pp)
    return(ee)    
}


