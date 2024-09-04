#' compute PCs.
#'
#' For internal use only. Performs Principal Componenent analysis.
#'
#' Three methods are implemented:
#'   * regular: a regular PCA ('prcomp')
#'   * topological: PCA using a pathway topology.
#'   * sparse: sparse PCA analysis implemented by 'elasticnet'
#'
#' @inheritParams topoCompPCs
#' @param method one of 'regular', 'topological' and 'sparse'
#' @param maxPCs the maximum number of PCs to consider
#'
#' @return a list with the following elements:
#'   \item{x}{the computed PCs}
#'   \item{sdev}{the standard deviation captured by the PCs}
#'   \item{loadings}{the loadings}
#'
#' @importFrom FactoMineR estim_ncp
#'
#' @export
#'
computePCs <- function(exp, shrink=FALSE, method=c("regular", "topological",
                                                   "sparse"),
                       cliques=NULL, maxPCs=3) {
  k<- min(FactoMineR::estim_ncp(exp,scale=FALSE,ncp.min=1)$ncp, maxPCs)
  switch(method[1],
         "regular"     = compPCs(exp=exp, shrink=shrink, k=k),
         "topological" = topoCompPCs(exp=exp, shrink=shrink, cliques=cliques, 
                                     k=k),
         "sparse"      = sparseCompPCs(exp=exp, shrink=shrink, k=k))
}

#' Topological PCA
#'
#' @param exp a matrix
#' @param shrink logical, whether to shrink or not.
#' @param cliques the pathway topology summarized in a list of cliques
#' @param k the number of components to use
#'
#' @return a list with the following elements:
#'   \item{x}{the computed PCs}
#'   \item{sdev}{the standard deviation captured by the PCs}
#'   \item{loadings}{the loadings}
#'
#' @importFrom qpgraph qpIPF
#' @export
#'
topoCompPCs <- function(exp, shrink, cliques, k) {
  if (is.null(cliques))
    stop("Cliques argument is needed")
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink)
  covmat <- makePositiveDefinite(covmat)$m1
  cliquesIdx <- lapply(cliques, function(c) match(c, row.names(covmat)))
  if (any(sapply(cliquesIdx, function(x) {any(is.na(x))})))
    stop("Some genes in cliques are not present in expression.")
  scovmat <- qpgraph::qpIPF(covmat, cliquesIdx)
  pcCov <- base::eigen(scovmat)
  eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

#' Sparse PCA
#'
#' @inheritParams topoCompPCs
#'
#' @return a list with the following elements:
#'   \item{x}{the computed PCs}
#'   \item{sdev}{the standard deviation captured by the PCs}
#'   \item{loadings}{the loadings}
#'
#' @rdname computePCs
#' @importFrom elasticnet spca
#' @export
#'
sparseCompPCs <- function(exp, shrink, k) {
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink)
  covmat <- makePositiveDefinite(covmat)$m1
  paraSingle <- min(round((NCOL(exp)/2)),5) ## Parametri fissi da valutare
  pcCov <- elasticnet::spca(covmat, K =k, para = rep(paraSingle,k),
                            type = "Gram", sparse = "varnum")
  eigenvector  <- pcCov$loadings[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

#' Regular PCA
#'
#' @inheritParams topoCompPCs
#'
#' @return a list with the following elements:
#'   \item{x}{the computed PCs}
#'   \item{sdev}{the standard deviation captured by the PCs}
#'   \item{loadings}{the loadings}
#'
#' @rdname computePCs
#' @export
#'
compPCs <- function(exp, shrink, k) {
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink) ## Consider collapse with the following line!
  covmat <- makePositiveDefinite(covmat)$m1
  scovmat<-covmat
  pcCov <- base::eigen(scovmat)
  eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}
