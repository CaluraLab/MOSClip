#' Get available Omics Summarizing Methods
#'
#' Gives a vector of the available methods to summarize omics.
#'
#' @return character vector with the implemented methods.
#'
#' @examples
#' availableOmicMethods()
#' 
#' @export
availableOmicMethods <- function() {
  return(c("summarizeToBinaryEvents",
           "summarizeToNumberOfEvents",
           "summarizeWithPca",
           "summarizeInCluster",
           "summarizeToBinaryDirectionalEvents",
           "summarizeToNumberOfDirectionalEvents"))
}

#' Summarize To Binary Events
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @param data a data matrix
#' @param features a vector with the features to analyze
#' @param name prefix of the covariates
#' @param binaryClassMin the minimum number of event to include the covariate
#' @param cliques the features organized in cliques. Only use for topology
#' 
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#'  \item{evenThr}{threshold fot event counting}

summarizeToBinaryEvents <- function(data, features, name="bin",
                                    binaryClassMin=10, cliques=NULL) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes) == 0)
    return(NULL)

  dataClique <- t(data[genes, , drop=FALSE])
  if (ncol(dataClique) == 0)
    return(NULL)

  collapsed <- apply(dataClique>0, 1, any, na.rm=TRUE)

  if (sum(collapsed) < binaryClassMin | 
      sum(collapsed) > NROW(dataClique)-binaryClassMin)
    return(NULL)

  collapsed <- data.frame(collapsed, row.names = names(collapsed), 
                          stringsAsFactors = FALSE)
  colnames(collapsed) <- name
  list(x=collapsed, usedGenes=genes, namesCov=name,
       method="binary", omicName=name, eventThr=1)
}

#' Summarize To Number of Binary Events
#'
#' Given a matrix it summarize to a 0 or 1
#'
#' @inheritParams summarizeToBinaryEvents
#' @param min_prop minimal proportion in classes
#'
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#'  \item{evenThr}{threshold fot event counting}
#'  \item{min_prop}{minimum proportion of samples to exclude to check the 
#'  variability of values}

summarizeToNumberOfEvents <- function(data, features, name="event", 
                                      min_prop=0.1, cliques=NULL){
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes) == 0)
    return(NULL)

  dataClique <- t(data[genes, , drop=FALSE])
  if (ncol(dataClique) == 0)
    return(NULL)

  collapsed <- apply(dataClique>0, 1, sum, na.rm=TRUE)

  keep <- check_minimal_proportion(collapsed, min_prop=min_prop)
  if (!keep)
    return(NULL)

  collapsed <- data.frame(collapsed, row.names = names(collapsed), 
                          stringsAsFactors = FALSE)
  colnames(collapsed) <- name
  list(x=collapsed, usedGenes=genes, namesCov=name, method="count", 
       omicName=name, eventThr = 1, min_prop=min_prop)
}

#' Summarize Using Cluster Analysis
#'
#' Given a matrix it summarize in classes
#'
#' The user can define a maximum of classes. The function
#' guess the optimal number of clusters using NbClust methods.
#'
#' @inheritParams summarizeToBinaryEvents
#' @param dictionary translate features (genes) into sets 
#' (row.names of the data)
#' @param max_cluster_number the maximum number of cluster to evaluate
#'
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{cls}{the genes in clusters}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#' 
#' @importFrom stats cutree dist hclust
#' @importFrom NbClust NbClust

summarizeInCluster <- function(data, features, name="clu",
                               dictionary=NULL, max_cluster_number=3, 
                               cliques=NULL) {
  if (is.null(data) || (ncol(data) == 0) || !(is.matrix(data)))
    return(NULL)

  if (is.null(dictionary)) {
    genes <- intersect(rownames(data), features)
    if (length(genes) == 0)
      return(NULL)
    used <- genes
    names(used) <- genes
    datamatClique <- t(data[genes, ,drop=FALSE])
  } else {
    genes <- intersect(names(dictionary), features)
    if (length(genes) == 0)
      return(NULL)
    used <- dictionary[genes]
    clusters <- unlist(dictionary[genes])
    clusters <- intersect(clusters, row.names(data))
    if (length(clusters) == 0)
      return(NULL)
    datamatClique <- t(data[clusters, , drop=FALSE])
  }

  if (ncol(datamatClique) == 0)
    return(NULL)

  ## CREATE CLUSTERS
  covs <- createOptiomalClusterClasses(datamatClique, name, 
                                       max_cluster_number = max_cluster_number)

  if (any(table(covs[[1]])<2)){
    warning("Not meaningful class separation\n")
    return(NULL)
  }

  collapse <- covs
  list(x=collapse, usedGenes=names(used), namesCov=names(covs), cls=used,
       method="cluster", omicName=name)
}

createOptiomalClusterClasses <- function(datamatClique, name,
                                         max_cluster_number,
                                         index_method="silhouette") {
  nb <- sinkNbClust(data=datamatClique, min.nc=2, max.nc=max_cluster_number,
                    method="ward.D2", index=index_method)
  covs <- data.frame(factor(nb$Best.partition), stringsAsFactors = TRUE)
  optimalCLusterNumber <- length(table(nb$Best.partition))
  names(covs) <- paste0(name, optimalCLusterNumber, "k")
  covs
}

#' @importFrom grDevices dev.off pdf
#'
sinkNbClust <- function(data, min.nc=2, max.nc=6, method="ward.D2",
                        index="silhouette"){
  if (index == "all")
    sink(file=tempfile()); pdf(file=NULL)

  nb <- NbClust(data=data, min.nc=min.nc, max.nc=max.nc, method=method,
                index=index)

  if (index == "all")
    sink(); dev.off()

  return(nb)
}


#' Summarize Using PCA
#'
#' Given a matrix it summarize to principal components. 
#' The user can specify the number of principal components. Default 3.
#'
#' @inheritParams summarizeToBinaryEvents
#' @param shrink shirnk or not the covariance matrix.
#' @param method either "regular", "sparse" or "topological"
#' @param cliques the features organized in cliques. Only use for topology.
#' @param maxPCs maximum number of pcs to consider
#' @param loadThr loading threshold
#'
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample (principal components)}
#'  \item{sdev}{standard deviation of the principal components}
#'  \item{loadings}{loadings of PCA}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#'
#' @importFrom stats sd

summarizeWithPca <- function(data, features, name="pca", shrink=FALSE,
                             method="regular", cliques=NULL, maxPCs=3,
                             loadThr=0.6) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes) == 0)
    return(NULL)

  dataClique <- t(data[genes, , drop=FALSE])
  if (ncol(dataClique) == 0)
    return(NULL)

  if (NCOL(dataClique) != 1) {
    pcs <- computePCs(dataClique, shrink=shrink, method=method,
                      cliques=cliques, maxPCs=maxPCs)
    colnames(pcs$x) <- paste0(name, colnames(pcs$x))
    names(pcs$sdev) <- paste0(name, names(pcs$sdev))
    colnames(pcs$loadings) <- paste0(name, colnames(pcs$loadings))
  } else {
    colnames(dataClique) <- paste0(name, "PC1")
    pcs <- list(x=dataClique, sdev=sd(dataClique), loadings=1)
  }

  pcs$usedGenes <- genes
  pcs$method <- "pca"
  pcs$namesCov <- colnames(pcs$x)
  pcs$omicName <- name
  pcs
}

#' Summarize With Directed Sum
#'
#' Given a matrix it summarize the positive and negative in two vectors, 
#' with counts of the events
#'
#' @inheritParams summarizeToNumberOfEvents
#' @param eventThr the absolute value to threshold an event
#'
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#'  \item{evenThr}{threshold fot event counting}
#'  \item{min_prop}{minimum proportion of samples to exclude to check the 
#'  variability of values}

summarizeToNumberOfDirectionalEvents <- function(data, features, name="dCount",
                                                 eventThr=2, min_prop=0.1, 
                                                 cliques=NULL) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes) == 0)
    return(NULL)

  dataClique <- t(data[genes, , drop=FALSE])
  if (ncol(dataClique) == 0)
    return(NULL)

  posDataClique <- extractPositivePortion(dataClique)
  negDataClique <- extractPositivePortion(dataClique, invert=TRUE)

  positive <- apply(posDataClique >= eventThr, 1, sum, na.rm=TRUE)
  negative <- apply(negDataClique >= eventThr, 1, sum, na.rm=TRUE)

  collapsed <- data.frame(positive=positive, negative=negative,
                          row.names = names(positive), stringsAsFactors = FALSE)
  colnames(collapsed) <- paste0(name, c("POS","NEG"))

  keep <- vapply(collapsed, function(x) {
    check_minimal_proportion(x, min_prop=min_prop)},
                 logical(1))
  collapsed <- collapsed[, keep, drop=FALSE]

  if (NCOL(collapsed) == 0)
    return(NULL)
  list(x=collapsed, usedGenes=genes, namesCov=names(collapsed),
       method="directedCount", omicName=name, eventThr=eventThr,
       min_prop=min_prop)
}

#' @importFrom stats quantile
check_minimal_proportion <- function(x, min_prop=0.1){
  min <- quantile(x, probs=c(min_prop))
  max <- quantile(x, probs=c(1-min_prop))
  if ((min == min(x)) && (max == min(x)))
    return(FALSE)

  if ((min == max(x)) && (max == max(x)))
    return(FALSE)

  TRUE
}

#' Summarize To Binary Directional Events
#'
#' Given a matrix it summarize the positive and negative to 0 or 1
#'  in two vectors
#'
#' @inheritParams summarizeToNumberOfDirectionalEvents
#' @inheritParams summarizeToBinaryEvents
#'
#' @return a list with summary of the omic:
#'  \item{x}{summary of the omic for each sample}
#'  \item{usedGenes}{genes list of genes used to calculate the summary}
#'  \item{namesCov}{names of the covariates}
#'  \item{method}{method used for the analysis}
#'  \item{omicName}{name of the omic}
#'  \item{evenThr}{threshold fot event counting}
#'  
summarizeToBinaryDirectionalEvents <- function(data, features, name="dirBin",
                                               binaryClassMin=10, eventThr=2,
                                               cliques=NULL) {
  if (is.null(data))
    return(NULL)

  genes <- intersect(row.names(data), features)
  if (length(genes) == 0)
    return(NULL)

  dataClique <- t(data[genes, , drop=FALSE])
  if (ncol(dataClique) == 0)
    return(NULL)

  posDataClique <- extractPositivePortion(dataClique)
  negDataClique <- extractPositivePortion(dataClique, invert=TRUE)

  positive <- apply(posDataClique >= eventThr, 1, any, na.rm=TRUE)
  negative <- apply(negDataClique >= eventThr, 1, any, na.rm=TRUE)

  collapsed <- data.frame(positive=positive, negative=negative,
                          row.names = names(positive), stringsAsFactors = FALSE)
  colnames(collapsed) <- paste0(name, c("POS","NEG"))


  keep <- vapply(collapsed, sum, as.numeric(1)) >= binaryClassMin | 
    vapply(collapsed, sum, as.numeric(1)) <= NROW(dataClique)-binaryClassMin
  collapsed <- collapsed[, keep, drop=FALSE]

  if (NCOL(collapsed) == 0)
    return(NULL)
  list(x=collapsed, usedGenes=genes, namesCov=names(collapsed),
       method="directedBinary", omicName=name, eventThr=eventThr)
}
