#' Create the list of covariates that are going to be tested
#'
#' @importFrom methods new
#' @importFrom survival Surv
#' @return list with
#'   1 reduced representation of the omics
#'   2 sdev
#'   3 loadings or eigenvector
#'   4 usedGenes
#'   5 method
#'   6 namesCov
#'   7 omicName
#'
createMOMView <- function(omicsObj, genes) {
  listCovariates <- lapply(seq_along(omicsObj@ExperimentList@listData), function(i) {
    test <- get(omicsObj@modelInfo[i])
    specificArgs <- omicsObj@specificArgs[[i]]
    args <- list(data=omicsObj@ExperimentList@listData[[i]], features=genes)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)
    do.call(test, args)
  })
  
  listCovariates[!vapply(listCovariates, is.null, logical(1))] 
}


#' Create Cox Object
#' 
#' Create the coxObj from the covariates used in the test
#'
#' @param colData colData from multiOmic object
#' @param moView modulesView or pathView from multiOmicsModules or multiOmicsPathway object
#'
#' @return data.frame, samples in the rows, covariates in the columns
createCoxObj <- function(colData, moView){
  
  additionalCovariates <- lapply(moView, function(mo) {mo$x})
  additionalCovariates <- do.call(cbind, additionalCovariates)
  
  if (is.null(additionalCovariates))
    return(NULL)
  
  if (!identical(row.names(colData), row.names(additionalCovariates)))
    stop("Mismatch in covariates and daysStatus annotations rownames.")
  
  coxObj <- data.frame(colData, additionalCovariates)
  return(coxObj)
}


#' Create Data Module
#' 
#' Extract sub-matrix for the genes of a module or pathway from data matrix of a specific omic
#'
#' @param omic modulesView or pathView object 
#' @param multiOmicObj object of class 'Omics'
#' 
#' @return matrix, genes in the rows, samples in the columns
createDataModule <- function(omic, multiOmicObj){
  if (omic$omicName == "met"){
    genes <- unname(unlist(omic$cls))
  } else {
    genes <- omic$usedGenes
  }
  assay <- assay(multiOmicObj, omic$omicName)
  dataModule <- assay[genes,]
  return(dataModule)
}