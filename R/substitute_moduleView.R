
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