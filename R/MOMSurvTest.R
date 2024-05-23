#' Create the list of covariates that are going to be tested
#'
#' @importFrom methods new
#' @importFrom survival Surv
#' @return list with
#'   1 reduced representation of the omics
#'   2 sdev
#'   3 loadings or eigenvector
#'   4 data module ## Consider adding here only genes)
#'   5 Method
#'   6 namesCov
#'   7 OmicName
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

  listCovariates[!sapply(listCovariates, is.null)] # remove sapply
}


#' @importFrom methods new
#' @importFrom survival Surv

MOMSurvTest <- function(genes, omicsObj,
                        survFormula = "Surv(days, status) ~",
                        autoCompleteFormula=T, robust=FALSE, include_from_annot=F) {

  # check if topological method has been used
  for (i in seq_along(omicsObj@ExperimentList@listData)) {
    if (omicsObj@modelInfo[i] == "summarizeWithPca") {
      if (omicsObj@specificArgs[[i]]$method=="topological") {
        stop("Topological: not valid method for module analysis.")
      }
    }
  }

  moView <- createMOMView(omicsObj, genes)
  formula = survFormula

  coxObj <- omicsObj@colData

  additionalCovariates <- lapply(moView, function(mo) {
    mo$x
  })

  moduleData <- lapply(moView, function(mo) mo$dataModule)
  usedGenes <- lapply(moView, function(mo) mo$usedGenes)

  additionalCovariates <- do.call(cbind, additionalCovariates)

  if (is.null(additionalCovariates))
    return(NULL)

  if (!identical(row.names(coxObj), row.names(additionalCovariates)))
    stop("Mismatch in covariates and daysStatus annotations rownames.")


  coxObj <- data.frame(coxObj, additionalCovariates)

  add_covs <- colnames(additionalCovariates)
  if (include_from_annot) {
    add_annot_covs <- colnames(coxObj)[!colnames(coxObj) %in% c("days", "status")]
    add_covs <- c(add_covs, add_annot_covs)
  }

  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(add_covs, collapse="+"))

  if (robust) {
    scox <- suppressWarnings(survivalcoxr(coxObj, formula)) ### Check warnings
  } else {
    scox <- suppressWarnings(survivalcox(coxObj, formula)) ### Check warnings
  }

  scox$moView <- moView # consider removing
  scox$formula <- formula
  scox$moduleData <- moduleData # consider removing
  scox$usedGenes <- usedGenes

  scox
}




