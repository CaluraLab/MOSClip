
#' @importFrom methods new
#' @importFrom survival Surv

MOMSurvTest <- function(genes, omicsObj,
                        survFormula = "Surv(days, status) ~",
                        autoCompleteFormula=T, robust=FALSE, include_from_annot=F) {

  # check if topological method has been used
  for (i in seq_along(omicsObj@ExperimentList@listData)) {
    if (omicsObj@modelInfo[i] == "summarizeWithPca") {
      if (!is.null(omicsObj@specificArgs[[i]]$method)) {
        if (omicsObj@specificArgs[[i]]$method=="topological") {
        stop("Topological: not valid method for module analysis.")
      }}
      else {
        message("Method for pca automatically set to sparse")
        omicsObj@specificArgs[[i]]$method="sparse"}
    }
  }

  moView <- createMOMView(omicsObj, genes)
  formula = survFormula
  coxObj <- createCoxObj(omicsObj@colData, moView)

  add_covs <- unlist(lapply(moView, function(mo) {mo$namesCov}))

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

  scox$moView <- moView

  scox
}




