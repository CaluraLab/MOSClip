
#' Two-classes glm test.
#'
#' @param data data
#' @param fullModelFormula complete model
#' @param nullModelFormula null model formula
#' @return Two-classes glm test results
#' @importFrom stats glm poisson pchisq deviance df.residual na.omit
glmTest <- function(data, fullModelFormula, nullModelFormula){

  glmRes <- glm(as.formula(fullModelFormula), family="binomial",
                data=na.omit(data))
  glmSummary <- summary(glmRes)
  zlist <- glmSummary$coefficients[,"Pr(>|z|)"][-1]
  names(zlist) <- row.names(glmSummary$coefficients)[-1]

  # test
  fullModel=glm(as.formula(fullModelFormula), family=poisson, data=data)
  nullModel=glm(as.formula(nullModelFormula), family=poisson, data=data)
  pvalue <- pchisq(deviance(nullModel)-deviance(fullModel),
                   df.residual(nullModel)-df.residual(fullModel),
                   lower.tail=FALSE)

  return(list(pvalue=pvalue, zlist=zlist))
}

#' @importFrom methods new
MOMglmTest <- function(genes, omicsObj, classAnnot,
                       baseFormula = "classes ~ ",
                       autoCompleteFormula=TRUE,
                       nullModel = "classes ~ 1") {

  # check if topological method has been used
  for (i in seq_along(omicsObj@ExperimentList@listData)) {
    if (omicsObj@modelInfo[i] == "summarizeWithPca") {
      if (!is.null(omicsObj@specificArgs[[i]]$method)) {
        if (omicsObj@specificArgs[[i]]$method=="topological") {
          stop("Invalid method for module analysis: topological")
        }}
      else if (is.null(omicsObj@specificArgs[[i]]$method)) {
        omicsObj@specificArgs[[i]]$method="sparse"}
    }
  }

  moView <- createMOMView(omicsObj, genes)

  additionalCovariates <- lapply(moView, function(mo) {
    mo$x
  })
  moduleData <- lapply(moView, function(mo) {
    mo$dataModule
  })

  additionalCovariates <- do.call(cbind, additionalCovariates)

  if (is.null(additionalCovariates))
    return(NULL)

  if (!identical(row.names(classAnnot), row.names(additionalCovariates))) {
    if (all(row.names(classAnnot) %in% row.names(additionalCovariates))) {
      res <- resolveAndOrder(list(classAnnot = classAnnot,
                                  additionalCovariates = additionalCovariates))
      classAnnot = res$classAnnot
      additionalCovariates = res$additionalCovariates
      }
    else {stop("Mismatch in covariates and annotations row names.") }}

  dataTest <- data.frame(classAnnot, additionalCovariates)

  # nullModelFormula <- paste0(baseFormula,"1")
  nullModelFormula <- nullModel

  dependentVar <- all.vars(as.formula(nullModelFormula))[1]
  if(!(dependentVar %in% colnames(dataTest)))
    stop(paste0("Data does not contain the model dependent variable: ",
                dependentVar))

  twoClasses <- unique(dataTest[,dependentVar])
  if(length(twoClasses) != 2)
    stop(paste0("Classes in column ",dependentVar," are not two: ",twoClasses))

  dataTestOut <- dataTest

  if(!all(twoClasses == c(0,1))) {
    dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
    dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
    dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
  }

  fullModelFormula = baseFormula
  if (autoCompleteFormula)
    fullModelFormula = paste0(baseFormula, paste(colnames(additionalCovariates),
                                                 collapse="+"))

  res <- suppressWarnings(glmTest(dataTest, fullModelFormula, nullModelFormula))

  res$moView <- moView
  res
}
