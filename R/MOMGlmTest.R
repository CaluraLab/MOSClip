#' Two-classes glm test.
#'
#' @param data data
#' @param fullModelFormula complete model
#' @param nullModelFormula null model formula
#' @return Two-classes glm test results
#' @importFrom stats glm poisson pchisq deviance df.residual na.omit

glmTest_old <- function(data, fullModelFormula, nullModelFormula) {
    glmRes <- glm(as.formula(fullModelFormula),
        family = "binomial", data = na.omit(data)
    )
    
    glmSummary <- summary(glmRes)
    zlist <- glmSummary$coefficients[, "Pr(>|z|)"][-1]
    names(zlist) <- row.names(glmSummary$coefficients)[-1]

    fullModel <- glm(as.formula(fullModelFormula),
        family = poisson, data = data
    )
    nullModel <- glm(as.formula(nullModelFormula),
        family = poisson, data = data
    )

    pvalue <- pchisq(deviance(nullModel) - deviance(fullModel),
        df.residual(nullModel) - df.residual(fullModel),
        lower.tail = FALSE
    )

    return(list(pvalue = pvalue, zlist = zlist))
}

#' Two-classes glm test with mixed effects.
#'
#' @param data data
#' @param fullModelFormula complete model
#' @param nullModelFormula null model formula
#' @param mixed use mixed effects model
#' @return Two-classes glm test results
#' @importFrom stats glm poisson pchisq deviance df.residual na.omit
#' @importFrom lme4 glmer

glmTest <- function(data, fullModelFormula, nullModelFormula, mixed = FALSE) {
  args <- list(formula = as.formula(fullModelFormula), family = "binomial", 
               data = na.omit(data))
  model <- ifelse(mixed == TRUE, "glmer", "glm")
  
  glmRes <- do.call(model, args)

  glmSummary <- summary(glmRes)
  
  coeffs <- glmSummary$coefficients
  #if (mixed == TRUE) {
  #  coeffs <- glmSummary$coefficients$cond
  #} else {
  #  coeffs <- glmSummary$coefficients
  #}
  zlist <- coeffs[, "Pr(>|z|)"][-1]
  names(zlist) <- row.names(coeffs)[-1]
  
  args <- list(formula = as.formula(fullModelFormula), family = "poisson", 
               data = data)
  fullModel <- do.call(model, args)
  
  args <- list(formula = as.formula(nullModelFormula), family = "poisson", 
               data = data)
  nullModel <- do.call(model, args)
  
  pvalue <- pchisq(deviance(nullModel) - deviance(fullModel),
                   df.residual(nullModel) - df.residual(fullModel),
                   lower.tail = FALSE
  )
  
  return(list(pvalue = pvalue, zlist = zlist))
}




checkCorrelation <- function(glm, thr){
  vif_res <- car::vif(glm)
  if (sum(vif_res > thr) >= 2){
    high_vif <- vif_res[vif_res > thr]
    message("high vif for covariates: ", 
            paste(names(high_vif), round(high_vif, 3), 
                  sep = "=", collapse = ", "))
  }
}


#' @importFrom methods new
MOMglmTest <- function(
    genes, omicsObj, classAnnot, baseFormula = "classes ~ ",
    autoCompleteFormula = TRUE,
    nullModel = "classes ~ 1", patientCol = "patient",
    mixed = FALSE) {
    # check if topological method has been used
    for (i in seq_along(omicsObj@ExperimentList@listData)) {
        if (omicsObj@modelInfo[i] == "summarizeWithPca") {
            if (!is.null(omicsObj@specificArgs[[i]]$method)) {
                if (omicsObj@specificArgs[[i]]$method == "topological") {
                    stop("Invalid method for module analysis: topological")
                }
            }
        }
    }

    moView <- createMOMView(omicsObj, genes)

    additionalCovariates <- lapply(moView, function(mo) {
        mo$x
    })

    additionalCovariates <- do.call(cbind, additionalCovariates)

    if (is.null(additionalCovariates)){
      # | NCOL(additionalCovariates) < 2) {
        return(NULL)
    } 

    if (!identical(row.names(classAnnot), row.names(additionalCovariates))) {
        if (all(row.names(classAnnot) %in% row.names(additionalCovariates))) {
            res <- resolveAndOrder(list(
                classAnnot = classAnnot,
                additionalCovariates = additionalCovariates
            ))
            classAnnot <- res$classAnnot
            additionalCovariates <- res$additionalCovariates
        } else {
            stop("Mismatch in covariates and annotations row names.")
        }
    }

    dataTest <- data.frame(classAnnot, additionalCovariates)

    nullModelFormula <- nullModel
    if (mixed == TRUE) {
      if (grepl(patientCol, nullModelFormula) == FALSE){
        effect <- paste0("(1|", patientCol, ")")
        nullModelFormula <- paste0(c(nullModel, effect), collapse = "+")
      }
    }

    dependentVar <- all.vars(as.formula(nullModelFormula))[1]

    if (!(dependentVar %in% colnames(dataTest))) {
        stop("Data does not contain one of the model dependent variables")
    }

    twoClasses <- unique(dataTest[, dependentVar])
    if (length(twoClasses) != 2) {
        stop(
            "Classes should be only two. ",
            "Check your dependent variables columns"
        )
    }

    dataTestOut <- dataTest

    if (!(all(twoClasses %in% c(0, 1)) & is.numeric(twoClasses))) {
        dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
        dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
        dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
    }

    fullModelFormula <- baseFormula
    if (autoCompleteFormula) {
      if (mixed == TRUE){
        patient <- omicsObj@colData[, patientCol]
        dataTest <- cbind(patient, dataTest)
        effect <- paste0("(1|", patientCol, ")")
        fullModelFormula <-  paste0(baseFormula, 
                                    paste(c(colnames(additionalCovariates), effect),
                                          collapse = "+"))
      }
      else{
        fullModelFormula <- paste0(baseFormula, 
                                   paste(colnames(additionalCovariates), 
                                                      collapse = "+"))
      }}

    res <- glmTest(dataTest, fullModelFormula, nullModelFormula, mixed)

    res$moView <- moView
    return(res)
}
