#' Cox Model Analysis
#'
#' Cox Analysis
#'
#' For internal use only
#'
#' @param coxObj data.frame: patients x covariates
#' @param formula formula to use
#'
#' @return A list with
#' \item{pvalue}{pvalue of the model}
#' \item{zlist}{pvalues of single covariates}
#' \item{coxObj}{the original coxObj passed to the function}
#'
#' @importFrom stats as.formula na.omit
#' @importFrom survival coxph Surv
#'
survivalcox <- function(coxObj, formula){
  originalCoxObj <- coxObj
  coxObj <- na.omit(coxObj)
  coxRes <- survival::coxph(as.formula(formula), na.omit(coxObj))
  coxSummary <- summary(coxRes)
  zlist <- coxSummary$coefficients[,"Pr(>|z|)"]
  names(zlist) <- row.names(coxSummary$coefficients)
  pvalue <- coxSummary$logtest["pvalue"]
  return(list(pvalue=pvalue, zlist=zlist, coxObj=originalCoxObj))
}

#' Cox Robust Model Analysis
#'
#' Cox Robust Analysis
#'
#' For internal use only
#'
#' @param coxObj data.frame: patients x covariates
#' @param formula formula to use
#'
#' @return A list with
#' \item{pvalue}{pvalue of the model}
#' \item{zlist}{pvalues of single covariates}
#' \item{coxObj}{the original coxObj passed to the function}
#'
#' @importFrom stats as.formula na.omit
#' @importFrom survival Surv
#' @importFrom coxrobust coxr
#'
survivalcoxr <- function(coxObj, formula){
  originalCoxObj <- coxObj
  coxObj <- na.omit(coxObj)
  coxRes <- coxrobust::coxr(as.formula(formula), na.omit(coxObj))
  coxSummary <- coxrsummary(coxRes)
  zlist <- coxSummary$robust[,"p"]
  names(zlist) <- row.names(coxSummary$robust)
  pvalue <- coxSummary$extWakd
  return(list(pvalue=pvalue, zlist=zlist, coxObj=originalCoxObj))
}

#' Summarize Cox Robust
#'
#' @param x a coxr.obj
#'
#' @return a list with wald test and robust and partial coefficients
#'
#' @rdname survivalcoxr
#'
#' @importFrom stats pchisq pnorm
#'
coxrsummary <- function(x) {
  sd.ple <- sqrt( diag( x$var.ple ) )
  sd <- sqrt(diag(x$var))
  df <- sum(!is.na(x$coef))

  tmp <- cbind(x$ple.coefficients, exp(x$ple.coefficients), sd.ple,
               2*(1-pnorm(abs(x$ple.coefficients/sd.ple))))

  dimnames(tmp)[[1]] <- names(x$coef)
  dimnames(tmp)[[2]] <- c("coef", "exp(coef)", "se(coef)", "p")

  partialLikelihood <- tmp
  wald <- 1 - pchisq(x$wald.test, df)

  tmp <- cbind(x$coef, exp(x$coef), sd, 2*(1-pnorm(abs(x$coef/sd))))
  dimnames(tmp)[[1]] <- names(x$coef)
  dimnames(tmp)[[2]] <- c("coef", "exp(coef)", "se(coef)", "p")

  robustEstimator <- tmp
  extWald <- 1 - pchisq(x$ewald.test, df)

  list(wald=wald, partial=partialLikelihood, extWakd=extWald,
       robust=robustEstimator)
}
