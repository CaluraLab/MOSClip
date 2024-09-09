#' Extract Summary Binary from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeToBinaryEvents" extract the most 
#' important features.
#'
#' @param omic a summarized omic
#' @param multiOmicObj `Omics` object
#' @param n maximum number of features to retrieve
#'
#' @return Meant for internal use only. The summary for omic summarized 
#' using binary events.
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted 
#' (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), cov=sum binary events)}
#' \item{covsConsidered}{the name of the considered omic}
extractSummaryFromBinary <- function(omic, multiOmicObj, n=3) {
  moduleMat <- createDataModule(omic, multiOmicObj)
  covs <- omic$namesCov

  impact <- lapply(covs, mostlyMutated, moduleMat=t(moduleMat),
                   name=omic$omicName,eventThr=omic$eventThr)
  mostlyImpacted <- lapply(impact, head, n=n)
  involved <- unique(unlist(lapply(mostlyImpacted, row.names)))
  sigModule <- moduleMat[involved, , drop=FALSE]

  discrete <- data.frame(lapply(
    omic$x, as.numeric), row.names=row.names(omic$x))
  list(sigModule=sigModule, discrete=discrete, subset=mostlyImpacted,
       covsConsidered=omic$namesCov)
}


#' Extract Summary Cluster from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeInCluster" extract the most 
#' important features.
#'
#' @inheritParams extractSummaryFromBinary
#'
#' @return summary for omic summarized using clusters
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted 
#' (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), metClust=topGenes)}
#' \item{pvalues}{Kruskal Wallis pvalues of the selected features}
#' \item{covsConsidered}{the name of the considered omic}
#' @importFrom utils head
extractSummaryFromCluster <-function(omic, multiOmicObj, n=3) {
  moduleMat <- createDataModule(omic, multiOmicObj)
  classes <- omic$x[,1]
  KMsigMat <- KWtest(t(moduleMat), classes)

  if (is.null(names(omic$cls)))
    stop("cls in omic needs to be a list or a named vector")

  g <- lapply(names(omic$cls), function(gene) {
    cbind(gene, omic$cls[[gene]])
  })
  g <- do.call(rbind, g)
  genes <- g[,1]
  names(genes) <- g[,2]

  involved <- head(KMsigMat, n)
  sigModule <- moduleMat[row.names(involved), , drop=FALSE]
  topGenes <- genes[row.names(involved)]
  topGenes <- tapply(names(topGenes), topGenes, paste, collapse=";")

  list(sigModule=sigModule, discrete=omic$x,
       subset=data.frame(row.names=names(topGenes), metClust=topGenes),
       pvalues=involved,
       covsConsidered=omic$namesCov)
}

#' @importFrom stats kruskal.test
KWtest <- function(moduleMat, classes) {
  kwTest <- apply(moduleMat, 2, function(gene) {
    kruskal.test(x=gene, g=classes)[c("p.value", "statistic")]
  })
  res <- do.call(rbind, kwTest)
  res[order(as.numeric(res[,"p.value"])),, drop=FALSE]
}

#' Extract Summary PCA from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeWithPca" extract the most important 
#' features.
#'
#' @inheritParams extractSummaryFromBinary
#' @param moduleCox the coxObj of the interesting module
#' @param analysis type of analysis: survival or twoClass
#' @param loadThr the thr value to select the most influent features according 
#' to the loading
#' @param atleast the minimum number of gene to retrieve
#' @param minprop the minimal proportion of cutp
#' @param multiOmicObj Omics object
#' @param analysis two-class or survival type
#'
#' @return summary for omic summarized using pca
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted 
#' (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), covariate)}
#' \item{covsConsidered}{the name of the considered omic}
extractSummaryFromPCA <- function(omic, multiOmicObj, moduleCox, analysis, 
                                  loadThr=0.6, atleast=1, minprop=0.1) {
  covs <- omic$namesCov
  lds <- omic$loadings
  discretePC <- createDiscreteClasses(coxObj=moduleCox, covs, analysis,
                                      minprop=minprop)
  topLoad <- extractHighLoadingsGenes(lds, thr=loadThr, atleast=atleast)
  moduleMat <- createDataModule(omic, multiOmicObj)
  sigModule <- moduleMat[row.names(topLoad), , drop=FALSE]

  list(sigModule=sigModule, discrete=discretePC, subset=topLoad,
       covsConsidered=covs)
}

#' @importFrom survminer surv_cutpoint surv_categorize
#' @importFrom stats median
createDiscreteClasses <- function(coxObj, covs, analysis,
                                  labels= c("low", "high"), minprop=0.1) {
  
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }
  
  check <- vapply(coxObj[, covs, drop=FALSE], check_minimal_proportion,
                  c(min_prop=minprop), logical(1))
  if (any(!check)){
    stop("minprop ", minprop, " is too high. Try a smaller one")
  }
  
  if (analysis == "survival") {
    sc <- surv_cutpoint(coxObj, time="days", event="status", variables = covs,
                        minprop=minprop)
    surv_categorize(sc, labels=labels)
  } else if (analysis == "twoClass") {
    covs_classified <- lapply(covs, function(cov) {
      median_value <- median(as.numeric(coxObj[[cov]]), na.rm = TRUE)
      ifelse(coxObj[[cov]] >= median_value, "high", "low")
    })
    sc <- coxObj
    sc[covs] <- covs_classified
    return(sc)
  } else {
    stop("Invalid type of analysis: ", analysis, " Check results object")
  }
}

#' @importFrom survminer surv_cutpoint surv_categorize
retrieveNumericClasses <- function(coxObj, covs, analysis) {
  diff <- setdiff(covs, colnames(coxObj))
  if (length(diff) != 0) {
    stop(paste0(paste(diff, collapse=", "), " not in coxObj."))
  }
  
  if (analysis == "survival") {
    coxObj[, c("days", "status", covs), drop=FALSE]
  } else if (analysis == "twoClass") {
    coxObj[, covs, drop=FALSE]
  } else {
    stop("Invalid type of analysis ", analysis, " Check results object")
  }
}

extractHighLoadingsGenes <- function(loadings, thr, atleast=1) {
  l <- lapply(colnames(loadings), function(pc) {
    ld <- loadings[, pc]
    genes <- names(which(abs(ld) >= thr))

    if (length(genes) == 0)
      genes <- names(ld[order(abs(ld), decreasing = TRUE)][seq_len(atleast)])

    data.frame(row.names=genes, component=rep(pc, length(genes)),
               stringsAsFactors = FALSE)
  })
  l <- collapse(l)
  do.call(rbind, l)
}

collapse <- function(list) {
  df <- data.frame(genes=unlist(lapply(list, function(x) row.names(x))),
                   components=unlist(lapply(list, function(x) x$component)),
                   stringsAsFactors = FALSE)
  tapply(seq_len(NROW(df)), df$genes, function(idx){
    paste(df$components[idx], collapse=";")
  }, simplify = FALSE)
}

#' Extract Summary Binary from MultiOmics Objects
#'
#' Given an omic summarized by "summarizeToNumberOfEvents" extract the most 
#' important features.
#'
#' @inheritParams extractSummaryFromPCA
#' @param labels the category labels
#' @param multiOmicObj Omics object
#' @param analysis two-class or survival type
#'
#' @return Meant for internal use only. The summary for omic summarized using 
#' counting of events.
#' \item{sigModule}{the original data for significant features}
#' \item{discrete}{the discrete version of the significant covariates converted 
#' (when needed) into the discrete version}
#' \item{subset}{data.frame(row.names=names(topGenes), covariates=covariate)}
#' \item{covsConsidered}{the name of the considered omic}
extractSummaryFromNumberOfEvents <- function(omic, multiOmicObj, moduleCox, 
                                             analysis, n=3, minprop=0.1,
                                             labels=c("few","many")) {
  covs <- omic$namesCov
  moduleMat <- createDataModule(omic, multiOmicObj)
  discreteClass <- createDiscreteClasses(coxObj=moduleCox, covs, analysis,
                                          labels=labels, minprop=minprop)
  numericClass <- retrieveNumericClasses(coxObj=moduleCox, covs, analysis)

  impact <- lapply(covs, mostlyMutated, moduleMat=t(moduleMat),
                   name=omic$omicName, eventThr = omic$eventThr)

  mostlyImpacted <- lapply(impact, head, n=n)
  names(mostlyImpacted) <- covs

  covNames <- lapply(covs, function(cov) {
    rep(cov, length(row.names(mostlyImpacted[[cov]])))
  })
  covNames <- do.call(c, covNames)
  mostlyImpactedGenes <- unlist(lapply(mostlyImpacted, row.names))
  sigModule <- moduleMat[mostlyImpactedGenes, , drop=FALSE]

  covariates <- tapply(covNames, mostlyImpactedGenes,  function(x) x)
  involved <- data.frame(covariates, stringsAsFactors = FALSE)

  list(sigModule=sigModule, discrete=discreteClass, subset=involved,
       covsConsidered=omic$namesCov, numericClass=numericClass)
}

mostlyMutated <- function(cov, moduleMat, name, eventThr=2) {
  direction <- gsub(name, "", cov)
  if (direction == "NEG"){
    invert <- TRUE
  } else {
    invert <- FALSE
  }
  priority <- extractPositivePortion(moduleMat, invert = invert)
  priority <- colSums(priority >= eventThr, na.rm=TRUE)
  priority <- data.frame(row.names = names(priority), priority)
  names(priority) <- cov
  priority[order(priority[[cov]], decreasing = TRUE),, drop=FALSE]
}

