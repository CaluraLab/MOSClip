#' Guess the most influent features from MultiOmics Survival or Two-class 
#' results.
#'
#' Given a pathway analyzed by `multiOmicsModuleSurvivalTest` or 
#' `multiOmicsTwoClassModuleTest`, it retrieves for each omic the most 
#' influent features.
#'
#' @param pathway `MultiOmicsModules` object from a pathway
#' @param moduleNumber the module number
#' @param loadThr the loading threshold to select genes (PCA only)
#' @param n the maximum number of genes to retrieve (cluster and binary only)
#' @param atleast the minimum number of features to select (PCA only)
#' @param min_prop_pca the minimal proportion to compute the PCA classes
#' @param min_prop_events the minimal proportion to compute the event classes
#' @param ... additional arguments passed to `get` function
#'
#' @return a list. Each item of the list corresponds to an omic that is 
#' summarized with the specific 'extractSummary' functions. Each item is the 
#' summary for an omic summarized using the setted method: pvalues are present 
#' only for cluster method.
#'
guessInvolvement <- function(
    pathway, moduleNumber, loadThr = 0.6, n = 3, atleast = 1,
    min_prop_pca = 0.1, min_prop_events = 0.1, ...
) {
    multiOmicObj <- get(pathway@multiOmicObj, ...)
    omics <- pathway@modulesView[[moduleNumber]]
    moduleCox <- createCoxObj(multiOmicObj@colData, moView = omics)
    analysis <- pathway@analysis

    lapply(
        omics, function(omic) {
            if (omic$method == "pca") {
                extractSummaryFromPCA(
                  omic, multiOmicObj, moduleCox, analysis, loadThr,
                  atleast, minprop = min_prop_pca
              )
            } else if (omic$method == "cluster") {
                extractSummaryFromCluster(omic, multiOmicObj, n)
            } else if (omic$method %in% c("binary", "directedBinary")) {
                extractSummaryFromBinary(omic, multiOmicObj, n)
            } else if (omic$method %in% c("count", "directedCount")) {
                extractSummaryFromNumberOfEvents(
                  omic, multiOmicObj, moduleCox, analysis, n = 3,
                  minprop = min_prop_events
              )
            } else {
                stop("Unsupported method.")
            }
        }
    )
}

#' Guess the most influent features from MultiOmics Survival or Two-class 
#' results.
#'
#' Given a pathway analyzed by `multiOmicsSurvivalPathwayTest` or 
#' `multiOmicsTwoClassPathwayTest`, it retrieves for each omic the most 
#' influent features.
#'
#' @inheritParams guessInvolvement
#'
#' @inherit guessInvolvement return
#'
guessInvolvementPathway <- function(
    pathway, loadThr = 0.6, n = 3, atleast = 1, min_prop_pca = 0.1,
    min_prop_events = 0.1, ...
) {

    multiOmicObj <- get(pathway@multiOmicObj, ...)
    omics <- pathway@pathView
    moduleCox <- createCoxObj(multiOmicObj@colData, moView = omics)
    analysis <- pathway@analysis

    lapply(
        omics, function(omic) {
            if (omic$method == "pca") {
                extractSummaryFromPCA(
                  omic, multiOmicObj, moduleCox, analysis, loadThr,
                  atleast, minprop = min_prop_pca
              )

            } else if (omic$method == "cluster") {
                extractSummaryFromCluster(omic, multiOmicObj, n)
            } else if (omic$method %in% c("binary", "directedBinary")) {
                extractSummaryFromBinary(omic, multiOmicObj, n)
            } else if (omic$method %in% c("count", "directedCount")) {
                extractSummaryFromNumberOfEvents(
                  omic, multiOmicObj, moduleCox, analysis, n = 3,
                  minprop = min_prop_events
              )
            } else {
                stop("Unsupported method.")
            }
        }
    )
}
