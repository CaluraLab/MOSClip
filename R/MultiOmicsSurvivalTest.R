#' Compute Multi Omics Survival in Pathways
#'
#' Performs topological survival analysis using an 'Omics' object.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param survFormula Formula to compute survival
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#' @param robust should be used the robust mode for cox
#' @param include_from_annot compute cox using additional covariates from annot
#'
#' @return MultiOmicsPathway object
#'
#' @importFrom graph nodes
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom methods new is
#' @importFrom survival Surv
#' @importFrom survClip survivalcox survivalcoxr
#' @export

multiOmicsSurvivalPathwayTest <- function(omicsObj, graph,
                                          survFormula = "Surv(days, status) ~",
                                          autoCompleteFormula=T,
                                          useThisGenes=NULL,
                                          pathName=NULL, robust=FALSE,
                                          include_from_annot=F) {
  
  if (is.null(pathName) && is(graph, "Pathway")) {
    pathName <- graph@title
  }
  
  graph <- convertPathway(graph, useThisGenes)
  genesToUse <- graph::nodes(graph)
  if (length(genesToUse)== 0)
    stop("There is no nodes on the graph.")
  
  moduleView <- lapply(seq_along(omicsObj@ExperimentList@listData), function(i) {
    test <- omicsObj@modelInfo[i]
    specificArgs <- omicsObj@specificArgs[[i]]
    
    cliques=NULL
    if (test=="summarizeWithPca") {
      genesToUse <- intersect(row.names(omicsObj@ExperimentList@listData[[i]]),
                              genesToUse)
      graph <- graph::subGraph(genesToUse, graph)
      cliques <- houseOfClipUtility::extractCliquesFromDag(graph)
    }
    
    args <- list(data=omicsObj@ExperimentList@listData[[i]],
                 features=genesToUse, cliques=cliques)
    
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)
    
    do.call(test, args)
  })
  
  moduleView <- moduleView[!sapply(moduleView, is.null)]
  
  covariates <- lapply(moduleView, function(mo) mo$x)
  
  moduleData <- lapply(moduleView, function(mo) mo$dataModule)
  
  covariates <- do.call(cbind, covariates)
  
  if (is.null(covariates))
    return(NULL)
  
  if (!identical(row.names(colData(omicsObj)), row.names(covariates)))
    stop("Mismatch in covariates and annot annotations rownames.")
  
  coxObj <- data.frame(omicsObj@colData, covariates)
  
  formula = survFormula
  
  add_covs <- colnames(covariates)
  if (include_from_annot) {
    add_annot_covs <- colnames(omicsObj@colData)[!colnames(omicsObj@colData) %in% c("days", "status")]
    add_covs <- c(add_covs, add_annot_covs)
  }
  
  if (autoCompleteFormula)
    formula = paste0(survFormula, paste(add_covs, collapse="+"))
  
  if (robust) {
    scox <- suppressWarnings(survClip::survivalcoxr(coxObj, formula)) ### Check warnings
  } else {
    scox <- suppressWarnings(survClip::survivalcox(coxObj, formula)) ### Check warnings
  }
  new("MultiOmicsPathway", pvalue=scox$pvalue, zlist=scox$zlist, coxObj=scox$coxObj,
      pathView=moduleView, formula=formula, title=pathName)
}

#' Compute Multi Omics Survival in Pathway Modules
#'
#' Performs survival analysis using an 'Omics' object. The pathway (graph) used is decomposed in modules (cliques) using graph theory.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param survFormula Formula to compute survival
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#' @param robust should be used the robust mode for cox
#' @param include_from_annot compute cox using additional covariates from annot
#'
#' @return MultiOmicsModules object
#'
#' @importFrom graph nodes
#' @importFrom houseOfClipUtility extractCliquesFromDag
#' @importFrom methods new is
#' @importFrom survival Surv
#' 
#' @export
multiOmicsSurvivalModuleTest <- function(omicsObj, graph,
                                         survFormula = "Surv(days, status) ~",
                                         autoCompleteFormula=T, useThisGenes=NULL,
                                         pathName=NULL, robust=FALSE, include_from_annot=F) {
  
  if (is(graph, "character"))
    stop("Module test can not handle gene list.")
  
  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title
  
  graph <- convertPathway(graph, useThisGenes)
  
  genes <- graph::nodes(graph)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  # create the modules
  cliques <- houseOfClipUtility::extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMSurvTest, omicsObj=omicsObj,
                    survFormula = survFormula,
                    autoCompleteFormula=autoCompleteFormula,
                    robust=robust, include_from_annot=include_from_annot)
  
  alphas   <- as.numeric(sapply(results, extractPvalues))
  zlist    <- lapply(results, function(x) x$zlist)
  momics   <- lapply(results, function(x) x$moView)
  coxObjs  <- lapply(results, function(x) x$coxObj)
  moduleData <- lapply(results, function(x) x$moduleData)
  modules  <- cliques
  formulas <- lapply(results, function(x) x$formula)
  
  names(alphas) <- NULL
  new("MultiOmicsModules", alphas=alphas, zlists=zlist, coxObjs=coxObjs,
      modulesView=momics, modules=modules, formulas=formulas,
      title=pathName)
}