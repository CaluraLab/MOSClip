#' Compute Multi Omics two-class in Pathways
#'
#' Performs topological two-class analysis using an 'Omics' object.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param classAnnot class annotation: name of the classes, row.names are samples
#' @param baseFormula model to be used for the test
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param nullModel the null model
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#'
#' @return MultiOmicsPathway object
#'
#' @importFrom graph nodes
#' @importFrom methods new is
#' @export
multiOmicsTwoClassesPathwayTest <- function(omicsObj, graph, classAnnot,
                                            baseFormula = "classes ~ ",
                                            autoCompleteFormula=T,
                                            useThisGenes=NULL,
                                            nullModel = "classes ~ 1",
                                            pathName=NULL) {

  if (is.null(pathName) && is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)
  genesToUse <- graph::nodes(graph)
  if (length(genesToUse)== 0)
    stop("There is no nodes on the graph.")

  moduleView <- lapply(seq_along(omicsObj@ExperimentList@listData),
                       function(i) { 
    test <- get(omicsObj@modelInfo[i]) 
    specificArgs <- omicsObj@specificArgs[[i]]

    cliques=NULL
    if (omicsObj@modelInfo[i]=="summarizeWithPca") {
      genesToUse <- intersect(row.names(omicsObj@ExperimentList@listData[[i]]),
                              genesToUse)
      graph <- graph::subGraph(genesToUse, graph)
      cliques <- extractCliquesFromDag(graph)
    } 
    args <- list(data=omicsObj@ExperimentList@listData[[i]],
                 features=genesToUse, cliques=cliques)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)
  })
  moduleView <- moduleView[!sapply(moduleView, is.null)]
  covariates <- lapply(moduleView, function(mo) {
    mo$x
  })

  moduleData <- lapply(moduleView, function(mo) {
    mo$dataModule
  })

  if (is.null(covariates))
    return(NULL)

  if (!identical(row.names(classAnnot), row.names(covariates)))
    stop("Mismatch in covariates and classes annotations rownames.")

  dataTest <- data.frame(classAnnot, covariates)


  if(!(dependentVar %in% colnames(dataTest)))
    stop(paste0("Data does not contain the model dependent variable: ",
                dependentVar))

  twoClasses <- unique(dataTest[,dependentVar])
  if(length(twoClasses) != 2)
    stop(paste0("Classes in column ",dependentVar," are not two: ",twoClasses))

  if(!all(twoClasses == c(0,1))) {
    dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
    dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
    dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
  }

  fullModelFormula = baseFormula
  if (autoCompleteFormula)
    fullModelFormula = paste0(baseFormula,
                             paste(colnames(covariates), collapse="+"))

  res <- suppressWarnings(glmTest(dataTest, fullModelFormula, nullModelFormula))

  new("MultiOmicsPathway", pvalue=res$pvalue, zlist=res$zlist, coxObj=res$data,
      pathView=moduleView, formula=fullModelFormula, title=pathName,
      analysis="twoClass")
}
# no graphNEL slot, add it?
###****************************fine Pathway***************************


#' Compute Multi Omics two-class in Pathway Modules
#'
#' Performs topological two-class analysis using an 'Omics' object on pathway module.
#'
#' @param omicsObj Object of class 'Omics'
#' @param graph a pathway in graphNEL, Pathway or geneset format.
#' @param classAnnot class annotation: name of the classes, row.names are samples
#' @param baseFormula model to be used for the test
#' @param autoCompleteFormula logical. If TRUE autocomplete the survFormula using all the available covariates
#' @param useThisGenes vector of genes used to filter pathways
#' @param nullModel the null model
#' @param pathName title of the pathway. If NULL and graph is "Pathway" graph@title is used as title
#'
#' @return MultiOmicsModule object
#'
#' @importFrom graph nodes
#' @importFrom methods new is
#' @export
multiOmicsTwoClassesModuleTest <- function(omicsObj, graph, classAnnot,
                                           baseFormula = "classes ~ ",
                                           autoCompleteFormula=TRUE,
                                           useThisGenes=NULL,
                                           nullModel = "classes ~ 1",
                                           pathName=NULL) {

  if (is(graph, "character"))
    stop("Module test can not handle gene list.")

  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)

  genes <- graph::nodes(graph)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and
         the node names in the graph.")

  # create the modules
  cliques <- extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMglmTest, omicsObj=omicsObj,
                    classAnnot=classAnnot,
                    baseFormula=baseFormula,
                    autoCompleteFormula=autoCompleteFormula,
                    nullModel=nullModel)

  alphas   <- as.numeric(sapply(results, extractPvalues))
  zlists    <- lapply(results, function(x) x$zlist)
  coxObjs <- lapply(results, function(x) x$dataTest)
  momics   <- lapply(results, function(x) x$moView)
  formulas <- lapply(results, function(x) x$formula)
  analysis <- "twoClass"

  names(alphas) <- NULL
  new("MultiOmicsModules", alphas=alphas, zlists=zlists, coxObjs=coxObjs,
      modulesView=momics, modules=cliques, formulas=formulas, title=pathName,
      analysis=analysis)
}
# NO graphNEL slot