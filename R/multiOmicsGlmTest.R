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

  baseFormula_input <- strsplit(baseFormula, " ")[[1]]
  if (!(baseFormula_input[1] %in% colnames(classAnnot))){
    stop("Invalid formula. Class column not found in classAnnot")
  } else if (length(baseFormula_input) == 1 | baseFormula_input[2] != "~") {
    stop("Invalid formula. Formula should be written as: 'classes ~'")
  }
  
  nullModel_input <- strsplit(nullModel, " ")[[1]]
  if (!(nullModel_input[1] %in% colnames(classAnnot))){
    stop("Invalid formula. Class column not found in classAnnot")
  } else if (length(nullModel_input) == 1) {
    stop("Formula is too short. Formula should be written as 'classes ~ 1")
  } else if (nullModel_input[2] != "~" | nullModel_input[3] != "1") {
    stop("Invalid formula. Formula should be written as: 'classes ~'")
  }
  
  if (is.null(pathName) && is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)
  genesToUse <- graph::nodes(graph)
  if (length(genesToUse)== 0)
    stop("There is no nodes on the graph.")

  moduleView <- lapply(seq_along(omicsObj@ExperimentList@listData), function(i) 
    {
    test <- get(omicsObj@modelInfo[i]) 
    specificArgs <- omicsObj@specificArgs[[i]]

    cliques=NULL
    if (omicsObj@modelInfo[i]=="summarizeWithPca") {
      genesToUse <- intersect(row.names(omicsObj@ExperimentList@listData[[i]]),
                              genesToUse)
      if (identical(genesToUse, character(0))) {
        stop(paste0("This data ", omicsObj@ExperimentList@listData[[i]],
                    "have no genes for this pathway"))
      }
      graph <- graph::subGraph(genesToUse, graph)
      cliques <- extractCliquesFromDag(graph)
    } 
    args <- list(data=omicsObj@ExperimentList@listData[[i]],
                 features=genesToUse, cliques=cliques)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)
    do.call(test, args)
  })
  
  moduleView <- moduleView[!sapply(moduleView, is.null)]
  covariates <- lapply(moduleView, function(mo) {mo$x})
  covariates <- do.call(cbind, covariates)

  if (is.null(covariates))
    return(NULL)

  #if (nrow(classAnnot) != nrow(covariates))
   # warning("Mismatch in the number of samples.")

  if (!identical(row.names(classAnnot), row.names(covariates))) {
    if (all(row.names(classAnnot) %in% row.names(covariates))) {
      res <- resolveAndOrder(list(classAnnot = classAnnot, 
                                  covariates = covariates))
      classAnnot = res$classAnnot
      covariates = res$covariates 
    } 
    else {stop("Mismatch in covariates and classes annotations row names.") }}
  
  dataTest <- data.frame(classAnnot, covariates)
  
  # nullModelFormula <- paste0(baseFormula,"1")
 nullModelFormula <- nullModel
  
  dependentVar <- all.vars(as.formula(nullModelFormula))[1]
  if(!(dependentVar %in% colnames(dataTest)))
    stop(paste0(
      "Data does not contain the model dependent variable: ", dependentVar))

  twoClasses <- unique(dataTest[,dependentVar])
  if(length(twoClasses) != 2)
    stop(paste0(
      "Classes in column ", dependentVar, " are not two: ", twoClasses))

  if(!all(twoClasses == c(0,1))) {
    dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
    dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
    dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
  }

  fullModelFormula = baseFormula
  if (autoCompleteFormula) # i'd remove that or change because we r not doing survival 4 2 class
    fullModelFormula = paste0(baseFormula,
                             paste(colnames(covariates), collapse="+"))

  res <- suppressWarnings(glmTest(dataTest, fullModelFormula, nullModelFormula))

  new("MultiOmicsPathway", 
      pvalue=res$pvalue, 
      zlist=res$zlist, 
      pathView=moduleView, 
      analysis="twoClass",
      multiOmicObj=deparse(substitute(omicsObj)),
      title=pathName)
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
  
  baseFormula_input <- strsplit(baseFormula, " ")[[1]]
  if (!(baseFormula_input[1] %in% colnames(classAnnot))){
    stop("Invalid formula. Class column not found in classAnnot")
  } else if (length(baseFormula_input) == 1 | baseFormula_input[2] != "~") {
    stop("Invalid formula. Formula should be written as: 'classes ~'")
  }
  
  nullModel_input <- strsplit(nullModel, " ")[[1]]
  if (!(nullModel_input[1] %in% colnames(classAnnot))){
    stop("Invalid null formula. Class column not found in classAnnot")
  } else if (length(nullModel_input) == 1) {
    stop("Null formula is too short. Formula should be written as 'classes ~ 1")
  } else if (nullModel_input[2] != "~" | nullModel_input[3] != "1") {
    stop("Invalid null formula. Formula should be written as: 'classes ~'")
  }

  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useThisGenes)
  genes <- graph::nodes(graph)
  if (length(genes) == 0)
    stop("There is no intersection between expression feature names and
         the node names in the graph.")

  # create the modules
  cliques <- extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMglmTest, omicsObj=omicsObj,
                    classAnnot=classAnnot,
                    baseFormula=baseFormula,
                    autoCompleteFormula=autoCompleteFormula,
                    nullModel=nullModel)
  
  #if (nrow(classAnnot) != nrow(multiOmics@colData))
   # warning("Mismatch in the number of samples")

  alphas   <- as.numeric(sapply(results, extractPvalues))
  zlists    <- lapply(results, function(x) x$zlist)
  momics   <- lapply(results, function(x) x$moView)
  analysis <- "twoClass"

  names(alphas) <- NULL
  new("MultiOmicsModules", 
      alphas=alphas, 
      zlists=zlists, 
      modulesView=momics, 
      modules=cliques, 
      analysis=analysis,
      multiOmicObj=deparse(substitute(omicsObj)), 
      title=pathName
      )
}
# NO graphNEL slot