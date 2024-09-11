#' Compute Multi Omics Two-Class in Pathways
#'
#' Performs topological two-class analysis using an `Omics` object.
#' 
#' @param omicsObj object of class `Omics`
#' @param graph a pathway as a `graphNEL` object.
#' @param classAnnot a `data.frame` with the class annotation. It is necessary
#' at least a column with the classes labels, and the row.names as the samples 
#' labels
#' @param baseFormula model formula to be used for the test. It should be
#' written as "classes ~ ", while "classes" being the column name for the class
#' labels
#' @param autoCompleteFormula a logical value. If TRUE. It autocompletes the
#' formula used to fit generalized lienar models function  using all the
#' available covariates (omics)
#' @param useTheseGenes (optional) vector of specific genes to be used
#' @param nullModel the null model formula. It should be written the same as 
#' the baseFormula, followed by " 1". (e.g. "classes ~ 1")
#' @param pathName (optional) title of the pathway. If NULL, `graph@title` is
#' used as title
#' @examples 
#' data("multiOmicsTwoClass")
#' data("reactSmall")
#' 
#' genesToUse <- row.names(multiOmicsTwoClass[[1]])
#' 
#' classAnnot <- data.frame("treatment" = c(rep("A", 25), rep("B", 25)),
#'                          row.names = colnames(multiOmicsTwoClass[[1]]))
#' 
#' MOP_twoClasses <- multiOmicsTwoClassesPathwayTest(
#'   multiOmicsTwoClass, reactSmall[[1]], classAnnot,
#'   baseFormula = "treatment ~ ", nullModel = "treatment ~ 1",
#'   useTheseGenes = genesToUse)
#'   
#' @return `MultiOmicsPathway` object
#'
#' @importFrom graph nodes
#' @importFrom methods new is
#' 
#' @export

multiOmicsTwoClassPathwayTest <- function(omicsObj, graph, classAnnot,
                                          baseFormula = "classes ~ ",
                                          autoCompleteFormula = TRUE,
                                          useTheseGenes = NULL,
                                          nullModel = "classes ~ 1",
                                          pathName = NULL) {

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
  
  if (nrow(classAnnot) != nrow(omicsObj@colData))
    stop("Mismatch in the number of samples")
  if (is.null(pathName) && is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useTheseGenes)
  genesToUse <- graph::nodes(graph)
  if (length(genesToUse) == 0)
    stop("There is no nodes on the graph.")

  moduleView <- lapply(seq_along(omicsObj@ExperimentList@listData), function(i) 
    {
    test <- get(omicsObj@modelInfo[i]) 
    specificArgs <- omicsObj@specificArgs[[i]]

    cliques <- NULL
    if (omicsObj@modelInfo[i] == "summarizeWithPca") {
      genesToUse <- intersect(row.names(omicsObj@ExperimentList@listData[[i]]),
                              genesToUse)
      if (identical(genesToUse, character(0))) {
        stop("Genes not found in some experiments of the omics object")
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
  
  moduleView <- moduleView[!vapply(moduleView, is.null, logical(1))]
  covariates <- lapply(moduleView, function(mo) {mo$x})
  covariates <- do.call(cbind, covariates)

  if (is.null(covariates))
    return(NULL)

  if (!identical(row.names(classAnnot), row.names(covariates))) {
    if (all(row.names(classAnnot) %in% row.names(covariates))) {
      res <- resolveAndOrder(list(classAnnot = classAnnot, 
                                  covariates = covariates))
      classAnnot <- res$classAnnot
      covariates <- res$covariates 
    } 
    else {stop("Mismatch in covariates and annotations row names.") }}
  
  dataTest <- data.frame(classAnnot, covariates)
  
  nullModelFormula <- nullModel
  
  dependentVar <- all.vars(as.formula(nullModelFormula))[1]
  if(!(dependentVar %in% colnames(dataTest)))
    stop("Data does not contain one of the model dependent variables")

  twoClasses <- unique(dataTest[,dependentVar])
  if(length(twoClasses) != 2)
    stop("Classes should be only two. Check your dependent variables columns")

  if(!all(twoClasses == c(0,1))) {
    dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
    dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
    dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
  }

  fullModelFormula <- baseFormula
  if (autoCompleteFormula) 
    fullModelFormula <- paste0(baseFormula,
                             paste(colnames(covariates), collapse ="+"))

  res <- glmTest(dataTest, fullModelFormula, nullModelFormula)

  new("MultiOmicsPathway", 
      pvalue=res$pvalue, 
      zlist=res$zlist, 
      pathView=moduleView, 
      analysis="twoClass",
      multiOmicObj=deparse(substitute(omicsObj)),
      title=pathName)
}

###****************************fine Pathway***************************


#' Computes Multi Omics Two-Class in Pathway Modules
#'
#' Performs topological two-class analysis using an `Omics` object. It
#' decomposes graphs (pathways) into modules.
#'
#' @inheritParams multiOmicsTwoClassPathwayTest
#'
#' @return `MultiOmicsModule` object
#'
#' @examples
#' data("multiOmicsTwoClass")
#' data("reactSmall")
#' 
#' genesToUse <- row.names(multiOmicsTwoClass[[1]])
#' 
#' classAnnot <- data.frame("treatment" = c(rep("A", 25), rep("B", 25)),
#'                          row.names = colnames(multiOmicsTwoClass[[1]]))
#'                          
#' MOM_twoclasses <- multiOmicsTwoClassesModuleTest(
#'   multiOmicsTwoClass, reactSmall[[1]], classAnnot,
#'   baseFormula = "treatment ~ ", nullModel = "treatment ~ 1",
#'   useTheseGenes = genesToUse)
#'   
#' @importFrom graph nodes
#' @importFrom methods new is
#' @export

multiOmicsTwoClassModuleTest <- function(omicsObj, graph, classAnnot,
                                         baseFormula = "classes ~",
                                         autoCompleteFormula=TRUE,
                                         useTheseGenes=NULL,
                                         nullModel = "classes ~ 1",
                                         pathName=NULL) {

  if (is(graph, "character"))
    stop("graph argument should be a graphNEL object, not a gene list.")
  
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
    stop("Null formula should be written as 'classes ~ 1")
  } else if (nullModel_input[2] != "~" | nullModel_input[3] != "1") {
    stop("Invalid null formula. Formula should be written as: 'classes ~'")
  }
  
  if (nrow(classAnnot) != nrow(omicsObj@colData))
    stop("Mismatch in the number of samples")
  if (is.null(pathName) & is(graph, "Pathway"))
    pathName <- graph@title

  graph <- convertPathway(graph, useTheseGenes)
  genes <- graph::nodes(graph)
  if (length(genes) == 0)
    stop("There is no intersection between expression feature names and ",
         "the node names in the graph.")
  
  cliques <- extractCliquesFromDag(graph)

  results <- lapply(cliques, MOMglmTest, omicsObj=omicsObj,
                    classAnnot=classAnnot,
                    baseFormula=baseFormula,
                    autoCompleteFormula=autoCompleteFormula,
                    nullModel=nullModel)

  alphas   <- as.numeric(vapply(results, extractPvalues, numeric(1)))
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
