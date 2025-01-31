#' Check formulas
#' 
#' Check input baseFormula and nullModel formula are correct.
#' 
#' @inheritParams multiOmicsTwoClassPathwayTest
#' 
#' @return NULL. Stop function if one of the formulas is not correct.
#' 

checkFormula <- function(baseFormula, nullModel, classAnnot, 
                         autoCompleteFormula, omicsObj, paired, patientCol){
  
  baseFormula_input <- strsplit(baseFormula, "\\s*~\\s*")[[1]]
  if (!(baseFormula_input[1] %in% colnames(classAnnot))) {
    stop("Invalid base formula. Class column not found in classAnnot")
  } else if (!grepl("~", baseFormula)) {
    stop("Invalid base formula. Formula should be written as: 'classes ~': ",
         "classes is the name of the class column in classAnnot")
  }
  
  nullModel_input <- strsplit(nullModel, "\\s*~\\s*")[[1]]
  if (!(nullModel_input[1] %in% colnames(classAnnot))) {
    stop("Invalid null formula. Class column not found in classAnnot")
  } else if (length(nullModel_input) == 1) {
    stop("Null formula is too short. ",
         "Formula should be written as 'classes ~ 1': ",
         "classes in the name of the class column in classAnnot")
  } else if (!grepl("~", nullModel)) {
    stop("Invalid null formula. ",
         "Formula should be written as: 'classes ~ 1': ",
         "classes is the name of the class column in classAnnot")
  }
  
  if (autoCompleteFormula == TRUE & length(baseFormula_input) != 1){
    stop("With autoCompleteFormula = TRUE baseFormula should be written as: ",
         "'classes ~'. Covariates will be automatically added, including ",
         "sample information for paired test.")
  }
  
  if (paired == TRUE){
    if (!(patientCol %in% colnames(omicsObj@colData))){
      stop(patientCol, " column for paired test not found in colData.")
    }
  }
}


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
#' written as 'classes ~ ', while 'classes' being the column name for the class
#' labels
#' @param autoCompleteFormula a logical value, default TRUE. 
#' If TRUE, it autocompletes the formula used to fit generalized linear 
#' models function  using all the available covariates (omics)
#' @param useTheseGenes (optional) vector of specific genes to be used
#' @param nullModel the null model formula. It should be written the same as
#' the baseFormula, followed by ' 1'. (e.g. 'classes ~ 1')
#' @param pathName (optional) title of the pathway. If NULL, `graph@title` is
#' used as title
#' @param paired whether samples are paired or not. Default FALSE. If TRUE, 
#' a specific column should be present in `colData` specifying paired samples
#' @param patientCol name of the column to be searched in `colData` specifying 
#' paired samples. Default "patient"
#' 
#' @examples
#' data("multiOmics")
#' data("reactSmall")
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' classAnnot <- data.frame(
#'     "treatment" = c(rep("A", 25), rep("B", 25)),
#'     row.names = colnames(multiOmics[[1]])
#' )
#'
#' MOP_twoClasses <- multiOmicsTwoClassPathwayTest(
#'     multiOmics, reactSmall[[1]], classAnnot,
#'     baseFormula = "treatment ~ ", nullModel = "treatment ~ 1",
#'     useTheseGenes = genesToUse
#' )
#'
#' @return `MultiOmicsPathway` object
#'
#' @importFrom graph nodes
#' @importFrom methods new is
#'
#' @export

multiOmicsTwoClassPathwayTest <- function(
    omicsObj, graph, classAnnot, baseFormula = "classes ~ ",
    autoCompleteFormula = TRUE, useTheseGenes = NULL, nullModel = "classes ~ 1",
    pathName = NULL, paired = FALSE, patientCol = "patient") {
    
    checkFormula(baseFormula, nullModel, classAnnot, autoCompleteFormula, 
                 omicsObj, paired, patientCol)

    if (nrow(classAnnot) != nrow(omicsObj@colData)) {
        stop("Mismatch in the number of samples")
    }
    if (is.null(pathName) && is(graph, "Pathway")) {
        pathName <- graph@title
    }

    graph <- convertPathway(graph, useTheseGenes)
    genesToUse <- graph::nodes(graph)
    if (length(genesToUse) == 0) {
        stop("There is no nodes on the graph.")
    }

    moduleView <- lapply(seq_along(omicsObj@ExperimentList@listData), 
                         function(i) {
        test <- get(omicsObj@modelInfo[i])
        specificArgs <- omicsObj@specificArgs[[i]]

        cliques <- NULL
        if (omicsObj@modelInfo[i] == "summarizeWithPca") {
            genesToUse <- intersect(
                row.names(omicsObj@ExperimentList@listData[[i]]),
                genesToUse
            )
            if (identical(genesToUse, character(0))) {
                stop("Genes not found in some experiments of the omics object")
            }
            graph <- graph::subGraph(genesToUse, graph)
            cliques <- extractCliquesFromDag(graph)
        }
        args <- list(
            data = omicsObj@ExperimentList@listData[[i]], features = genesToUse,
            cliques = cliques
        )
        if (!is.null(specificArgs)) {
            args <- c(args, specificArgs)
        }
        do.call(test, args)
    })

    moduleView <- moduleView[!vapply(moduleView, is.null, logical(1))]
    covariates <- lapply(moduleView, function(mo) {
        mo$x
    })
    covariates <- do.call(cbind, covariates)

    if (is.null(covariates)) {
        return(NULL)
    }

    if (!identical(row.names(classAnnot), row.names(covariates))) {
        if (all(row.names(classAnnot) %in% row.names(covariates))) {
            res <- resolveAndOrder(list(classAnnot = classAnnot, 
                                        covariates = covariates))
            classAnnot <- res$classAnnot
            covariates <- res$covariates
        } else {
            stop("Mismatch in covariates and annotations row names.")
        }
    }

    dataTest <- data.frame(classAnnot, covariates)

    nullModelFormula <- nullModel
    if (paired == TRUE) {
      if (grepl(patientCol, nullModelFormula) == FALSE){
        nullModelFormula <- paste0(c(nullModel, patientCol), collapse = "+")
      }
    }

    dependentVar <- all.vars(as.formula(nullModelFormula))[1]
    if (!(dependentVar %in% colnames(dataTest))) {
        stop("Data does not contain one of the model dependent variables")
    }
    
    dataTest <- dataTest[, c(dependentVar, colnames(covariates))]

    twoClasses <- unique(dataTest[, dependentVar])
    if (length(twoClasses) != 2) {
        stop("Classes should be only two. ", 
             "Check your dependent variables columns")
    }

    if (!(all(twoClasses %in% c(0, 1)) & is.numeric(twoClasses))) {
        dataTest[dataTest[, dependentVar] == twoClasses[1], dependentVar] <- 0
        dataTest[dataTest[, dependentVar] == twoClasses[2], dependentVar] <- 1
        dataTest[, dependentVar] <- as.numeric(dataTest[, dependentVar])
    }
    
    fullModelFormula <- baseFormula
    if (autoCompleteFormula) {
        if (paired == TRUE){
          patient <- omicsObj@colData[, patientCol]
          fullModelFormula <-  paste0(baseFormula, 
                                      paste(c(colnames(covariates), patientCol),
                                            collapse = "+"))
          dataTest <- cbind(patient, dataTest)
        }
        else{
          fullModelFormula <- paste0(baseFormula, paste(colnames(covariates), 
                                                        collapse = "+"))
    }}

    res <- glmTest(dataTest, fullModelFormula, nullModelFormula)

    new("MultiOmicsPathway",
        pvalue = res$pvalue, zlist = res$zlist, pathView = moduleView,
        analysis = "twoClass", multiOmicObj = deparse(substitute(omicsObj)), 
        title = pathName
    )
}

### ****************************fine Pathway***************************


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
#' data("multiOmics")
#' data("reactSmall")
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' classAnnot <- data.frame(
#'     "treatment" = c(rep("A", 25), rep("B", 25)),
#'     row.names = colnames(multiOmics[[1]])
#' )
#'
#' MOM_twoclasses <- multiOmicsTwoClassModuleTest(
#'     multiOmics, reactSmall[[1]], classAnnot,
#'     baseFormula = "treatment ~ ", nullModel = "treatment ~ 1",
#'     useTheseGenes = genesToUse
#' )
#'
#' @importFrom graph nodes
#' @importFrom methods new is
#' @export

multiOmicsTwoClassModuleTest <- function(
    omicsObj, graph, classAnnot, baseFormula = "classes ~",
    autoCompleteFormula = TRUE, useTheseGenes = NULL, nullModel = "classes ~ 1",
    pathName = NULL, paired = FALSE, patientCol = "patient") {
    if (is(graph, "character")) {
        stop("graph argument should be a graphNEL object, not a gene list.")
    }
  
    checkFormula(baseFormula, nullModel, classAnnot, autoCompleteFormula, 
                 omicsObj, paired, patientCol)

    if (nrow(classAnnot) != nrow(omicsObj@colData)) {
        stop("Mismatch in the number of samples")
    }
    if (is.null(pathName) & is(graph, "Pathway")) {
        pathName <- graph@title
    }

    graph <- convertPathway(graph, useTheseGenes)
    genes <- graph::nodes(graph)
    if (length(genes) == 0) {
        stop("There is no intersection between expression feature names and ", 
             "the node names in the graph.")
    }

    cliques <- extractCliquesFromDag(graph)

    results <- lapply(cliques, MOMglmTest,
        omicsObj = omicsObj, classAnnot = classAnnot,
        baseFormula = baseFormula, autoCompleteFormula = autoCompleteFormula, 
        nullModel = nullModel, paired = paired, patientCol = patientCol
    )

    alphas <- as.numeric(vapply(results, extractPvalues, numeric(1)))
    zlists <- lapply(results, function(x) x$zlist)
    momics <- lapply(results, function(x) x$moView)
    analysis <- "twoClass"

    names(alphas) <- NULL
    new("MultiOmicsModules",
        alphas = alphas, zlists = zlists, modulesView = momics,
        modules = cliques, analysis = analysis, 
        multiOmicObj = deparse(substitute(omicsObj)),
        title = pathName
    )
}
