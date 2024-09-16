setClassUnion("characterOrNULL", c("character", "NULL"))


#' Omics class initializer function
#' 
#' makeOmics creates the `Omics` object necessary to perform most of the 
#' analyses of this package. It contains all the omics data in the format of a 
#' `ExperimentList`, the clinical data, and all the information necessary for 
#' the dimensionality reduction step.
#' 
#' @param experiments A \code{list} or \link{ExperimentList} of all
#' combined experiments
#' @param colData A \code{\linkS4class{DataFrame}} or \code{data.frame} of
#' characteristics for all biological units
#' @param sampleMap A \code{DataFrame} or \code{data.frame} of assay names,
#' sample identifiers, and colname samples
#' @param metadata An optional argument of "ANY" class (usually list) for
#' content describing the experiments
#' @param drops A \code{list} of unmatched information
#' (included after subsetting)
#' @param modelInfo A list with length equal to length(data) that are modelInfo 
#' to process each dataset
#' @param specificArgs a list with length equal to length(data) to set 
#' additional parameters specific of the modelInfo
#' 
#' @return an `Omics` class object
#' 
#' @examples
#' data(ovarianDataset)
#' 
#' myColData <- data.frame(status = sample(c(0,1), 50, replace = TRUE),
#'                         days = sample(c(0, 500), 50, replace = TRUE),
#'                         row.names = colnames(ovarianDataset$exp))
#' 
#' myOmicsObj <- makeOmics(experiments = ovarianDataset,
#'                         colData = myColData,
#'                         modelInfo = c(
#'                           "summarizeWithPca",
#'                           "summarizeInCluster",
#'                           "summarizeToNumberOfEvents",
#'                           "summarizeToNumberOfDirectionalEvents"),
#'                         specificArgs = list(
#'                           pcaArgs = list(name = "exp", shrink = "FALSE",
#'                                          method = "sparse", maxPCs = 3),
#'                           clusterArgs = list(name = "met",
#'                                              max_cluster_number = 3),
#'                           countEvent = list(name = "mut", min_prop = 0.05),
#'                           cnvAgv = list(name = "cnv", min_prop = 0.05)))
#' 
#' @importFrom S4Vectors DataFrame
#'
#' @export

makeOmics <- function(experiments = ExperimentList(),
                      colData = S4Vectors::DataFrame(),
                      sampleMap = S4Vectors::DataFrame(
                        assay = factor(),
                        primary = character(),
                        colname = character()),
                      metadata = list(),
                      drops = list(),
                      modelInfo = character(),
                      specificArgs = list())
   {
   if (missing(sampleMap)) {
      if (missing(colData)){
         MAE <- MultiAssayExperiment(experiments)
      } else {
            MAE <- MultiAssayExperiment(experiments, colData)
         }

   } else if (missing(colData)){
      colData <- S4Vectors::DataFrame(
         row.names <- unique(sampleMap[["primary"]])
         )
      MAE <- MultiAssayExperiment(experiments, colData, sampleMap)
   }

  else { MAE <- MultiAssayExperiment(experiments, colData, sampleMap) }

   if (length(MAE@ExperimentList) != length(modelInfo)){
      message(
      "Data and relative methods to analyze them must be equal in length.")
      return()
      }
   if (length(MAE@ExperimentList) != length(specificArgs)){
      message("data and specificArgs must be equal in length.")
      return()
      }

   match <- !(modelInfo %in% availableOmicMethods())
   if (any(match)) {
      message(paste(paste(modelInfo[match], collapse=", "),
                    "modelInfo not found in method. Try availableOmicMethods."))
      return()
   }

   samplesNumber <- vapply(MAE@ExperimentList, ncol, integer(1))
   if (length(unique(samplesNumber))!=1){
      message("Mismatch in sample numbers")
      return()
      }

   samples <- lapply(MAE@ExperimentList, colnames)
   if (length(samples) > 1) {
      ref <- samples[[1]]
      cmps <- vapply(seq(from=2, to=length(samples)), function(i) {
         identical(ref, samples[[i]])
      }, logical(1))
      if (!all(cmps)){
        message("Samples order mismatch")
        return()
        }}
   duplo <- vapply(MAE@ExperimentList, function(data) {
     any(duplicated(row.names(data)))}, 
     logical(1))
   if (any(duplo)) {
      message(paste0("Duplicated row.names found in omics ",
                     paste(which(duplo), collapse = ", ")))
      return()
   }

   data_validation <- lapply(MAE@ExperimentList@listData,
                             function(exp) all(exp == 0 | is.na(exp)))
   
   if (any(data_validation == TRUE)) {
     stop(paste0("Some experiment matrices contains only 0 or NA values: ",
                 paste(names(data_validation)[which(data_validation == TRUE)],
                       collapse = ", ")))
   }
   
   new("Omics",
       ExperimentList = MAE@ExperimentList,
       colData = MAE@colData,
       sampleMap = MAE@sampleMap,
       metadata = MAE@metadata,
       modelInfo = modelInfo,
       specificArgs = specificArgs)
}



#' A generic functions showing parameter associated with each omics
#' @param object an object of class `Omics`
#' 
#' @return NULL
#' 
#' @examples
#' data(multiOmics)
#' 
#' showOmics(multiOmics)
#' 
#' @export
setGeneric("showOmics", function(object) standardGeneric("showOmics"))

#' @export
#' @describeIn Omics shows model parameters
#' @param object an object of class `Omics`
#' @return NULL

setMethod("showOmics",  signature(object = "Omics"),
          function(object) {
            nm <- names(assays(object))
            if (is.null(nm))
              nm <- seq_len(length(nm))
            for (i in seq_along(nm)) {
              cat(paste0("Data \"", nm[i], "\" to be process with \"",
                         object@modelInfo[[i]],"\". "))
              if (is.null(object@specificArgs[[i]])) {
                cat("Default parameters\n")
              } else {
                cat("Arguments:\n")
                arguments <- vapply(
                  seq_along(object@specificArgs[[i]]), function(argI) {
                    values <- object@specificArgs[[i]][[argI]]
                    if (is.list(values)) {
                      values <- values[seq_len(2)]
                      values <- paste(paste(names(values),
                                            unlist(values), sep=" -> "),
                                      collapse = " ,")
                    values <- paste0(values, " ...")
                    }
                    paste(names(object@specificArgs[[i]])[argI],
                          values, sep=" :")
                    }, character(1))
                cat(paste(arguments, collapse ="\n"), "\n")
              }
              }
            })

#' A generic function showing pathway's module info
#' @param object an object of class `MultiOmicsModules`
#' 
#' @return NULL
#' 
#' @export
setGeneric("showModule", function(object)
   standardGeneric("showModule"))

#' @describeIn MultiOmicsModules shows module info
#' @param object an object of class `MultiOmicsModules`
#' @return NULL
#' @examples
#' data(multiomics)
#' data(reactSmall)
#' 
#' genesToUse <- row.names(multiOmics[[1]])
#' 
#' MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
#'   survFormula="Surv(days, status) ~", autoCompleteFormula = TRUE,
#'   useTheseGenes = genesToUse)
#'   
#' showModule(MOM_survival)
#' 
#' @export
setMethod("showModule",
          signature = "MultiOmicsModules",
          definition = function(object){
            sthis <- seq_len(min(length(object@alphas), 3))
            sthis <- order(object@alphas)[sthis]

            sigCliquesIdx <- which(object@alphas <= 0.05)

            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }

            for (i in sthis) {
              cat(paste0("Module ",i, ": pvalue ", object@alphas[i], "\n"))
              covs <- names(which(object@zlists[[i]] <=0.05))
              if (length(covs)!=0)
                cat("The following covariates are implicated:\n",
                    paste(covs, collapse=", "),"\n")
              cat("Module is composed by the followings:\n")
              cat(paste(object@modules[[i]], collapse=", "))
              cat("\n-+-\n")
            }

            if (length(sthis) < length(sigCliquesIdx)) {
              cat(paste0("There are other ",
                         length(sigCliquesIdx)-length(sthis),
                         " cliques with pvalue <= 0.05"))
            }

            invisible(NULL)
          })

#' A generic function showing pathway info
#' @param object an object of class `MultiOmicsPathway`
#' 
#' @return NULL
#' 
#' @examples
#' data(multiomics)
#' data(reactSmall)
#' 
#' genesToUse <- row.names(multiOmics[[1]])
#' 
#' MOP_survival <- multiOmicsSurvivalPathwayTest(multiOmics, reactSmall[[1]],
#'   survFormula="Surv(days, status) ~", autoCompleteFormula = TRUE,
#'   useTheseGenes = genesToUse
#'   
#' showPathway(MOP_survival)
#' 
#' @export
setGeneric("showPathway", function(object)
   standardGeneric("showPathway") )

#' @export
#' @describeIn MultiOmicsPathway shows module info
#' @param object an object of class `MultiOmicsPathway`
#' @return NULL
setMethod("showPathway",
          signature = "MultiOmicsPathway",
          definition = function(object){
            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }
            cat(paste0("Pathway overall pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })

#' A generic function to convert pathway
#' @param graph a graphNEL object
#' @param useThisGenes list of genes to be used
setGeneric("convertPathway", function(graph, useThisGenes)
   standardGeneric("convertPathway"))

#' @importFrom graph subGraph
#' @importFrom graphite pathwayGraph
#' @importClassesFrom graphite Pathway
setMethod("convertPathway",
          signature("Pathway", "characterOrNULL"),
          function(graph, useThisGenes) {
             graph <- graphite::pathwayGraph(graph)
             if (!is.null(useThisGenes)) {
                usableGenes <- intersect(useThisGenes, graph::nodes(graph))
                graph <- graph::subGraph(usableGenes, graph)
             }
             graph
          })

#' @importFrom graph subGraph
#' @importClassesFrom graph graphNEL
setMethod("convertPathway",
          signature("graphNEL", "characterOrNULL"),
          function(graph, useThisGenes) {
             if (!is.null(useThisGenes)) {
                usableGenes <- intersect(useThisGenes, graph::nodes(graph))
                graph <- graph::subGraph(usableGenes, graph)
             }
             graph
          })

#' @importFrom graph subGraph
setMethod("convertPathway",
          signature("character", "characterOrNULL"),
          function(graph, useThisGenes) {
             graph <- new("graphNEL", nodes = graph, edgeL = list())
             if (!is.null(useThisGenes)) {
                usableGenes <- intersect(useThisGenes, graph::nodes(graph))
                graph <- graph::subGraph(usableGenes, graph)
             }
             graph
          })
