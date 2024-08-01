setClassUnion("characterOrNULL", c("character", "NULL"))


#' Omics class initializer function
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
#' @param modelInfo a list with length equal to length(data) that are modelInfo to process each dataset.
#' @param specificArgs a list with length equal to length(data) to set additional parameters specific of the modelInfo.
#' @param pathResult a list containing pathways result
#' @export
#' @importFrom S4Vectors DataFrame
makeOmics <- function(experiments = ExperimentList(),
                  colData = S4Vectors::DataFrame(),
                  sampleMap = S4Vectors::DataFrame(
                     assay = factor(),
                     primary = character(),
                     colname = character()),
                  metadata = list(),
                  drops = list(),
                  modelInfo= character(),
                  specificArgs= list())
   {
   if (missing(sampleMap)) {
      if (missing(colData)){
         MAE <- MultiAssayExperiment(experiments)
      } else {
            MAE <- MultiAssayExperiment(experiments,colData)
         }

   } else if (missing(colData)){
      colData <- S4Vectors::DataFrame(
         row.names = unique(sampleMap[["primary"]])
         )
      MAE <- MultiAssayExperiment(experiments, colData, sampleMap)
   }
  
  else { MAE <- MultiAssayExperiment(experiments, colData, sampleMap) }

   if (length(MAE@ExperimentList) != length(modelInfo)){
      message("Data and relative methods to analyze them must be equal in 
              length.")
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

   samplesNumber <- sapply(MAE@ExperimentList, ncol)
   if (length(unique(samplesNumber))!=1){
      message("Mismatch in sample numbers")
      return()
      }

   samples <- lapply(MAE@ExperimentList, colnames)
   if (length(samples) > 1) {
      ref <- samples[[1]]
      cmps <- sapply(seq(from=2, to=length(samples)), function(i) {
         identical(ref, samples[[i]])
      })
   }
   if (!all(cmps)){
      message("Samples order mismatch")
      return()
   }
   duplo <- sapply(MAE@ExperimentList, function(data) {
      any(duplicated(row.names(data)))
   })
   if (any(duplo)) {
      message(paste0("Duplicated row.names found in omics ",
                     paste(which(duplo), collapse = ", ")))
      return()
   }

   new("Omics",
       ExperimentList = MAE@ExperimentList,
       colData = MAE@colData,
       sampleMap = MAE@sampleMap,
       metadata = MAE@metadata,
       modelInfo = modelInfo,
       specificArgs = specificArgs)
}



#' A generic functions showing parameter asociated with each omics
#' @param object an object of class Omics
#' @export
setGeneric("showOmics", function(object) standardGeneric("showOmics"))

#' @export
#' @describeIn Omics shows model parameters
#' @param object an object of class Omics

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
                arguments <- sapply(
                  seq_along(object@specificArgs[[i]]), function(argI) {
                    values <- object@specificArgs[[i]][[argI]]
                    if (is.list(values)) {
                      values <- values[1:2]
                      values <- paste(paste(names(values),
                                            unlist(values), sep=" -> "),
                                      collapse = " ,")
                    values <- paste0(values, " ...")
                    }
                    paste(names(object@specificArgs[[i]])[argI],
                          values, sep=" :")
                    })
                cat(paste(arguments, collapse ="\n"), "\n")
              }
              }
            })

#' A generic function showing pathway's module info
#' @param object an object of class MultiOmicsModules
#' @export
setGeneric("showModule", function(object)
   standardGeneric("showModule") )

#' @export
#' @describeIn MultiOmicsModules shows module info
#' @param object an object of class MultiOmicsModules
setMethod("showModule",
          signature = "MultiOmicsModules",
          definition = function(object){
            sthis <- seq_len(min(length(object@alphas), 3))
            sthis <- order(object@alphas)[sthis]

            sigCliquesIdx = which(object@alphas <= 0.05)

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
#' @param object an object of class MultiOmicsPathway
#' @export
setGeneric("showPathway", function(object)
   standardGeneric("showPathway") )

#' @export
#' @describeIn MultiOmicsPathway shows module info
#' @param object an object of class MultiOmicsPathway
setMethod("showPathway",
          signature = "MultiOmicsPathway",
          definition = function(object){
            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }
            cat(paste0("Pathway overall pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })



#' @export
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
