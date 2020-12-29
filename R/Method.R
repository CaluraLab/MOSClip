setClassUnion("characterOrNULL", c("character", "NULL"))


setClass(Class ="MultiOmicsModules", 
         slots = c(alphas  = "numeric",
                   zlists  = "list",
                   coxObjs = "list",
                   modulesView  = "list",
                   modules     = "list",
                   formulas = "list",
                   title = "characterOrNULL"))



setClass(Class ="MultiOmicsPathway", 
         slots = c(pvalue = "numeric",
                   zlist = "numeric",
                   coxObj = "data.frame",
                   pathView = "list",
                   formula = "character",
                   title = "characterOrNULL"))


setClass(Class ="SinglePath", package = "biocmosclip",
         slots = c(global = "MultiOmicsPathway",
                   modules = "MultiOmicsModules"))

#' This class is the storage for the different omic datasets that we need to analyze.
#'
#' @slot modelInfo a list with length equal to length(data) that are modelInfo to process each dataset.
#' @slot specificArgs a list with length equal to length(data) to set additional parameters specific of the modelInfo.
#'
#' @name Omics-class
#' @rdname Omics-class
#' 
#' @export
#' @import MultiAssayExperiment
Omics <- setClass(Class = "Omics", package = "biocmosclip",
         
         slots = c(modelInfo = "list",
                   specificArgs = "list",
                   pathResult = "list"),
         contains = "MultiAssayExperiment"
)



#' @export
Omics <- function(experiments = ExperimentList(), colData = S4Vectors::DataFrame(),
               sampleMap = S4Vectors::DataFrame(
                  assay = factor(),
                  primary = character(),
                  colname = character()),
               metadata = list(),
               drops = list(),
               modelInfo= character(),
               specificArgs= list(),
               pathResult= list())
      {
         MAE <- MultiAssayExperiment(experiments)
         cat("nnnn\n")
         new("Omics",
             ExperimentList = MAE@ExperimentList,
             colData = MAE@colData,
             sampleMap = MAE@sampleMap,
             metadata = MAE@metadata,
             modelInfo = modelInfo,
             specificArgs = specificArgs,
             pathResult = pathResult)
      }
   
      
        # 
        # if (missing(modelInfo)){
        #    modelInfo <- vector("list", length(assays(.Object)))
        # }
        # if (missing(specificArgs)) {
        #    specificArgs <- vector("list", length(assays(.Object)))
        # }
        # 
        # # .Object@colData@rownames <- character()
        # # .Object@SampleMap@listData <- list(assay = factor(), primary = character(),colname = character())
        # .Object@modelInfo <- modelInfo
        # .Object@specificArgs <- specificArgs
        # .Object@pathResult <- list(a = matrix(1,3,3))
        # 
    
  


#' @export
setGeneric("showOmics", function(object) standardGeneric("showOmics"))

#' @export
setMethod("showOmics",  signature(object = "Omics"),
          function(object) {
            nm <- names(assays(object))
            if (is.null(nm))
              nm <- seq_len(length(nm))
            for (i in seq_along(nm)) {
              cat(paste0("Data \"", nm[i], "\" to be process with \"", object@modelInfo[[i]],"\". "))
              if (is.null(object@specificArgs[[i]])) {
                cat("Default parameters\n")
              } else {
                cat("Arguments:\n")
                arguments <- sapply(seq_along(object@specificArgs[[i]]), function(argI) {
                  values <- object@specificArgs[[i]][[argI]]
                  if (is.list(values)) {
                    values <- values[1:2]
                    values <- paste(paste(names(values), unlist(values), sep=" -> "), collapse = " ,")
                    values <- paste0(values, " ...")
                  }
                  paste(names(object@specificArgs[[i]])[argI], values, sep=" :")
                })
                cat(paste(arguments, collapse ="\n"), "\n")
              }
              
            }
          })

#' @export
setGeneric("check", function(object) standardGeneric("check") )

#' @export
setMethod("check",  signature(object = "Omics"),
   check_Omics <- function(object) {
      if (length(object@ExpermentList@listData) != length(object@modelInfo)){
      msg <- "Data and relative methods to analyze them must be equal in length."
      return(msg)
   }
  
   if (length(object@ExpermentList@listData) != length(object@specificArgs)){
      msg <- "data and specificArgs must be equal in length."
      return(msg)
   }
  
   match <- !(object@modelInfo %in% availableOmicMethods())
   if (any(match)) {
      msg <- paste(paste(object@modelInfo[match], collapse=", "),
                   "methods not found. Try availableOmicMethods.")
      return(msg)
   }
  
   samplesNumber <- sapply(object@ExpermentList@listData, ncol)
   if (length(unique(samplesNumber))!=1)
      return("Mismatch in sample numbers")
  
   samples <- lapply(object@ExpermentList@listData, colnames)
   if (length(samples) > 1) {
      ref <- samples[[1]]
      cmps <- sapply(seq(from=2, to=length(samples)), function(i) {
        identical(ref, samples[[i]])
      })
   }
   
   if (!all(cmps))
    return("Samples order mismatch")
  
   duplo <- sapply(object@ExpermentList@listData, function(data) {
      any(duplicated(row.names(data)))
   })
   
   if (any(duplo)) {
      return(paste0("Duplicated row.names found in omics ",
                    paste(which(duplo), collapse = ", ")))
   }
  ###### PERCHE?
  return(TRUE)
})

#' @export
setGeneric("showModule", function(object) 
   standardGeneric("showModule") )

#' @export
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
                cat("The following covariates are implicated:\n",paste(covs, collapse=", "),"\n")
              cat("Module is composed by the followings:\n")
              cat(paste(object@modules[[i]], collapse=", "))
              cat("\n-+-\n")
            }
            
            if (length(sthis) < length(sigCliquesIdx)) {
              cat(paste0("There are other ", length(sigCliquesIdx)-length(sthis), " cliques with pvalue <= 0.05"))
            }
            
            invisible(NULL)
          })

#' @export
setGeneric("showPathway", function(object) 
   standardGeneric("showPathway") )

#' @export
setMethod("showPathway",
          signature = "MultiOmicsPathway",
          definition = function(object){
            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }
            cat(paste0("Pathway overall pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })


setClassUnion("characterOrNULL", c("character", "NULL"))

setGeneric("convertPathway", function(graph, useThisGenes) standardGeneric("convertPathway"))

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