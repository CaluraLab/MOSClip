setClassUnion("characterOrNULL", c("character", "NULL"))


setClass("MultiOmicsModules", package = "MOSClip",
         slots = c(alphas  = "numeric",
                   zlists  = "list",
                   coxObjs = "list",
                   modulesView  = "list",
                   modules     = "list",
                   formulas = "list",
                   title = "characterOrNULL"))



setClass("MultiOmicsPathway", package = "MOSClip",
         slots = c(pvalue = "numeric",
                   zlist = "numeric",
                   coxObj = "data.frame",
                   pathView = "list",
                   formula = "character",
                   title = "characterOrNULL"))


setClass("SinglePath", package = "MOSClip",
         slots = c(global = "MultiOmicsPathway",
                   modules = "MultiOmicsModules"))

#'export
#'
#' This class is the storage for the different omic datasets that we need to analyze.
#'
#' @slot data a MultiAssayExperiment containing all the data.
#' @slot modelInfo a list with length equal to length(data) that are modelInfo to process each dataset.
#' @slot specificArgs a list with length equal to length(data) to set additional parameters specific of the modelInfo.
#'
#' @name Omics-class
#' @rdname Omics-class
#' 
#' @exportClass Omics
#' @import MultiAssayExperiment

Omics<-setClass("Omics", package = "MOSClip",
         contains = "MultiAssayExperiment",
         slots = c(modelInfo = "list",
                   specificArgs = "list",
                   pathResult = "list")
        )







