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
                  
                  slots = c(modelInfo = "character",
                            specificArgs = "list",
                            pathResult = "list"),
                  contains = "MultiAssayExperiment"
)