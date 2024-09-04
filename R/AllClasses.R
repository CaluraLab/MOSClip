setClassUnion("characterOrNULL", c("character", "NULL"))


#' Multi Omics Modules.
#'
#' This class organize the results of the Multi Omics Module Test analysis, in
#' which corresponds to one pathway decomposed into modules. 
#'
#' @slot alphas numeric vector of the pvalues of all the modules.
#' @slot zlists a list of all the zs of the covariates.
#' @slot coxObjs a list with all the data.frame used for coxph of each module.
#' @slot modulesView a list of module information: for each omic the module data, the method used and the covariate analyzed.
#' @slot modules a list woth the genes that belong to the module.
#' @slot formulas a list, for each module the character of the formula used in the coxph.
#' @slot title the name of the pathway.
#' @slot analysis The type of analysis done: survival or two-class.
#' @name MultiOmicsModules-class
#' @rdname MultiOmicsModules-class
#'
setClass(Class ="MultiOmicsModules",
         slots = c(alphas  = "numeric",
                   zlists  = "list",
                   modulesView  = "list",
                   modules     = "list",
                   analysis = "characterOrNULL",
                   multiOmicObj = "characterOrNULL",
                   title = "characterOrNULL"))

#' Multi Omics Pathway.
#'
#' This class organize the results of the Multi Omics Module Survival Test analysis.
#'
#' @slot pvalue numeric, the pvalues of the whole module.
#' @slot zlist a vector of all the zs of the covariates.
#' @slot coxObj a data.frame used for the coxph model.
#' @slot pathView a list, for each omic the pathway data, the method used and the covariate analyzed.
#' @slot title the name of the pathway.
#' @slot analysis The type of analysis done: survival or two-class.
#' @name MultiOmicsPathway-class
#' @rdname MultiOmicsPathway-class
#'
setClass(Class ="MultiOmicsPathway",
         slots = c(pvalue = "numeric",
                   zlist = "numeric",
                   pathView = "list",
                   analysis = "characterOrNULL",
                   multiOmicObj = "characterOrNULL",
                   title = "characterOrNULL"))

# setClass(Class ="SinglePath", package = "biocmosclip",
#          slots = c(global = "MultiOmicsPathway",
#                    modules = "MultiOmicsModules"))

#' Omics.
#' 
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
                            specificArgs = "list"),
                  contains = "MultiAssayExperiment"
)
