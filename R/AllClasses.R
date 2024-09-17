setClassUnion("characterOrNULL", c("character", "NULL"))


#' Multi Omics Modules.
#'
#' This class organizes the results of the Multi Omics Module Test analysis, in
#' which corresponds to one pathway decomposed into modules. 
#'
#' @slot alphas a numeric vector of the pvalues of all the modules.
#' @slot zlists a list of numeric vectors with the zs of the covariates 
#' for each module.
#' @slot modulesView a list of module information: for each omic, 
#' the name of the omic, the genes used, the method, the name of the 
#' covariates analyzed and other specific information based on the omic.
#' @slot modules a list with the genes that belong to the module.
#' @slot title the name of the pathway.
#' @slot analysis the type of analysis done: survival or two-class.
#' @name MultiOmicsModules-class
#' @rdname MultiOmicsModules-class
#'
setClass(
    Class = "MultiOmicsModules", slots = c(
        alphas = "numeric", zlists = "list", modulesView = "list",
        modules = "list", analysis = "characterOrNULL",
        multiOmicObj = "characterOrNULL", title = "characterOrNULL"
    )
)

#' Multi Omics Pathway.
#'
#' This class organize the results of the Multi-Omics Pathway Survival Test 
#' analysis.
#'
#' @slot pvalue a numeric vector of the pvalues of the pathways.
#' @slot zlist a numeric vector with the zs of all the covariates. 
#' @slot pathView a list of pathway information: for each omic, 
#' the name of the omic, the genes used, the method, the name of the 
#' covariates analyzed and other specific information based on the omic.
#' @slot title the name of the pathway.
#' @slot analysis the type of analysis done: survival or two-class.
#' @name MultiOmicsPathway-class
#' @rdname MultiOmicsPathway-class
#'
setClass(
    Class = "MultiOmicsPathway", slots = c(
        pvalue = "numeric", zlist = "numeric", pathView = "list",
        analysis = "characterOrNULL", multiOmicObj = "characterOrNULL",
        title = "characterOrNULL"
    )
)


#' Omics.
#' 
#' This class is the storage for the different omic datasets that we need to 
#' analyze. It is based on `MultiAssayExperiment`.
#'
#' @slot modelInfo a list with length equal to length(data) that are modelInfo 
#' to process each dataset.
#' @slot specificArgs a list with length equal to length(data) to set 
#' additional parameters specific of the modelInfo.
#'
#' @name Omics-class
#' @rdname Omics-class
#'
#' @import MultiAssayExperiment
Omics <- setClass(
    Class = "Omics", package = "MOSClip",
    slots = c(modelInfo = "character", specificArgs = "list"),
    contains = "MultiAssayExperiment"
)
