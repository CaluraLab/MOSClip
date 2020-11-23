setMethod("initialize",
          signature=signature(.Object="Omics"),
          
          function(.Object, modelInfo, specificArgs,pathResult) {
            .Object <- callNextMethod()
            if (missing(modelInfo))
              modelInfo <- vector("list", length(assays(.Object)))
            if (missing(specificArgs)) {
              specificArgs <- vector("list", length(assays(.Object)))
            }
            .Object@modelInfo <- modelInfo
            .Object@specificArgs <- specificArgs
            .Object@pathResult <- list()
            .Object
          }
)

#'Wrapper of Omics
#' @name Omics
#' @param ... insert the slots. see \code{Slots}
#' @rdname Omics-class
#' @export
# Omics <- function(...) new("Omics", ...)


setMethod("show", signature = "Omics",
          definition = function(object) {
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

setMethod("check", signature = "Omics",
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







setMethod("show",
          signature = "MultiOmicsModules",
          definition = function(object) {
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


setMethod("show",
          signature = "MultiOmicsPathway",
          definition = function(object) {
            if (!is.null(object@title)) {
              cat("\"",object@title, "\"\n", sep = "")
            }
            cat(paste0("Pathway overall pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })