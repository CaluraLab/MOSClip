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

setClass("Omics", 
         slots = c(data = "MultiAssayExperiment",
                   modelInfo = "list",
                   specificArgs = "list")
        )


setMethod("initialize",
          signature=signature(.Object="Omics"),
          function(.Object, data, modelInfo, specificArgs) {
            if (missing(data)) {
              cat("missing\n")
              data <- MultiAssayExperiment()
            }
            if (missing(modelInfo))
              modelInfo <- vector("list", length(assays(data)))
            if (missing(specificArgs)) {
              specificArgs <- vector("list", length(assays(data)))
            }
            .Object@data <- data
            .Object@modelInfo <- modelInfo
            .Object@specificArgs <- specificArgs
            .Object
          })


#' Wrapper of Omics
#' @name Omics
#' @param ... insert the slots. see \code{Slots}
#' @rdname Omics-class
#' @export
Omics <- function(...) new("Omics", ...)

setMethod("show",
          signature = "Omics",
          definition = function(object) {
            nm <- names(assays(data))
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