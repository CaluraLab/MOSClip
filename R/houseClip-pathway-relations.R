#' Download Reactome Pathway Relations
#' 
#' Download Pathway Relations from Reactome. The file is retrieved from the 
#' [url](https://reactome.org/download/current/ReactomePathwaysRelation.txt)
#'
#' @param url the location of the file. Can be local. If NULL pick the package 
#' reactome file.
#' @param speciesAbbr species acronim
#'
#' @return A data frame with 2 columns: 
#'  \item{parent}{The Reactome pathway ID of the parent pathway.}
#'  \item{child}{The Reactome pathway ID of the child pathway.}
#' 
#' 
#' @importFrom utils read.table data
#'
#' @examples
#' \donttest{
#' url = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
#' downloadPathwayRelationFromReactome(url, speciesAbbr = "HSA")
#' }
#' @export
downloadPathwayRelationFromReactome <- function(url=NULL, speciesAbbr = "HSA") {
  if (is.null(url)) {
    url <- system.file("extdata", "ReactomePathwaysRelation.txt", 
                       package = "biocmosclip", mustWork = TRUE)
  }
  df <- read.table(url, sep="\t", header=FALSE, quote="\"", 
                   stringsAsFactors=FALSE, check.names = FALSE)
  colnames(df) <- c("parent", "child")
  df <- df[grepl(speciesAbbr, df$parent) & grepl(speciesAbbr, df$child), ,
           drop=FALSE]
  row.names(df) <- NULL
  df
}

#' Convert id to pathway name
#'
#' For internal use only. Retrieves name from pathway id.
#'
#' You must provide a namedVect to be used as translator.
#'
#' @param idList a list of pathway id
#' @param namedVect a named vector
#'
#' @return a character vector with the names
#'
id2name <- function(idList, namedVect) {
  stopifnot(!is.null(names(namedVect)))
  lapply(idList, function(x) {
    unlist(unname(namedVect[x]))
  })
}
