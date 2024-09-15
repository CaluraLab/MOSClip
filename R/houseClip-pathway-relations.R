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
#' downloadPathwayRelationFromReactome()
#' 
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


#' Retrieves pathways relatives
#'
#' For internal use only. Retrieves relatives given a pathway id.
#'
#' Pathway Hierarchy is needed as igraph object.
#'
#' @param pathway a pathway id
#' @param hierarchyGraph a igraph with pathway hierarchy
#' @param ord how far you need to go backward
#' @param plot plot relatives. For checking purpose
#'
#' @return a character vector with the relatives
#'
#' @importFrom checkmate assertClass
#' @importFrom igraph V V<- as_ids make_ego_graph ego distances
#' @export
#'
getPathFathers <- function(pathway, hierarchyGraph, ord=3, plot=FALSE) {
  checkmate::assertClass(hierarchyGraph, "igraph")
  
  if (!(pathway %in% names(V(hierarchyGraph)))){
    warning("Id ", pathway, " is not in the hierarchy.")
    return(pathway)
  }
  
  mm <- make_ego_graph(hierarchyGraph, ord, nodes = pathway, mode="in")
  mmlist <- ego(hierarchyGraph, ord, nodes = pathway, mode="in")
  
  if (plot)
    plot(mm[[1]])
  
  chain <- as_ids(mmlist[[1]])
  parents <- chain[-1]
  if (length(parents)==0)
    return(chain)
  
  dis <- distances(mm[[1]])
  idx <- which.max(dis[pathway, ])
  return(colnames(dis)[idx])
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
