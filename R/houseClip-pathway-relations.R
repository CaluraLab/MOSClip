#' Download Reactome Pathway Relations
#'
#' @param url the location of the file. Can be local. If NULL pick the package reactome file.
#' @param speciesAbbr species acronim
#'
#' @return a data frame
#' @importFrom utils read.table data
#'
#' @examples
#' \dontrun{
#' url = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
#' downloadPathwayRelationFromReactome(url, speciesAbbr = "HSA")
#' }
#' @export
downloadPathwayRelationFromReactome <- function(url=NULL, speciesAbbr = "HSA") {
  if (is.null(url)) {
    url = system.file("ReactomePathwaysRelation.txt", package = "biocmosclip",
                      mustWork = TRUE)
    # ReactomePathwaysRelation <- NULL
    # data(ReactomePathwaysRelation)
    # df <- ReactomePathwaysRelation
  }
  df <- read.table(url, sep="\t", header=F, quote="\"", stringsAsFactors = F, check.names = F)
  colnames(df) <- c("parent", "child")
  df <- df[grepl(speciesAbbr, df$parent) & grepl(speciesAbbr, df$child), , drop=F]
  row.names(df) <- NULL
  df
}

#' This is a copy of the file
#' "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
#' downloaded Dec 2018.
#'
#' A brand new file can be downloaded using:
#'  downloadPathwayRelationFromReactome
#'
#' @name ReactomePathwaysRelation
#' @docType data
#' @keywords data
NULL

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
getPathFathers <- function(pathway, hierarchyGraph, ord=3, plot=F) {
  checkmate::assertClass(hierarchyGraph, "igraph")

  if (!(pathway %in% names(V(hierarchyGraph)))){
    warning(paste0("Id ", pathway, " is not in the hierarchy."))
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
#' @export
#'
id2name <- function(idList, namedVect) {
  stopifnot(!is.null(names(namedVect)))
  lapply(idList, function(x) {
    unlist(unname(namedVect[x]))
  })
}
