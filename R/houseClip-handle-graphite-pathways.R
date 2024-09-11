#' Map Pathways ID from Graphite
#'
#' For internal use only. Retrieve pathway id and names from Pathways object.
#'
#' @param pathways a PathwayList object
#' @param pathwayNames in not NULL, a subset of pathway to extract
#'
#' @return a data frame, id and pathway name
#' 
#' @importFrom checkmate assertClass
#'
mapPathwaysIDfromGraphite <- function(pathways, pathwayNames=NULL) {
  assertClass(pathways, "PathwayList")
  if (is.null(pathwayNames)) {
    pathwayNames <- names(pathways)
  }
  l <- lapply(pathwayNames, function(p) {
    if (!(p %in% names(pathways))) {
      warning("No id found for ", p )
      return(NULL)
    }
    data.frame(id=pathways[[p]]@id, pname=p, stringsAsFactors = FALSE)
  })
  do.call(rbind, l)
}
