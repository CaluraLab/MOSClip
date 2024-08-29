#' Summarize Omics Covaraites By Min Pvalue
#'
#' For internal use only. for each line extrac 'col' and get the minimum.
#'
#' @param col columns to extract from the line
#' @param mat the matrix to be summarized (were to extract lines and 'col')
#'
#' @return a summarized version of the matrix.
#'
#' @examples
#'   summarizeOmicsResByMinPvalue(2:3, mat=matrix(c(1,2,4,1,2,5), nrow=2))
#'
#' @export
#'
summarizeOmicsResByMinPvalue <- function(col, mat) {
  apply(as.matrix(mat[, col, drop=F]), 1, minOrNA)
}


#' Minimum or NA
#'
#' For internal use only. Get back minimum or NA.
#'
#' @param x a numeric
#'
#' @return a numeric. The minimum or NA
#'
#' @examples
#'   minOrNA(c(1,5,0.1,NA))
#'   minOrNA(c(NA,NA,NA))
#'
#' @export
#'
minOrNA <- function(x) {
  if (all(is.na(x)))
    return(NA)
  min(x, na.rm=T)
}


#' Compute frequences from list
#'
#' For internal use only. Compute frequences
#'
#' @param elementsIntersections a names list
#'
#' @return a 'data.frame' with the frequences
#'
#' @export
#'
computeFreqs <- function(elementsIntersections) {
  freques <- lapply(names(elementsIntersections), function(x) {
    counts <- table(elementsIntersections[[x]])
    data.frame(category=names(counts), frequencies=as.numeric(counts),
               class=x, stringsAsFactors = F)
  })
  do.call(rbind, freques)
}


#' Plot frequences
#'
#' For internal use only. Plot supertest frequences
#'
#' @param frequencies a data.frame created from 'computeFreqs'
#' @param manualColors optional vector of colors
#' @param minSize the minimal fontsize. Maximal frequencies will be added for each class
#' @param maxSize the maximal fontsize dimension, all values above are clipped
#' @param width the number of character to wrap the labels
#' @param relMagnificationOfLegend the relative magnification of the text of the legend
#' @param lineSize the thickness of the lines
#'
#' @return NULL
#' @examples
#' df <- data.frame(category=c("talk", "too","mutch", "dear"),
#'   frequencies=c(1,2,1,3),
#'   class=rep("Mut",4), stringsAsFactors = FALSE)
#' plotFrequencies(df)
#'
#' @importFrom ggplot2 ggplot aes element_text coord_polar geom_point
#'   geom_polygon geom_path theme theme_bw labs element_blank
#'   ggplot_gtable ggplot_build rel element_line
#' @importFrom grid grid.draw grid.newpage
#' @importFrom stats reorder
#' @export
#'
plotFrequencies <- function(frequencies, manualColors=NULL, minSize=4,
                            maxSize=20, width=20, relMagnificationOfLegend=0.5,
                            lineSize=1){
  category <- NULL
  if (!all(colnames(frequencies) %in% c("category", "frequencies", "class")))
    stop("Frequences dataframe must contain columns category, frequencies and 
         class")
  if (!is.null(manualColors)) {
    g <- ggplot2::ggplot(frequencies, aes(y = frequencies,
                                          x = factor(category),
                                          group = class, colour=class))
    g <- g + ggplot2::scale_colour_manual(values=manualColors)
  } else {
    g <- ggplot2::ggplot(frequencies, aes(y = frequencies,
                                          x = factor(category),
                                          group = class, colour = class))
  }
  size <- tapply(seq_along(frequencies$frequencies),
                 factor(frequencies$category), function(idx) 
                   max(frequencies$frequencies[idx]))
  size <- as.numeric(size)+minSize
  size[size > maxSize] <- maxSize
  g <- g + ggplot2::coord_polar() +
    ggplot2::geom_point(stat='identity') +
    ggplot2::geom_polygon(fill=NA)+
    ggplot2::geom_path(linewidth=lineSize) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = element_blank(),
                   axis.line.x = element_blank(),
                   axis.line.y = element_blank()) +
    ggplot2::theme(panel.grid=ggplot2::element_line(linewidth = lineSize*0.5),
                   axis.text.x=element_text(size=size, colour="black"),
                   axis.text=element_text(colour="black")) +
    ggplot2::scale_x_discrete(labels=function(x) lapply(
      strwrap(x, width = width, simplify = FALSE), paste, collapse="\n"))
  g <- g + theme(legend.position.inside = c(1,1), legend.justification=c(0, 1),
                 legend.text=element_text(
                   size=ggplot2::rel(relMagnificationOfLegend)))
  grid::grid.newpage()
  gt <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(g))
  gt$layout$clip <- "off"
  grid::grid.draw(gt)
}
