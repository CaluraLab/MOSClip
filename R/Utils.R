conversionToSymbols <- function(idsGraphiteStyle, orgDbi="org.Hs.eg.db") {
  if (!requireNamespace(orgDbi))
    return(idsGraphiteStyle)
  
  if(is.null(idsGraphiteStyle))
    return(NULL)
  
  if (!(length(grep(":", idsGraphiteStyle)) == length(idsGraphiteStyle) &
        length(unique(do.call(rbind,(strsplit(idsGraphiteStyle, ":")))
                      [,1])) == 1 ))
    return(idsGraphiteStyle)
  
  typeId <- unique(do.call(rbind,(strsplit(idsGraphiteStyle, ":")))[,1])
  originals <- gsub(paste0(typeId,":"), "", idsGraphiteStyle)
  symbols <- select(get(orgDbi), keys=originals,
                    columns = c("SYMBOL"), keytype=typeId)$SYMBOL
  symbols[is.na(symbols)] <- originals[is.na(symbols)]
  as.character(symbols)
}


formatAnnotations <- function(listOfMostlyInvolvedGenesInOmics, sortBy) {
  involved=listOfMostlyInvolvedGenesInOmics
  samplesList <- row.names(involved[[1]]$discrete)
  annotationFull <- lapply(seq_along(involved), function(i) {
    covNames <- involved[[i]]$covsConsidered
    annotations <- involved[[i]]$discrete[samplesList, covNames, drop=F]
    annotations
  })
  annotationFull <- do.call(cbind, annotationFull)
  if (is.null(sortBy)) {
    annotationFull <- annotationFull[
      order(annotationFull[, 1],annotationFull[, ncol(annotationFull)]), ,
      drop=F]
  } else {
    ord <- getMultiColOrder(annotationFull, sortBy)
    annotationFull <- annotationFull[ord, , drop=F]
  }
  annotationFull
}

sortAnnotations <- function(annotations, sortBy) {
  if (is.null(sortBy))
    return(annotations)
  
  missing <- setdiff(sortBy, colnames(annotations))
  if (length(missing)!=0)
    stop(paste0(paste(missing, collapse = ", "), ": covariates not found"))
  
  ord <- getMultiColOrder(annotations, sortBy)
  annotations[ord, , drop=F]
}

getMultiColOrder <- function(df, sortBy) {
  if (!is.data.frame(df))
    stop("df must be a data frame.")
  
  columns <- paste0("df$", sortBy)
  columns <- paste(columns, collapse = ", ")
  exp <- paste0("order(", columns, ")")
  eval(parse(text=exp))
}

extractPvalues <- function(x) {
  p <- x$pvalue
  if (is.null(p)) {
    return(NA)
  } else {
    return(p)
  }
}

matchArguments <- function(dots, defaults) {
  if (length(defaults)==0)
    return(dots)
  
  defaults[names(defaults) %in% names(dots)] <- NULL
  c(defaults, dots)
}

guessOmic <- function(covs) {
  unique(sub(omicsRegexp, "", covs, perl=TRUE, ignore.case=FALSE))
}

guessOmics <- function(covs) {
  sub(omicsRegexp, "", covs, perl=TRUE, ignore.case=FALSE)
}

guessOmicsColors <- function(omics) {
  uomics <- unique(omics)
  if (length(uomics) < 7){
    MOcols <- names(MOSpalette)[seq_along(uomics)]
  } else {
    MOcols <- RColorBrewer::brewer.pal(length(uomics), "Set3")  
  }
  names(MOcols) <- uomics
  MOcols
}

mapColor <- function(omic, MOcolors) {
  color <- MOSpalette[MOcolors[omic]]
  if (is.na(color))
    color <- MOcolors[omic]
  unname(color)
}

createColors <- function(omics, MOcolors) {
  sapply(unique(omics), function(o) mapColor(o, MOcolors))
}

matchAnnotations <- function(d1, d2){
  if (nrow(d1)!=nrow(d2))
    stop("Annotations have different row numbers")
  
  diff <- setdiff(row.names(d2), row.names(d1))
  if (length(diff) != 0)
    stop(paste0("We found samples", paste(diff, collapse = ", "),
                "that do not match MOM annotation"))
  
  d2 <- d2[row.names(d1), , drop=F]
  d2
}

getContinousPalette <- function(palette, n) {
  switch(palette,
         red = redShades(n),
         green = greenShades(n),
         blue = blueShades(n),
         yellow = yellowShades(n),
         violet = violetShades(n),
         teal = tealShades(n))
}

extractPositivePortion <- function(data, invert=FALSE) {
  .data <- data
  if (invert) {
    .data[data > 0] <- 0
    .data <- abs(.data)
  } else {
    .data[data < 0] <- 0
  }
  .data
}


#' Create the list of covariates that are going to be tested
#'
#' @importFrom methods new
#' @importFrom survival Surv
#' @return list with
#'   1 reduced representation of the omics
#'   2 sdev
#'   3 loadings or eigenvector
#'   4 usedGenes
#'   5 method
#'   6 namesCov
#'   7 omicName
#'
createMOMView <- function(omicsObj, genes) {
  listCovariates <- lapply(seq_along(omicsObj@ExperimentList@listData), 
                           function(i) {
    test <- get(omicsObj@modelInfo[i])
    specificArgs <- omicsObj@specificArgs[[i]]
    args <- list(data=omicsObj@ExperimentList@listData[[i]], features=genes)
    if (!is.null(specificArgs))
      args <- c(args, specificArgs)
    do.call(test, args)
  })
  
  listCovariates[!vapply(listCovariates, is.null, logical(1))] 
}


#' Create Cox Object
#' 
#' Create the coxObj from the covariates used in the test
#'
#' @param colData colData from multiOmic object
#' @param moView modulesView or pathView from multiOmicsModules or 
#' multiOmicsPathway object
#'
#' @return data.frame, samples in the rows, covariates in the columns
createCoxObj <- function(colData, moView){
  
  additionalCovariates <- lapply(moView, function(mo) {mo$x})
  additionalCovariates <- do.call(cbind, additionalCovariates)
  
  if (is.null(additionalCovariates))
    return(NULL)
  
  if (!identical(row.names(colData), row.names(additionalCovariates)))
    stop("Mismatch in covariates and daysStatus annotations rownames.")
  
  coxObj <- data.frame(colData, additionalCovariates)
  return(coxObj)
}


#' Create Data Module
#' 
#' Extract sub-matrix for the genes of a module or pathway from data matrix of 
#' a specific omic
#'
#' @param omic modulesView or pathView object 
#' @param multiOmicObj object of class 'Omics'
#' 
#' @return matrix, genes in the rows, samples in the columns
createDataModule <- function(omic, multiOmicObj){
  if (omic$omicName == "met"){
    genes <- unname(unlist(omic$cls))
  } else {
    genes <- omic$usedGenes
  }
  
  if (isFALSE(omic$omicName %in% names(multiOmicObj@ExperimentList))) {
    stop(paste0("omicName not found in ExperimentList.\n",
                "Names of experiments in ExperimentList should match",
                "the name arguments given in specificArgs"))
  }
  assay <- assay(multiOmicObj, omic$omicName)
  dataModule <- assay[genes, , drop = FALSE]
  return(dataModule)
}


#' Shows the MOSClip palette.
#'
#' This function shows the MOSClip palette. 
#' Each omic should be coupled to a color panel, this match will be preserved 
#' in plots.
#'
#' @examples
#' showMOSpalette()
#' 
#' @importFrom graphics axis title
#'
#' @export
showMOSpalette <- function(){
  plot(c(1:6,1:6,1:6,1:6), c(rep(1,6),rep(2,6), rep(3,6), rep(4,6)),
       col=apply(MOSpaletteSchema,2,c), pch=19, cex=13, xlim=c(0.5,6.5) , 
       ylim=c(0,4.8), xlab="", ylab="", axes=FALSE)
  axis(3, at=1:6, labels=rownames(MOSpaletteSchema), las=1, 
       cex.axis=0.8, lwd=0, pos=4.5, font=4)
  title("MOSClip palette")
}