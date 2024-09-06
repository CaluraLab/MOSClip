#' Performs a Exact test - analysis of omics intersection
#'
#' This function performs a exact test implementing a theorical framework using
#' the SuperExactTest package. It calculates the statistical distributions of
#' multi omics set intersections. It can be used with both a MultiOmicsModules 
#' or MultiOmicsPathway class objects.
#'
#' @param multiPathwayReportData data.frame, the output of the 
#' \code{\link{multiPathwayReport}} or \code{\link{multiPathwayModuleReport}}
#' functions.
#' @param pvalueThr numeric value. Overall pvalue cut-off to be used
#' @param zscoreThr numeric value. Covariates coefficient cut-off to be used.
#' @param resampligThr numeric value. Filters the modules according to the
#' number of success in the resampling procedure, takes only the modules above
#' this threshold.
#' @param plot character indicating the layout for plotting. It is one of 
#' \code{circular}, \code{landscape} or \code{noplot}. By default,
#' plot="circular", if plot="noplot" no plot will be provided.
#' @param sort.by character indicating how to sort the
#' intersections in the plot. It is one of "set" (by omics), "size"
#' (by intersection size), "degree" (by number of intersected omics),
#' and "p-value".
#' @param excludeColumns a vector of characters listing the columns of
#' \code{multiPathwayReportData} object to be excluded by the analysis.
#' In the case \code{multiPathwayReportData} derives from
#' \code{\link{multiPathwayModuleReport}}
#' you should set \code{excludeColumns = c("pathway","module")}.
#' @param color.on color that represent the active omics in the sector
#' @param color.off color that represent the omics mnot considered in the
#' sector
#'
#' @details This function calculates intersection sizes between multiple set of
#' pathways or modules and performs statistical test of the intersections using
#' the total amout of analyzed pathways or modules as background. The super
#' exact test of this function was described in Wang et al 2015.
#'
#' @return a data.frame containing all the numeric information of the plot
#' included the pathways shared by different omics.
#'
#' @references
#' Minghui Wang, Yongzhong Zhao, and Bin Zhang (2015).
#' Efficient Test and Visualization of Multi-Set Intersections.
#' Scientific Reports 5: 16923.
#'
#' @importFrom SuperExactTest supertest
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
runSupertest <- function(multiPathwayReportData, pvalueThr=0.05,
                         zscoreThr=0.05,
                         resampligThr = NULL,
                         plot=c('circular','landscape','noplot'),
                         sort.by=c('set','size','degree','p-value'),
                         excludeColumns=NULL,
                         color.on = "#f6bb42", color.off = "#D3D3D3"){

  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)

  checkPvalueThresholdFormat(pvalueThr, "pvalueThr")
  checkPvalueThresholdFormat(zscoreThr, "zscoreThr")

  plot <- plot[1]
  if(!(plot %in% c('circular','landscape','noplot')))
    stop("Plot argument should be either circular, landscape or noplot." )

  sort.by <- sort.by[1]
  if(plot != "noplot" & (!(sort.by %in% c('set','size','degree','p-value'))))
    stop("sort.by argument should be either set, size, degree or p-value.")

  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[
    multiPathwayReportData[,"pvalue"] <= pvalueThr,]

  if (!is.null(resampligThr)) {
    if (is.null(multiPathwayReportDataSig$resamplingCount))
      stop("resamplingCount column not found. Try addResamplingCounts")
    multiPathwayReportDataSig <- multiPathwayReportDataSig[
      multiPathwayReportDataSig[,"resamplingCount"] >= resampligThr,]
  }

  MOlistPval <- pvalueSummary(multiPathwayReportDataSig,
                              excludeColumns = excludeColumns, as.list=TRUE)

  MOlistPathSig <- lapply(MOlistPval, function(pp) {
    names(which(pp <= zscoreThr))})

  msetSupertest <- SuperExactTest::supertest(MOlistPathSig, n=universeSize)

  if(plot != "noplot"){
    plot(msetSupertest,
         color.on = color.on, color.off = color.off,
         heatmapColor = rev(pvalueShades),
         sort.by = sort.by, Layout = plot)
  }

  invisible(summary(msetSupertest)$Table)
}

#' Find Pathway Fathers
#' 
#' Given the hierarchy of the pathways, this formula finds the fathers of the
#' respective pathway (e.g. pathway: "PI3K Cascade"; father:
#' "Signaling Pathways")
#'
#' @param pathways vector of pathway names
#' @param graphiteDB graphite DB object
#' @param hierarchy a graph object with the pathway hierarchy
#'
#' @return a vector of the fathers' names
#'
#' @importFrom igraph V
#
#' @export

annotePathwayToFather <- function(pathways, graphiteDB, hierarchy) {
  ord <- length(igraph::V(hierarchy))*2
  pathway2id <- mapPathwaysIDfromGraphite(graphiteDB) # codici
  pathwayDict <- pathway2id$pname
  names(pathwayDict) <- pathway2id$id

  ids <- mapPathwaysIDfromGraphite(graphiteDB, pathways)$id
  path2fathers <- lapply(ids, getPathFathers, hierarchy,
                         ord=ord)
  names(path2fathers) <- ids
  ids2father <- id2name(path2fathers, pathwayDict)
  unlist(ids2father)
}

#' Compute Omics Intersections
#'
#' Finds the modules that have any intersection among the available omics
#'
#' @inheritParams runSupertest
#'
#' @return a list of pathway modules present for every intersection of omics
#' present
#'
#' @importFrom reshape melt
#' @export

computeOmicsIntersections <- function(multiPathwayReportData, pvalueThr=0.05,
                                      zscoreThr=0.05, resampligThr = NULL,
                                      excludeColumns=NULL){

  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)

  checkPvalueThresholdFormat(pvalueThr, "pvalueThr")
  checkPvalueThresholdFormat(zscoreThr, "zscoreThr")

  universeSize <- NROW(multiPathwayReportData)
  multiPathwayReportDataSig <- multiPathwayReportData[
    multiPathwayReportData[,"pvalue"] <= pvalueThr,]
  
  if (nrow(multiPathwayReportDataSig) == 0) {
    stop("There is no significant modules. Try changing the pvalue threshold")
  }

  if (!is.null(resampligThr)) {
    if (is.null(multiPathwayReportDataSig$resamplingCount))
      stop("resamplingCount column not found. Try addResamplingCounts")
    multiPathwayReportDataSig <- multiPathwayReportDataSig[
      multiPathwayReportDataSig[,"resamplingCount"] >= resampligThr,]
  }

  MOlistPval <- pvalueSummary(multiPathwayReportDataSig,
                              excludeColumns = excludeColumns, as.list=TRUE)

  MOlistPathSig <- lapply(MOlistPval, function(pp) {
    names(which(pp <= zscoreThr))})
  
  if (0 %in% lengths(MOlistPathSig)){
    stop("One or more omics do not have a significant z score for any of ",
         "their modules. Try increasing the z score threshold.")
  }

  df <- reshape::melt(MOlistPathSig)
  p2o <- tapply(seq_len(NROW(df)), df[,1], function(idx) {
    paste(df[idx,2], collapse = ";")
  })
  tapply(seq_along(p2o), p2o, function(idx) {
    names(p2o[idx])
  })
}

#' Remove Module Number From Pathway Name
#' 
#' Function to remove the suffix corresponding to the module number of the
#' pathway name. Necessary step for \code{\link{annotePathwayToFather}} and 
#' \code{\link{plotFrequencies}}
#'
#' @inheritParams annotePathwayToFather
#'
#' @return list of pathway names without the module number
#'
#' @export

stripModulesFromPathways <- function(pathways) {
  sub("\\.[0-9]+", "",pathways, perl = TRUE)
}

#' Remove pathways with no id in graphite
#'
#' This function id risky to use. It performs pairing without
#' checking proper pairs.
#'
#' @param omicClass the name of an omic class (or omic combination)
#' @param omicsClasses2pathways vector of pathway names
#' @param omicsClasses2fathers vector of pathway names
#'
#' @return data frame
#'
#' @export
#'
match_pathway_to_fathers <- function(omicClass, omicsClasses2pathways,
                                     omicsClasses2fathers) {
  if (!(omicClass %in% names(omicsClasses2pathways)))
    stop("Omic not found in class2pathways")
  if (!(omicClass %in% names(omicsClasses2fathers)))
    stop("Omic not found in class2pathways")

  if (is.null(omicsClasses2fathers[[omicClass]])){
    valid_set <- character()
  } else {
    valid_set <- intersect(omicsClasses2pathways[[omicClass]],
                           omicsClasses2fathers[[omicClass]])
  }
  leave_out <- setdiff(omicsClasses2pathways[[omicClass]], valid_set)
  if (length(leave_out)>0)
    warning(paste0(
      "The following pathways were excluded because no id was found; ",
      paste(leave_out, collapse = ", ")))

  select_path <- omicsClasses2pathways[[omicClass]] %in% valid_set
  select_father <- omicsClasses2fathers[[omicClass]] %in% valid_set
  data.frame(path=omicsClasses2pathways[[omicClass]][select_path],
             father=omicsClasses2fathers[[omicClass]][select_father],
             stringsAsFactors = FALSE)
}


#' Compute pvalue Summary
#'
#' @inheritParams runSupertest
#' @param as.list return a list rather than a data.frame
#'
#' @return a list
#'
#' @importFrom stats na.omit
#' @export
#'
pvalueSummary <- function(multiPathwayReportData, excludeColumns = NULL,
                          as.list = FALSE){
  checkReportFormat(multiPathwayReportData)
  checkColumnsExclusion(multiPathwayReportData, excludeColumns)

  columnsNotExcluded <- colnames(multiPathwayReportData)[
    !(colnames(multiPathwayReportData) %in% excludeColumns)]
  multiPathwayReportData <- multiPathwayReportData[,columnsNotExcluded,
                                                   drop=FALSE]


  colClasses <- vapply(multiPathwayReportData, class, character(1))
  if(any(unique(colClasses) != "numeric")){
    notNumericColumns <- colnames(multiPathwayReportData)[
      colClasses != "numeric"]
    stop(paste0("Data malformed. ",
                "The following columns are not numeric.
                Consider using excludeColumns argument: ",
                paste(notNumericColumns, collapse = ", ")))
  }

  covarColumns <- colnames(multiPathwayReportData) != "pvalue"
  multiPathwayReportDataSig <- multiPathwayReportData[,covarColumns,
                                                      drop=FALSE]

  malformedColumns <- apply(multiPathwayReportDataSig, 2,
                            function(col) 
                              any(na.omit(col) > 1 | na.omit(col) < 0))

  if (any(malformedColumns)) {
    stop(paste0("Data malformed. The following columns are not pvalues
                (values greater than 1 or lower than 0).
                Consider using excludeColumns argument: ",
                paste(colnames(multiPathwayReportDataSig)[malformedColumns],
                      collapse = ", ")
    ))
  }

  covars <- colnames(multiPathwayReportDataSig)
  covars2omics <- guessOmics(covars)

  MOlistPval <- tapply(colnames(multiPathwayReportDataSig),
                       covars2omics,
                       summarizeOmicsResByMinPvalue,
                       mat=multiPathwayReportDataSig,
                       simplify = FALSE)
  if (as.list)
    return(MOlistPval)
  do.call(cbind, MOlistPval)
}


checkReportFormat <-function(multiPathwayReportData) {
  if(!(any("pvalue" %in% colnames(multiPathwayReportData))))
    stop("Data malformed. There is not a overall pvalue column.")

  if(is.null(grep(omicsRegexp, colnames(multiPathwayReportData))))
    stop("Data malformed. Covariates names not in colnames")
}

checkColumnsExclusion <- function(multiPathwayReportData, excludeColumns) {
  if (is.null(excludeColumns))
    return()

  if (any(!(excludeColumns %in% colnames(multiPathwayReportData))))
    stop("Data malformed. excludeColumns not present in data")

  if("pvalue" %in% excludeColumns)
    stop("You cannot exclude the overall pvalue column.")
}

checkPvalueThresholdFormat <- function(thr, name="thr") {
  if(!is.numeric(thr))
    stop("pValue threshold should be numeric")
  else if((thr > 1) | (thr < 0))
    stop("pValue threshold should be a number between 0 and 1")
}
