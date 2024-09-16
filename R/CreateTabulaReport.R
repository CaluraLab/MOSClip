#' Summarize pathways' info from a list of MultiOmicsPathway objects (MOP)
#'
#' Given the list of MOPs, it creates the table.
#'
#' @param multiPathwayList a list of `MultiOmicsPathway` objects resulting from 
#' a multi-omics pathway test.
#' @param priority_to a vector with the covariates (omic name) that should go 
#' first.
#'
#' @return a data.frame, pathways in rows, overall pvalue of the coxph, 
#' followed by covariates pvalues, in columns.
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#' 
#' genesToUse <- row.names(multiOmics[[1]])
#' 
#' MOP_list <- lapply(reactSmall, function(g) {
#'    print(g@title)  #  to see which pathways are being calculated
#'    set.seed(1234)
#'    fcl = multiOmicsSurvivalPathwayTest(multiOmics, g,
#'                                        survFormula="Surv(days, status) ~",
#'                                        autoCompleteFormula = TRUE,
#'                                        useTheseGenes = genesToUse)
#'    fcl
#' })
#' 
#' pathwaysSummary <- multiPathwayReport(MOP_list)
#' 
#'
#' @export
multiPathwayReport <- function(multiPathwayList, priority_to=NULL){
  if (!is(multiPathwayList, "list") ||
      any(vapply(multiPathwayList, class, character(1)) != 
          "MultiOmicsPathway"))
    stop("A list of pathway results are expected.")

  pvalues <- vapply(multiPathwayList, function(p) {as.numeric(p@pvalue)}, 
                    numeric(1))

  zs <- sort(unique(unlist(lapply(multiPathwayList, function(p) {
    names(p@zlist)
  }))))

  zMat <- do.call(rbind, lapply(multiPathwayList, function(p) {
    fixedCols <- rep(NA, length(zs))
    names(fixedCols) <- zs
    fixedCols[names(p@zlist)] <- p@zlist
    fixedCols
  }))

  ord <- order(pvalues)
  df <- cbind(row.names=names(pvalues)[ord], pvalue=pvalues[ord],
              data.frame(zMat[ord, , drop=FALSE]))
  order_by_covariates(df, 1, priority_to)
}


#' Provides a Table of the Modules Test Results
#'
#' Summarizes the results of a multi omics module test given a list of
#' MultiOmicsModules objects
#'
#' @param multiPathwayModuleList a list of `MultiOmicsModules` objects 
#' resulting from a multi-omics module test.
#' @param priority_to a vector with the covariates (the omics names) 
#' that should appear first in the dataframe columns
#'
#' @return a data.frame class object. Rows correspond to the modules, and the
#' columns to the overall and covariates pvalues of the test.
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#' 
#' genesToUse <- row.names(multiOmics[[1]])
#' 
#'  MOM_list <- lapply(reactSmall[1:2], function(g) {
#'    print(g@title)  #  to see which pathways are being calculated
#'    set.seed(1234)
#'    fcl = multiOmicsSurvivalModuleTest(multiOmics, g,
#'                                        survFormula="Surv(days, status) ~",
#'                                        autoCompleteFormula = TRUE,
#'                                        useTheseGenes = genesToUse)
#'    fcl
#' })
#' 
#'  moduleSummary <- multiPathwayModuleReport(MOM_list)
#' 
#' @export

multiPathwayModuleReport <- function(multiPathwayModuleList, priority_to=NULL) {
  if (!is(multiPathwayModuleList, "list") ||
      any(vapply(multiPathwayModuleList, class, character(1)) != 
          "MultiOmicsModules"))
    stop("A list of pathway modules results are expected.")

  n_temp <- names(multiPathwayModuleList)

  multiMatrixRes <- lapply(n_temp, function(name,list){
    summary <- formatModuleReport(list[[name]]);
    data.frame(pathway=name, module=row.names(summary), summary,
               row.names=NULL, stringsAsFactors = FALSE) }, 
    multiPathwayModuleList)
  resDF <- mergeAll(multiMatrixRes)
  resDF <- resDF[order(resDF$pvalue), ]
  rownames(resDF) <- apply(resDF,1, function(r) paste(r["pathway"],
                                                      r["module"],sep="."))

  resDF <- order_by_covariates(resDF, 3, priority_to)

  return(resDF)
}



formatModuleReport <- function(smObj){
  alphas <- smObj@alphas
  z  <- smObj@zlists
  idxs <- order(alphas)

  zcols <- sort(unique(unlist(lapply(z, function(x){
    names(x)
  }))))

  colDescription <- do.call(rbind, lapply(idxs, function(i){
    additionalCols <- rep(NA, length(zcols))
    names(additionalCols) <- zcols
    additionalCols[names(z[[i]])] <- z[[i]]
    additionalCols
  }))

  cbind(row.names=idxs, pvalue=alphas[idxs], data.frame(colDescription))
}

mergeAll <- function(list) {
  allColumnsNames <- sort(unique(unlist(lapply(list, function(o) {
    colnames(o)
  }))))

  matrix <- do.call(rbind,
                    lapply(list, function(o) {
                      fixedCols <- matrix(NA, NROW(o), length(allColumnsNames))
                      colnames(fixedCols) <- allColumnsNames
                      for (col in colnames(o)) {
                        fixedCols[, col] <- o[,col]
                      }
                      fixedCols
                    }))
  removeCols <- match(c("pathway", "module","pvalue"), colnames(matrix))
  numericMat <- matrix[,-removeCols, drop=FALSE]
  data.frame(pathway=matrix[, "pathway"], module=matrix[, "module"],
             pvalue=as.numeric(matrix[, "pvalue"]),
             apply(numericMat, 2, as.numeric), stringsAsFactors = FALSE)
}

order_by_covariates <- function(dataF, skip_first_cols, priority_to=NULL) {


  if (is.null(priority_to))
    return(dataF)

  covariates <- dataF[, seq_len(ncol(dataF)-skip_first_cols)+skip_first_cols,
                      drop=FALSE]
  fixed <- dataF[, seq_len(skip_first_cols), drop=FALSE]
  omics_cov <- guessOmics(colnames(covariates))
  to_sort <- unique(omics_cov) %in% priority_to
  priority_to <- c(priority_to, unique(omics_cov)[!to_sort])
  cov_factor <- factor(omics_cov, levels=priority_to)
  cov_factor <- droplevels(cov_factor)
  covariates <- covariates[, order(cov_factor), drop=FALSE]
  cbind(fixed, covariates)
}
