#' Estimate Single Covariance Matrix
#'
#' For internal use only. Estimate Covariance from one matrix
#'
#' @param expr a numeric matrix
#' @inheritParams estimateCov
#'
#' @return a covariance matrix
#' @importFrom stats cov
#' @importFrom corpcor cov.shrink
#' @export
#'
estimateExprCov <- function(expr, shrink) {
  if (shrink)
    unclass(corpcor::cov.shrink(expr, verbose=FALSE))
  else
    cov(expr)
}

#' Make positive and definite covariance matrix
#'
#' @param m1 matrix 1
#' @param m2 matrix 2
#' @param m3 matrix 3
#' @param threshold threshold of difference
#'
#' @return list with
#'   \item{m1}{the matrix m1 positive and definite}
#'   \item{m2}{the matrix m2 positive and definite}
#'   \item{m3}{the matrix m3 positive and definite}
#'   \item{correction}{the magneturde of the correction}
#'   \item{value}{the value}
#'
#' @export
makePositiveDefinite<-function(m1, m2=NULL, m3=NULL, threshold=0.1){
  
  if(is.null(m2) & is.null(m3)) {
    eig<-min(round(eigen(m1)$values,2))
    
  } else {
    if (any(is.null(c(m2,m3)))) stop("Both m2 and m3 are needed")
    
    if (!(ncol(m1)==nrow(m1)) | !(ncol(m2)==nrow(m2)) | (!(ncol(m3)==nrow(m3))))
      stop("the matrices are not square")
    
    if (!all(rownames(m1)==colnames(m1)) | !all(rownames(m2)==colnames(m2)) | !all(rownames(m3)==colnames(m3)))
      stop("column and row elements must be equal")
    
    if( !all(dim(m1)==dim(m2)) | !all(dim(m1)==dim(m3)))
      stop("Matrices have different sizes")
    
    if (!all(colnames(m1)==colnames(m2)) | !all(colnames(m1)==colnames(m3)))
      stop("Matrices hav different names")
    
    eig<-min(round(c(eigen(m1)$values,eigen(m2)$values,eigen(m3)$values),2))
  }
  
  if(!(eig>0)) {
    n<-ncol(m1)
    
    m1n <- m1+(threshold-eig)*diag(n)
    
    m2n <- if(!is.null(m2)) m2+(threshold-eig)*diag(n) else NULL
    m3n <- if(!is.null(m3)) m3+(threshold-eig)*diag(n) else NULL
    
    correction<- TRUE
    value<-threshold-eig
    
  } else {
    m1n<-m1; m2n<-m2; m3n<-m3
    correction<-FALSE
    value<-NULL
  }
  
  return(list(m1=m1n, m2=m2n, m3=m3n, correction=correction, value=value))
}