# Creates fake matrices
dummy_mutation_like_dataset <- function(genes=NULL, seed=1234) {
  set.seed(seed)
  g_len <- if (!is.null(genes)) length(genes) else 5
  len = 200
  mut = matrix(sample(c(0,0,0,1), len*g_len, replace=T), ncol=len, nrow=g_len)
  colnames(mut) <- paste0("P_", seq_len(NCOL(mut)))
  if (is.null(genes)) {
    row.names(mut) <- paste0("gene_", seq_len(NROW(mut)))
  } else {row.names(mut) <- genes}
  return(mut)}

dummy_cnv_like_dataset <- function(genes=NULL, seed=1234) {
  set.seed(seed)
  if (!is.null(genes)) {
    num1 <- floor(length(genes)/2)
    genesSize <- c(num1, length(genes)-num1)
  } else {genesSize <- c(3,2)}
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- sample(c(-1,0,0,0,1,1), patientsRatio[[1]][1]*genesSize[1],replace = T)
  highValues <- sample(c(-2,0,0,0,0,2,2), patientsRatio[[1]][2]*genesSize[1], replace = T)
  lowValues2nd <- sample(c(-2,0,0,0,0,0,2), patientsRatio[[2]][1]*genesSize[2], replace = T)
  highValues2nd <- sample(c(-1,0,0,0,0,0,1), patientsRatio[[2]][2]*genesSize[2], replace = T)
  
  fake_cnv <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1], nrow=genesSize[1]),
                          matrix(highValues, ncol=patientsRatio[[1]][2], nrow=genesSize[1])),
                    cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1], nrow=genesSize[2]),
                          matrix(highValues2nd, ncol=patientsRatio[[2]][2], nrow=genesSize[2])))
  
  if (is.null(genes)) {
    row.names(fake_cnv) <- paste0("gene_", seq_len(NROW(fake_cnv)))}
  else {row.names(fake_cnv) <- genes}
  colnames(fake_cnv) <- paste0("P_", seq_len(NCOL(fake_cnv)))
  return(fake_cnv)}

dummy_methylation_like_dataset <- function(genes=NULL, seed=1234) {
  if (!is.null(genes)) {
    num1 <- floor(length(genes)/2)
    genesSize <- c(num1, length(genes)-num1)
  } else {genesSize <- c(3,2)}
  patientsRatio <- list(c(150,50), c(75,125))
  
  set.seed(seed)
  lowValues <- runif(patientsRatio[[1]][1]*genesSize[1], 0, 0.2)
  highValues <- runif(patientsRatio[[1]][2]*genesSize[1], 0.7, 0.99)
  lowValues2nd <- runif(patientsRatio[[2]][1]*genesSize[2], 0.6, 1)
  highValues2nd <- runif(patientsRatio[[2]][2]*genesSize[2], 0.2, 0.4)
  
  fake_data <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1], nrow=genesSize[1]),
                           matrix(highValues, ncol=patientsRatio[[1]][2], nrow=genesSize[1])),
                     cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1], nrow=genesSize[2]),
                           matrix(highValues2nd, ncol=patientsRatio[[2]][2], nrow=genesSize[2]))
  )
  
  
  if (is.null(genes)) {
    row.names(fake_data) <- paste0("gene_", seq_len(NROW(fake_data)))}
  else {row.names(fake_data) <- genes}
  colnames(fake_data) <- paste0("P_", seq_len(NCOL(fake_data)))
  return(fake_data)
}

dummy_methylation_like_flat_dataset <- function(genes=NULL, seed=1234) {
  flat_data <- matrix(c(rep(0.3,191), rep(0.4,9)), ncol=200, nrow=1,
                      dimnames = list(paste0("gene_", seq_len(1)),
                                      paste0("P_", seq_len(200))))
  
  row.names(flat_data) <- paste0("gene_", seq_len(NROW(flat_data)))
  colnames(flat_data) <- paste0("P_", seq_len(NCOL(flat_data)))
  return(flat_data)
}

dummy_expression_like_dataset <- function(genes=NULL, seed=1234) {
  set.seed(seed)
  if (!is.null(genes)) {
    num1 <- floor(length(genes)/2)
    genesSize <- c(num1, length(genes)-num1)
  } else {genesSize <- c(3,2)}
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- rnorm(patientsRatio[[1]][1]*genesSize[1], 5, 1)
  highValues <- rnorm(patientsRatio[[1]][2]*genesSize[1], 8, 1)
  lowValues2nd <- rnorm(patientsRatio[[2]][1]*genesSize[2], 9, 1)
  highValues2nd <- rnorm(patientsRatio[[2]][2]*genesSize[2], 2, 1)
  
  fake_exp <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1], nrow=genesSize[1]), 
                          matrix(highValues, ncol=patientsRatio[[1]][2], nrow=genesSize[1])),
                    cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1], nrow=genesSize[2]), 
                          matrix(highValues2nd, ncol=patientsRatio[[2]][2], nrow=genesSize[2]))
  )
  
  if (is.null(genes)) {
    row.names(fake_exp) <- paste0("gene_", seq_len(NROW(fake_exp)))}
  else {row.names(fake_exp) <- genes}
  colnames(fake_exp) <- paste0("P_", seq_len(NCOL(fake_exp)))
  return(fake_exp)
}

create_met_cluster_dict <- function(dataset) {
  dict = as.list(row.names(dataset))
  names(dict) <- row.names(dataset)
  return(dict)
}

dummy_colData <- function(len=10, type="survival"){
  num1 <- sample(1:(len-1), 1)
  num2 <- len - num1
  if (type=="two-classes") {
    dummy_colData <- data.frame(classes = c(rep("A", num1), rep("B", num2)))}
  else if (type=="survival") {
    dummy_colData <- data.frame(status = sample(c(0,1), size=len, replace=TRUE),
                                days = sample(5:5000, size=len, replace=TRUE))
  }
  rownames(dummy_colData) <- paste0("P_", seq_len(len))
  return(dummy_colData)}


fake_mo <- function(genes=NULL, type="survival",
                    omics=c("exp", "met", "mut", "cnv"), modelInfo=NULL,
                    specificArgs=NULL){
  exp <- dummy_expression_like_dataset(genes=genes)
  mut <- dummy_mutation_like_dataset(genes=genes)
  met <- dummy_methylation_like_dataset(genes=genes)
  cnv <- dummy_cnv_like_dataset(genes=genes)
  exp_list <- ExperimentList(exp=exp, mut=mut, met=met, cnv=cnv)
  cdata <- dummy_colData(len=dim(exp)[2], type)
  modelInfo = c(exp="summarizeWithPca", met="summarizeInCluster",
                mut="summarizeToNumberOfEvents", cnv="summarizeToNumberOfDirectionalEvents")
  specificArgs = list(exp=list(name="exp"), met=list(name="met"),
                      mut=list(name="mut"), cnv=list(name="cnv"))
  mo <- makeOmics(experiments=exp_list[omics],
                  colData=cdata,
                  modelInfo = modelInfo[omics],
                  specificArgs = specificArgs[omics])
  return(mo)
}

