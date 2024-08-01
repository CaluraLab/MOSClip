# Creates data for two class data

dummy_methylation <- function(seed=1234) {

  set.seed(seed)
  
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  genesSize <- c(8,5)
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- runif(patientsRatio[[1]][1]*genesSize[1], 0, 0.2)
  highValues <- runif(patientsRatio[[1]][2]*genesSize[1], 0.7, 0.99)
  lowValues2nd <- runif(patientsRatio[[2]][1]*genesSize[2], 0.6, 1)
  highValues2nd <- runif(patientsRatio[[2]][2]*genesSize[2], 0.2, 0.4)
  
  fake_data <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1],
                                  nrow=genesSize[1]),
                           matrix(highValues, ncol=patientsRatio[[1]][2],
                                  nrow=genesSize[1])),
                     cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1],
                                  nrow=genesSize[2]),
                           matrix(highValues2nd, ncol=patientsRatio[[2]][2],
                                  nrow=genesSize[2]))
  )
  
  row.names(fake_data) <- genes
  colnames(fake_data) <- paste0("p_", seq_len(NCOL(fake_data)))
  return(fake_data)
}

dummy_expression <- function(seed=1234) {
  
  set.seed(seed)
  
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  genesSize <- c(8,5)
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- rnorm(patientsRatio[[1]][1]*genesSize[1], 5, 1)
  highValues <- rnorm(patientsRatio[[1]][2]*genesSize[1], 8, 1)
  lowValues2nd <- rnorm(patientsRatio[[2]][1]*genesSize[2], 9, 1)
  highValues2nd <- rnorm(patientsRatio[[2]][2]*genesSize[2], 2, 1)
  
  fake_exp <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1],
                                 nrow=genesSize[1]), 
                          matrix(highValues, ncol=patientsRatio[[1]][2],
                                 nrow=genesSize[1])),
                    cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1],
                                 nrow=genesSize[2]), 
                          matrix(highValues2nd, ncol=patientsRatio[[2]][2],
                                 nrow=genesSize[2]))
  )
  
  row.names(fake_exp) <- genes
  colnames(fake_exp) <- paste0("p_", seq_len(NCOL(fake_exp)))
  return(fake_exp)
}

dummy_cnv <- function(seed=1234) {
  
  set.seed(seed)
  
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  genesSize <- c(8,5)
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- sample(c(-1,0,0,0,1,1), patientsRatio[[1]][1]*genesSize[1],
                      replace = T)
  highValues <- sample(c(-2,0,0,0,0,2,2), patientsRatio[[1]][2]*genesSize[1],
                       replace = T)
  lowValues2nd <- sample(c(-2,0,0,0,0,0,2), patientsRatio[[2]][1]*genesSize[2],
                         replace = T)
  highValues2nd <- sample(c(-1,0,0,0,0,0,1), patientsRatio[[2]][2]*genesSize[2],
                          replace = T)
  
  fake_cnv <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1],
                                 nrow=genesSize[1]),
                          matrix(highValues, ncol=patientsRatio[[1]][2],
                                 nrow=genesSize[1])),
                    cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1],
                                 nrow=genesSize[2]),
                          matrix(highValues2nd, ncol=patientsRatio[[2]][2],
                                 nrow=genesSize[2]))
  )
  
  row.names(fake_cnv) <- genes
  colnames(fake_cnv) <- paste0("p_", seq_len(NCOL(fake_cnv)))
  return(fake_cnv)
}

fake_mut <- dummy_cnv()
fake_mut[abs(fake_mut) < 2] <- 0
fake_mut[abs(fake_mut) >= 2] <- 1

dummy_colData <- data.frame(patient_id = paste0("p_", seq_len(200)),
                            treatment = c(rep("A", 100), rep("B", 100)))
rownames(dummy_colData) <- paste0("p_", seq_len(200))

createOmics <- function(){
  dummy_twoclass_omics <-  Omics(
    experiments = ExperimentList(expr =  dummy_expression(),
                                 met = dummy_methylation(),
                                 cnv = dummy_cnv(),
                                 mut = fake_mut),
    colData = dummy_colData,
    modelInfo = c("summarizeWithPca", "summarizeInCluster",
                  "summarizeToNumberOfEvents",
                  "summarizeToNumberOfDirectionalEvents"),
    specificArgs = list(pcaArgs = list(name = "exp", shrink = "FALSE",
                                       method = "sparse", maxPCs = 3),
                        clusterArgs = list(name = "met",
                                           max_cluster_number = 3),
                        countEvent = list(name = "mut", min_prop = 0.05),
                        cnvAgv = list(name = "cnv", min_prop = 0.05)
    )
  )
}

createClassAnnot <- function(){
  dummy_classAnnot <- data.frame(classes = c(rep("class1", 100),
                                             rep("class2", 100)),
                                 row.names = paste0("p_",
                                                    seq_len(200)))
}


