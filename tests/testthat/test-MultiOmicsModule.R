
dummy_methylation_like_dataset <- function(seed=1234) {
  genesSize <- c(3,2)
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

  
  row.names(fake_data) <- paste0("gene_", seq_len(NROW(fake_data)))
  colnames(fake_data) <- paste0("p_", seq_len(NCOL(fake_data)))
  return(fake_data)
}

dummy_methylation_like_flat_dataset <- function(seed=1234) {
  flat_data <- matrix(c(rep(0.3,191), rep(0.4,9)),ncol=200, nrow=1,
                      dimnames = list(paste0("gene_", seq_len(1)),
                                      paste0("p_", seq_len(200))))
  
  row.names(flat_data) <- paste0("gene_", seq_len(NROW(flat_data)))
  colnames(flat_data) <- paste0("p_", seq_len(NCOL(flat_data)))
  return(flat_data)
}

dummy_expression_like_dataset <- function(seed=1234) {
  genesSize <- c(3,2)
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
  
  row.names(fake_exp) <- paste0("gene_", seq_len(NROW(fake_exp)))
  colnames(fake_exp) <- paste0("p_", seq_len(NCOL(fake_exp)))
  return(fake_exp)
}

dummy_cnv_like_dataset <- function(seed=1234) {
  set.seed(seed)
  genesSize <- c(3,2)
  patientsRatio <- list(c(150,50), c(75,125))
  
  lowValues <- sample(c(-1,0,0,0,1,1), patientsRatio[[1]][1]*genesSize[1],replace = T)
  highValues <- sample(c(-2,0,0,0,0,2,2), patientsRatio[[1]][2]*genesSize[1], replace = T)
  lowValues2nd <- sample(c(-2,0,0,0,0,0,2), patientsRatio[[2]][1]*genesSize[2], replace = T)
  highValues2nd <- sample(c(-1,0,0,0,0,0,1), patientsRatio[[2]][2]*genesSize[2], replace = T)
  
  fake_cnv <- rbind(cbind(matrix(lowValues, ncol=patientsRatio[[1]][1], nrow=genesSize[1]),
                          matrix(highValues, ncol=patientsRatio[[1]][2], nrow=genesSize[1])),
                    cbind(matrix(lowValues2nd, ncol=patientsRatio[[2]][1], nrow=genesSize[2]),
                          matrix(highValues2nd, ncol=patientsRatio[[2]][2], nrow=genesSize[2]))
  )
  
  row.names(fake_cnv) <- paste0("gene_", seq_len(NROW(fake_cnv)))
  colnames(fake_cnv) <- paste0("p_", seq_len(NCOL(fake_cnv)))
  return(fake_cnv)
}

create_met_cluster_dict <- function(dataset) {
  dict = as.list(row.names(dataset))
  names(dict) <- row.names(dataset)
  return(dict)
}

fake_data <- dummy_methylation_like_dataset()
dict <- create_met_cluster_dict(fake_data)

test_that("clusterWithDict", {
  name = "cls"
  cls <- summarizeInCluster(fake_data, features=row.names(fake_data), name = name)
  clsDict <- summarizeInCluster(fake_data, features=row.names(fake_data), dictionary = dict, name = name)
  expect_equal(cls$namesCov,paste0(name, "3k"))
  expect_identical(cls$x, clsDict$x)
})


test_that("emptyInput", {
  expect_null(summarizeInCluster(data=NULL, features=row.names(fake_data)))
})


flat_data <- dummy_methylation_like_flat_dataset()

test_that("singleProfileInput2class", {
  name = "cls"
  cls <- summarizeInCluster(flat_data, features=row.names(flat_data), name = name)
  expect_equal(cls$namesCov,paste0(name, "2k"))
})

test_that("notIntersectionOfGenes", {
  expect_null(summarizeInCluster(flat_data, features="gene_2"))
})


fake_exp <- dummy_expression_like_dataset()
test_that("computePCs", {
  pcs <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="regular", cliques=NULL, maxPCs=10, loadThr=0.6)
  pcsS <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="sparse", cliques=NULL, maxPCs=10, loadThr=0.6)
  pcsT <- summarizeWithPca(fake_exp, features=row.names(fake_exp), name="pca", shrink=FALSE, method="topological",
                           cliques=list(row.names(fake_exp),
                                        row.names(fake_exp)[5]), maxPCs=10, loadThr=0.6)
  expect_equal(dim(pcs$x),c(200,2))
  expect_equal(dim(pcsS$x),c(200,2))
  expect_equal(dim(pcsT$x),c(200,2))
})


fake_cnv <- dummy_cnv_like_dataset()

test_that("cnv", {
  summary_cnvbin <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv))
  summary_cnv <- summarizeToNumberOfDirectionalEvents(fake_cnv, row.names(fake_cnv))
  
  expect_identical(summary_cnvbin$x$dirBinPOS, unname(apply(fake_cnv > 1, 2, any, na.rm=T)))
  expect_identical(summary_cnvbin$x$dirBinNEG, unname(apply(fake_cnv < -1, 2, any, na.rm=T)))
  expect_identical(summary_cnv$x$dCountPOS, unname(apply(fake_cnv > 1, 2, sum, na.rm=T)))
  expect_identical(summary_cnv$x$dCountNEG, unname(apply(fake_cnv < -1, 2, sum, na.rm=T)))
})


fake_mut <- fake_cnv
fake_mut[abs(fake_cnv) < 2] <- 0
fake_mut[abs(fake_cnv) >= 2] <- 1

test_that("mut", {
  summary_mutbin <- summarizeToBinaryEvents(fake_mut, row.names(fake_mut))
  summary_mut <- summarizeToNumberOfEvents(fake_mut, row.names(fake_mut))
  
  expect_identical(summary_mut$x$event, unname(apply(fake_mut == 1, 2, sum, na.rm=T)))
  expect_identical(summary_mutbin$x$bin, unname(apply(fake_mut == 1, 2, any, na.rm=T)))
})