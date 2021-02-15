genes = 4 
set.seed(1234)
mut = matrix(sample(c(0,0,0,1), 150*genes, replace=T), ncol=150, nrow=genes)
colnames(mut) <- paste0("P_", seq_len(NCOL(mut)))
row.names(mut) <- paste0("gene_", seq_len(NROW(mut)))

x = data.frame(mut = apply(mut > 0, 2, any), stringsAsFactors = F)
discrete <- x
discrete$mut <- as.numeric(x$mut)

mutationCount <- rowSums(mut)
ordMutCount <- mutationCount[order(mutationCount, decreasing = T)]
omic <- list(x=x,
             dataModule=mut,
             namesCov="mut",
             method="binary",
             omicName="mut",
             eventThr=1
)

test_that("extractyInfo_binary", {
  binInfo <- extractSummaryFromBinary(omic, n=3)
  expect_identical(row.names(binInfo$sigModule), head(names(ordMutCount),3))
  expect_identical(binInfo$discrete, discrete)
})


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
fake_cnv <- dummy_cnv_like_dataset()
omic <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv))

cnvCountPos <- rowSums(omic$dataModule >=2)
cnvCountNeg <- rowSums(omic$dataModule <=-2)
ordCnvCountPos <- cnvCountPos[order(cnvCountPos, decreasing = T)]
ordCnvCountNeg <- cnvCountNeg[order(cnvCountNeg, decreasing = T)]
ordCnvGene <- unique(c(head(names(ordCnvCountPos), 3), head(names(ordCnvCountNeg), 3)))

discrete2class <- data.frame(lapply(omic$x, as.numeric))
row.names(discrete2class) <- row.names(omic$x)

test_that("extractyInfo_dirBinary", {
  duobleBinInfo <- extractSummaryFromBinary(omic, n=3)
  expect_identical(row.names(duobleBinInfo$sigModule), ordCnvGene)
  expect_identical(duobleBinInfo$discrete, discrete2class)
})