test_that("extractyInfo_binary", {
  mut <- dummy_mutation_like_dataset()
  expect_null(summarizeToNumberOfDirectionalEvents(mut, row.names(mut)))
  x = data.frame(mut = apply(mut > 0, 2, any), stringsAsFactors = F)
  omic <- list(x=x,
               usedGenes=row.names(mut),
               namesCov="mut",
               method="binary",
               omicName="mut",
               eventThr=1)
  discrete <- x
  discrete$mut <- as.numeric(x$mut)
  mutationCount <- rowSums(mut)
  ordMutCount <- mutationCount[order(mutationCount, decreasing = T)]
  mo <- makeOmics(experiments = ExperimentList(mut=mut, mut2=mut), 
                  colData = dummy_colData(len=dim(mut)[2]), 
                  modelInfo = c("summarizeToBinaryEvents", "summarizeToNumberOfEvents"), 
                  specificArgs = list(binEvents = list(name = "mut"),
                                      countEvent = list(name = "mut", min_prop = 0.05)))
  binInfo <- extractSummaryFromBinary(omic, mo, n=3)
  expect_identical(row.names(binInfo$sigModule), head(names(ordMutCount),3))
  expect_identical(binInfo$discrete, discrete)
})


test_that("extractyInfo_dirBinary", {
  fake_cnv <- dummy_cnv_like_dataset()
  omic <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv), name="cnv")
  mo <- makeOmics(experiments = ExperimentList(cnv=fake_cnv, mut=dummy_mutation_like_dataset()), 
                  colData = dummy_colData(len=dim(fake_cnv)[2]), 
                  modelInfo = c("summarizeToBinaryDirectionalEvents", "summarizeToNumberOfEvents"),
                  specificArgs = list(countEvent=list(name="cnv"), countEvent2=list(name="mut")))
  dataModule <- createDataModule(omic, mo)
  cnvCountPos <- rowSums(dataModule >=2)
  cnvCountNeg <- rowSums(dataModule <=-2)
  ordCnvCountPos <- cnvCountPos[order(cnvCountPos, decreasing = T)]
  ordCnvCountNeg <- cnvCountNeg[order(cnvCountNeg, decreasing = T)]
  ordCnvGene <- unique(c(head(names(ordCnvCountPos), 3), head(names(ordCnvCountNeg), 3)))
  discrete2class <- data.frame(lapply(omic$x, as.numeric))
  row.names(discrete2class) <- row.names(omic$x)
  duobleBinInfo <- extractSummaryFromBinary(omic, mo, n=3)
  expect_identical(row.names(duobleBinInfo$sigModule), ordCnvGene)
  expect_identical(duobleBinInfo$discrete, discrete2class)
})





  