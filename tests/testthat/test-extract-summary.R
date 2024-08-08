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
  summ <- summarizeToBinaryEvents(mut, features = row.names(mut), name = "mut")
  expect_identical(summ$x, x)
  discrete <- x
  discrete$mut <- as.numeric(x$mut)
  mutationCount <- rowSums(mut)
  ordMutCount <- mutationCount[order(mutationCount, decreasing = T)]
  mo <- fake_mo(omics="mut")
  binInfo <- extractSummaryFromBinary(omic, mo, n=3)
  expect_identical(row.names(binInfo$sigModule), head(names(ordMutCount),3))
  expect_identical(binInfo$discrete, discrete)
})


test_that("extractyInfo_dirBinary", {
  fake_cnv <- dummy_cnv_like_dataset()
  omic <- summarizeToBinaryDirectionalEvents(fake_cnv, row.names(fake_cnv), name="cnv")
  mo <- fake_mo(omics="cnv")
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


test_that("extractyInfo_Cluster", {
  fake_met <- dummy_methylation_like_dataset()
  omic <- summarizeInCluster(fake_met, row.names(fake_met), name="met")
  mo <- fake_mo(omics="met")
  clsInfo <- extractSummaryFromCluster(omic, mo, n=3)
  dataModule <- createDataModule(omic, mo)
  KMsigMat <- KWtest(t(dataModule), omic$x[,1])
  involved <- head(KMsigMat, 3)
  expect_setequal(row.names(clsInfo$subset), row.names(involved))
})
  

test_that("extractyInfo_Pca", {
  fake_exp <- dummy_expression_like_dataset()
  omic <- summarizeWithPca(fake_exp, row.names(fake_exp), name="exp", maxPCs = 1)  
  mo <- fake_mo(omics="exp", type="survival")
  moview <- createMOMView(mo, row.names(fake_exp))
  cox <- createCoxObj(mo@colData, moview)
  pcaInfo <- extractSummaryFromPCA(omic, mo, cox, "survival", loadThr = 0.6, atleast = 1)
  lds <- omic$loadings[order(abs(omic$loadings), decreasing = T),]
  topgenes <- names(lds[abs(lds) >= 0.6])
  if (length(topgenes)==0) {topgenes <- names(lds)[1]}
  expect_setequal(topgenes, row.names(pcaInfo$subset))
  expect_true(length(row.names(pcaInfo$subset)) >= 1)
  pcaInfo <- extractSummaryFromPCA(omic, mo, cox, "survival", loadThr = 1, atleast = 0)
  expect_null(pcaInfo$subset)
})  


test_that("createDiscreteClasses", {
  mo <- fake_mo(omics="exp", type="survival")
  moview <- createMOMView(mo, row.names(assay(mo, "exp")))
  cox <- createCoxObj(mo@colData, moview)
  expect_error(createDiscreteClasses(cox, covs="mut", analysis="survival"), 
               regexp=" not in coxObj")
  covs <- colnames(cox)[-c(1,2)]
  expect_error(createDiscreteClasses(cox, covs=covs, analysis="class"), 
               regexp="Type of analysis .* not valid. Check results object")
  sc <- createDiscreteClasses(cox, covs, analysis = "survival")
  expect_setequal(colnames(sc), colnames(cox))
  expect_setequal(unique(unlist(apply(sc[,covs], 2, unique, simplify = FALSE))), c("low", "high"))
  expect_s3_class(sc, "surv_categorize")
  sc <- createDiscreteClasses(cox, covs, analysis = "twoClass")
  expect_s3_class(sc, "data.frame")
  cox[,"expPC2"] <- rep(0, 200)
  expect_error(createDiscreteClasses(cox, covs=covs, analysis="class"), 
               regexp="minprop 0.1 is too high. Try a smaller one")
})


test_that("extractyInfo_NumberOfEvents", {
  fake_mut <- dummy_mutation_like_dataset()
  omic <- summarizeToNumberOfEvents(fake_mut, row.names(fake_mut), name="mut")
  mo <- fake_mo(omics="mut", type="survival")
  moview <- createMOMView(mo, row.names(fake_mut))
  cox <- createCoxObj(mo@colData, moview)
  countInfo <- extractSummaryFromNumberOfEvents(omic, mo, cox, analysis = "survival")
  mut_sum <- rowSums(fake_mut)
  mut_sum <- mut_sum[order(mut_sum, decreasing = TRUE)]
  expect_setequal(names(mut_sum)[1:3], row.names(countInfo$subset))
})  



