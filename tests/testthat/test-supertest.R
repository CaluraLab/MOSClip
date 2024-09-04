test_that("superTest works", {
  whole_reactome <- readRDS(test_path("fixtures", "reactome-entrez.rds"))
  dummy_react_bigger <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
  dummy_omics <- createFakeOmics("two-classes")
  dummy_annot <- dummy_colData("two-classes")
  
  twoCModuleObj <- lapply(dummy_react_bigger, function(g) {
    print(g@title)
    fcl = multiOmicsTwoClassesModuleTest(dummy_omics,
                                         g,
                                         dummy_annot,
                                         baseFormula = "classes ~",
                                         nullModel = "classes ~ 1")
    fcl
  })
  
  moduleSummary <- multiPathwayModuleReport(twoCModuleObj)
  
  moduleSum_notnum <- moduleSummary
  moduleSum_notnum$mut <- rep("abc", nrow(moduleSum_notnum))
  expect_error(computeOmicsIntersections(moduleSum_notnum,
                                         pvalueThr = 0.05, 
                                         zscoreThr = 1,
                                         excludeColumns = 
                                           c("pathway", "module")),
               "columns are not numeric")
  pvalSum <- pvalueSummary(moduleSummary,
                           excludeColumns =  c("pathway", "module"),
                           as.list = FALSE)
  expect_identical(class(pvalSum), c("matrix", "array"))
  expect_false(all(vapply(pvalSum, is.null, logical(1))))
  moduleSum_notpval <- moduleSummary
  moduleSum_notpval$mut <- sample(2:5, size = nrow(moduleSum_notpval),
                                     replace = TRUE)
  moduleSum_notpval$mut <- as.numeric(moduleSum_notpval$mut)
  expect_error(computeOmicsIntersections(moduleSum_notpval,
                                         pvalueThr = 0.05, 
                                         zscoreThr = 1,
                                         excludeColumns = 
                                           c("pathway", "module")),
               "columns are not pvalues")  
  
  mi <- runSupertest(moduleSummary, pvalueThr = 0.05, zscoreThr = 1,
                     excludeColumns = c("pathway", "module"), plot="noplot")
  expect_identical(class(mi), "data.frame")
  expect_false(all(vapply(mi, is.null, logical(1))))
  
  expect_error(
    runSupertest(moduleSummary, pvalueThr = 0.05, zscoreThr = 1,
                 excludeColumns = c("pathway", "module")),
               NA)
  expect_error(
    runSupertest(moduleSummary, pvalueThr = 0.05, zscoreThr = 1,
                 excludeColumns = c("pathway", "module"), plot = "landscape"),
               NA)
  expect_error(
    runSupertest(
      moduleSummary, pvalueThr = 0.05, zscoreThr = 1, resampligThr = 60,
      excludeColumns = c("pathway", "module"), plot = "noplot"),
    "resamplingCount column not found")
  
  fake_resamplingCount <- sample(0:100, size = nrow(moduleSummary),
                                 replace = TRUE)
  names(fake_resamplingCount) <- row.names(moduleSummary)
  
  moduleSummary_resampled <- addResamplingCounts(
    moduleSummary, fake_resamplingCount)
  
  mi_resampled <- runSupertest(moduleSummary_resampled, pvalueThr = 0.05,
                               zscoreThr = 1, resampligThr = 60,
                               excludeColumns = c(
                                 "pathway", "module", "resamplingCount"),
                               plot="noplot")
  expect_identical(class(mi_resampled), "data.frame")
  expect_false(all(vapply(mi_resampled, is.null, logical(1))))
  
  expect_error(
    runSupertest(moduleSummary_resampled, pvalueThr = 0.05, zscoreThr = 1,
                 resampligThr = 60,
                 excludeColumns = c("pathway", "module", "resamplingCount")),
    NA)
  expect_error(
    runSupertest(moduleSummary_resampled, pvalueThr = 0.05, zscoreThr = 1,
                 resampligThr = 60,
                 excludeColumns = c("pathway", "module", "resamplingCount"),
                 plot = "landscape"),
    NA)
})

