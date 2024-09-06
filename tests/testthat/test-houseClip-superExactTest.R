test_that("superExactTest works", {
  skip()
  whole_reactome <- readRDS(test_path("fixtures", "reactome-entrez.rds"))
  dummy_react_bigger <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  pathHierarchy <- downloadPathwayRelationFromReactome()
  pathHierarchyGraph <- igraph::graph_from_data_frame(d = pathHierarchy,
                                                      directed = TRUE)
  dummy_omics <- fake_mo(genes=genes, type="two-classes")
  dummy_annot <- dummy_colData(type="two-classes")
  
  twoCModuleObj <- lapply(dummy_react_bigger, function(g) {
    fcl = multiOmicsTwoClassesModuleTest(dummy_omics,
                                         g,
                                         dummy_annot,
                                         baseFormula = "classes ~",
                                         nullModel = "classes ~ 1")
    fcl
  })

  moduleSummary <- multiPathwayModuleReport(twoCModuleObj)
  
  expect_error(
    computeOmicsIntersections(moduleSummary,pvalueThr = 0.05, zscoreThr = 0.05,
                               excludeColumns =  c("pathway", "module")),
    "significant z score")
  expect_error(
    computeOmicsIntersections(moduleSummary, pvalueThr = 0.05, zscoreThr = 1,
                              excludeColumns =c("pathway", "module")),
    NA)
  expect_error(
    computeOmicsIntersections(moduleSummary, pvalueThr = 1.8e-20, zscoreThr = 1,
                              excludeColumns =c("pathway", "module")),
    paste0("One or more omics do not have a significant z score for any of ", 
           "their modules. Try increasing the z score threshold."))
  expect_error(
    computeOmicsIntersections(
      moduleSummary, pvalueThr = 0.05, zscoreThr = 1, resampligThr = 60,
      excludeColumns = c("pathway", "module")),
    "resamplingCount column not found")
  
  omicsClasses2pathways <- computeOmicsIntersections(moduleSummary,
                                                     pvalueThr = 0.05, 
                                                     zscoreThr = 1,
                                                     excludeColumns = 
                                                       c("pathway", "module"))
  omicsClasses2pathways <- lapply(omicsClasses2pathways,
                                  stripModulesFromPathways)
  
  omicsClasses2fathers <- lapply(
    omicsClasses2pathways, annotePathwayToFather,
    graphiteDB = whole_reactome, 
    hierarchy = pathHierarchyGraph)
  
  expect_identical(class(omicsClasses2fathers), "list")
  expect_false(all(vapply(omicsClasses2fathers, is.null, logical(1))))
  
  freqDataframe <- computeFreqs(omicsClasses2fathers)
  
  expect_error(
    plotFrequencies(
      freqDataframe, minSize = 4, maxSize = 12, width = 10, lineSize = 3),
    NA)
  
  colors <- c(`exp;met;mut` = "darkorchid")
  expect_error(
    plotFrequencies(
      freqDataframe, minSize = 4, maxSize = 12, width = 10, lineSize = 3,
      manualColors = colors),
    NA)
  
  fake_resamplingCount <- sample(0:100, size = nrow(moduleSummary),
                                 replace = TRUE)
  names(fake_resamplingCount) <- row.names(moduleSummary)
  
  moduleSummary_resampled <- addResamplingCounts(
    moduleSummary, fake_resamplingCount)
  
  omicsClasses2pathways_r <- computeOmicsIntersections(
    moduleSummary_resampled, pvalueThr = 0.05, zscoreThr = 1, resampligThr = 60,
    excludeColumns = c("pathway", "module", "resamplingCount"))
  omicsClasses2pathways_r <- lapply(omicsClasses2pathways_r,
                                  stripModulesFromPathways)
  omicsClasses2fathers_r <- lapply(
    omicsClasses2pathways_r, annotePathwayToFather,
    graphiteDB = whole_reactome, 
    hierarchy = pathHierarchyGraph)
  freqDataframe <- computeFreqs(omicsClasses2fathers_r)
  
  expect_error(
    plotFrequencies(
      freqDataframe, minSize = 4, maxSize = 12, width = 10, lineSize = 3),
    NA)
  
  colors <- c(`exp;met;mut` = "darkmagenta")
  expect_error(
    plotFrequencies(
      freqDataframe, minSize = 4, maxSize = 12, width = 10, lineSize = 3,
      manualColors = colors),
    NA)
})
