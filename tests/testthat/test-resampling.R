test_that("resamplingTwoClasses works", {
  graphs <- readRDS(test_path("fixtures", "graph_list.rds"))
  gSubset<-graphs[c(2,4,6,10)]
  nodes <- sapply(gSubset, nodes)
  genes <- unlist(unname(lapply(nodes, function(x) sample(x, 5))))
  genes <- genes[!duplicated(genes)]
  mo <- fake_mo(genes=genes, type="two-classes")
  res <- resamplingTwoClasses(mo, classAnnot = mo@colData, pathdb = graphs,
                              pathwaySubset = names(gSubset), nperm=5)
  expect_length(res, 5)
  expect_true(all(sapply(res, is.data.frame)))
  dims <- lapply(res, dim)
  expect_true(all(sapply(dims, function(x) all(x==dims[[1]]))))
})


test_that("resamplingPathwayTwoClasses works", {
  graphs <- readRDS(test_path("fixtures", "graph_list.rds"))
  gSubset<-graphs[c(2,4,6,10)]
  nodes <- sapply(gSubset, nodes)
  genes <- unlist(unname(lapply(nodes, function(x) sample(x, 5))))
  genes <- genes[!duplicated(genes)]
  mo <- fake_mo(genes=genes, type="two-classes")
  res <- resamplingPathwayTwoClasses(mo, classAnnot = mo@colData,
                                     pathdb = graphs,
                                     pathwaySubset = names(gSubset), nperm=5)
  expect_length(res, 5)
  expect_true(all(sapply(res, is.data.frame)))
  dims <- lapply(res, dim)
  expect_true(all(sapply(dims, function(x) all(x==dims[[1]]))))
})


test_that("resamplingSurvival works", {
  graphs <- readRDS(test_path("fixtures", "graph_list.rds"))
  gSubset<-graphs[c(2,4,6,10)]
  nodes <- sapply(gSubset, nodes)
  genes <- unlist(unname(lapply(nodes, function(x) sample(x, 5))))
  genes <- genes[!duplicated(genes)]
  mo <- fake_mo(genes=genes, type="survival")
  res <- resamplingSurvival(mo, pathdb = graphs,
                            pathwaySubset = names(gSubset), nperm=5)
  expect_length(res, 5)
  expect_true(all(sapply(res, is.data.frame)))
  dims <- lapply(res, dim)
  expect_true(all(sapply(dims, function(x) all(x==dims[[1]]))))
})


test_that("resamplingPathwaySurvival works", {
  graphs <- readRDS(test_path("fixtures", "graph_list.rds"))
  gSubset<-graphs[c(2,4,6,10)]
  nodes <- sapply(gSubset, nodes)
  genes <- unlist(unname(lapply(nodes, function(x) sample(x, 5))))
  genes <- genes[!duplicated(genes)]
  mo <- fake_mo(genes=genes, type="survival")
  res <- resamplingPathwaySurvival(mo, pathdb = graphs,
                                   pathwaySubset = names(gSubset), nperm=5)
  expect_length(res, 5)
  expect_true(all(sapply(res, is.data.frame)))
  dims <- lapply(res, dim)
  expect_true(all(sapply(dims, function(x) all(x==dims[[1]]))))
})


test_that("selectStablePathwaysModules works", {
  graphs <- readRDS(test_path("fixtures", "graph_list.rds"))
  gSubset <- graphs[c(2,4,6,10)]
  nodes <- sapply(gSubset, nodes)
  genes <- unlist(unname(lapply(nodes, function(x) sample(x, 5))))
  genes <- genes[!duplicated(genes)]
  mo <- fake_mo(genes=genes, type="survival")
  modules <- lapply(gSubset, function(g) {
    multiOmicsSurvivalModuleTest(
      mo, g, useThisGenes = row.names(mo@ExperimentList$exp))})
  summ <- multiPathwayModuleReport(modules)
  res <- resamplingSurvival(mo, pathdb = graphs,
                            pathwaySubset = names(gSubset), nperm=5)
  stableModulesSummary <- selectStablePathwaysModules(perms = res, 
                                                      moduleSummary = summ, 
                                                      success = 2)
  resampligSuccessCount <- getPathwaysModulesSuccess(perms = res, 
                                                     moduleSummary = summ, 
                                                     thr=0.5)
  moduleSummary <- addResamplingCounts(summ, resampligSuccessCount)
  expect_identical(colnames(moduleSummary)[NCOL(moduleSummary)], 
                   "resamplingCount")
})
