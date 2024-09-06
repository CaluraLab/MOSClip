test_that("Survival Results object for pathways works", {
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  dummy_react <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
  dummy_Omics <- fake_mo(genes=genes, type="survival")
  dummy_survcolData <- dummy_colData(type="survival")

  res.survPathway <- multiOmicsSurvivalPathwayTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = TRUE, include_from_annot = TRUE)

  expect_s4_class(res.survPathway, "MultiOmicsPathway")
  expect_true(!is.null(res.survPathway@pvalue) &
                is.numeric(res.survPathway@pvalue))
  expect_true(!is.null(res.survPathway@zlist) & is.numeric(
    res.survPathway@zlist))
  expect_true(is.list(res.survPathway@pathView))
  expect_true(all(vapply(res.survPathway@pathView, function(x) !is.null(x),
                         logical(1))))
  
  expect_error(multiOmicsSurvivalPathwayTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = TRUE, include_from_annot = TRUE, 
    useThisGenes = c("test1", "test2")), "There is no nodes on the graph.")
  
  expect_s4_class(multiOmicsSurvivalPathwayTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T, robust = TRUE), "MultiOmicsPathway")
})


test_that("Survival Results object  for modules works", {
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  dummy_react <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
  dummy_Omics <- fake_mo(genes=genes, type="survival")
  dummy_survcolData <- dummy_colData(type="survival")

  res.survModule <- multiOmicsSurvivalModuleTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T)

  expect_s4_class(res.survModule, "MultiOmicsModules")
  expect_true(
    !is.null(res.survModule@alphas) & is.numeric(res.survModule@alphas))
  expect_true(!is.null(res.survModule@zlists) & is.list(res.survModule@zlists))

  graph <- convertPathway(dummy_react[[1]], NULL)
  cliques <- extractCliquesFromDag(graph)

  expect_true(is.list(res.survModule@modulesView))
  expect_identical(length(cliques), length(res.survModule@modulesView))
  expect_true(all(vapply(res.survModule@modulesView, function(x) !is.null(x),
                         logical(1))))
  
  res.survModule <- multiOmicsSurvivalModuleTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T, robust = TRUE)
  
  expect_s4_class(res.survModule, "MultiOmicsModules")
  expect_true(
    !is.null(res.survModule@alphas) & is.numeric(res.survModule@alphas))
  expect_true(!is.null(res.survModule@zlists) & is.list(res.survModule@zlists))
  
  graph <- convertPathway(dummy_react[[1]], NULL)
  cliques <- extractCliquesFromDag(graph)
  
  expect_true(is.list(res.survModule@modulesView))
  expect_identical(length(cliques), length(res.survModule@modulesView))
  expect_true(all(vapply(res.survModule@modulesView, function(x) !is.null(x),
                         logical(1))))
  
  expect_error(multiOmicsSurvivalModuleTest(
    dummy_Omics, genes, survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T), "Module test cannot handle gene list")
  
  expect_error(multiOmicsSurvivalModuleTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T, useThisGenes = c("test1", "test2")), 
    paste0("There is no intersection between expression feature ",
    "names and the node names in the graph"))
  
  dummy_Omics@specificArgs$exp$method = "topological"
  expect_error(res.survModule <- multiOmicsSurvivalModuleTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T), "Invalid method for module analysis: topological")
})


test_that("Survival Test works with additional covariates", {
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  dummy_react <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
  dummy_omics<-fake_mo(genes=genes, type="survival")
  moview<-createMOMView(dummy_omics, genes)
  values <- seq(-1, 1, length.out = 21)
  fake_column <- data.frame(fake=sample(values, size=200, replace=TRUE), 
                               row.names = row.names(dummy_omics@colData))
  dummy_omics@colData <- DataFrame(cbind(fake_column, dummy_omics@colData))
  createCoxObj(dummy_omics@colData, moview)
  res.survModule <- multiOmicsSurvivalModuleTest(
    dummy_omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T, include_from_annot = TRUE)
  expect_s4_class(res.survModule, "MultiOmicsModules")
  expect_true(
    !is.null(res.survModule@alphas) & is.numeric(res.survModule@alphas))
  expect_true(!is.null(res.survModule@zlists) & is.list(res.survModule@zlists))
  })
