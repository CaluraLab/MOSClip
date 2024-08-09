test_that("Survival Results object works", {
  
  dummy_react <- readRDS(test_path("fixtures", "reactSmallDummy.rds"))
  dummy_Omics <- createFakeOmics("survival")
  dummy_survcolData <- dummy_colData("survival")
  
  res.survPathway <- multiOmicsSurvivalPathwayTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T, include_from_annot = TRUE)
  
  expect_s4_class(res.survPathway, "MultiOmicsPathway")
  expect_true(!is.null(res.survPathway@pvalue) &
                is.numeric(res.survPathway@pvalue))
  expect_true(!is.null(res.survPathway@zlist) & is.numeric(
    res.survPathway@zlist))
  expect_true(is.list(res.survPathway@pathView))
  expect_true(all(vapply(res.survPathway@pathView, function(x) !is.null(x),
                         logical(1))))
  
  res.survModule <- multiOmicsSurvivalModuleTest(
    dummy_Omics, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T)
  
  expect_s4_class(res.survModule, "MultiOmicsModules")
  expect_true(!is.null(res.survModule@alphas) & is.numeric(res.survModule@alphas))
  expect_true(!is.null(res.survModule@zlists) & is.list(res.survModule@zlists))
  
  graph <- convertPathway(dummy_react[[1]], NULL)
  cliques <- extractCliquesFromDag(graph)
  
  expect_true(is.list(res.survModule@modulesView))
  expect_identical(length(cliques), length(res.survModule@modulesView))
  expect_true(all(vapply(res.survModule@modulesView, function(x) !is.null(x),
                         logical(1))))
})
