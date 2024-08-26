test_that("Two Class Pathway Test works", {
  genes <- paste0("ENTREZID:", c(10000, 90427, 5366, 23368, 637, 3002, 79792,
                                 840, 5414, 5599, 8655, 140735, 708))
  dummy_react <- readRDS(test_path("fixtures", "reactSmallDummy.rds"))
  dummy_Omics <- fake_mo(genes=genes, type="two-classes")
  dummy_annot <- dummy_colData(type="two-classes")

  twoCPathwayObj <- multiOmicsTwoClassesPathwayTest(dummy_Omics,
                                                    dummy_react[[1]],
                                                    dummy_annot,
                                                    baseFormula = "classes ~",
                                                    nullModel = "classes ~ 1")

  expect_s4_class(twoCPathwayObj, "MultiOmicsPathway")
  expect_true(!is.null(twoCPathwayObj@pvalue) &
                is.numeric(twoCPathwayObj@pvalue))
  expect_true(!is.null(twoCPathwayObj@zlist) & is.numeric(twoCPathwayObj@zlist))
  expect_true(is.list(twoCPathwayObj@pathView))
  expect_true(all(vapply(twoCPathwayObj@pathView, function(x) !is.null(x),
                         logical(1))))

  twoCModuleObj <- multiOmicsTwoClassesModuleTest(dummy_Omics,
                                                  dummy_react[[1]],
                                                  dummy_annot,
                                                  baseFormula = "classes ~",
                                                  nullModel = "classes ~ 1")

  expect_s4_class(twoCModuleObj, "MultiOmicsModules")
  expect_true(!is.null(twoCModuleObj@alphas) & is.numeric(twoCModuleObj@alphas))
  expect_true(!is.null(twoCModuleObj@zlists) & is.list(twoCModuleObj@zlists))

  graph <- convertPathway(dummy_react[[1]], NULL)
  cliques <- extractCliquesFromDag(graph)

  expect_true(is.list(twoCModuleObj@modulesView))
  expect_identical(length(cliques), length(twoCModuleObj@modulesView))
  expect_true(all(vapply(twoCModuleObj@modulesView, function(x) !is.null(x),
                         logical(1))))
})
