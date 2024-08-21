env <- test_env()
dummy_omicsSurv <- createFakeOmics("survival")
assign("dummy_omicsSurv", dummy_omicsSurv, envir = env)


test_that("Plots with two class object work", {
  # does not work w get()
  dummy_react <- readRDS(test_path("fixtures", "graph.rds"))
  dummy_omics <- createFakeOmics("two-classes")
  dummy_annot <- dummy_colData("two-classes")

  twoCPathwayObj <- multiOmicsTwoClassesPathwayTest(dummy_omics,
                                                    dummy_react[[1]],
                                                    dummy_annot,
                                                    baseFormula = "classes ~",
                                                    nullModel = "classes ~ 1")
  res.pathwayHeat <- plotPathwayHeat(
    twoCPathwayObj, sortBy = c("expPC2", "mut", "classes"),
    paletteNames = c(exp = "red", met = "green", mut = "blue", cnv = "yellow"),
    additionalAnnotations = dummy_annot,
    additionalPaletteNames = list(classes = "teal"),
    nrowsHeatmaps = 2, withSampleNames = F)

  expect_false(is.null(res.pathwayHeat))
  expect_identical(class(res.pathwayHeat), c("gg", "ggplot"))

  expect_error(plotPathwayHeat(
    twoCPathwayObj, sortBy = c("expPC2", "mut", "classes"),
    paletteNames = c(exp = "red", met = "green"),
    additionalAnnotations = dummy_annot,
    additionalPaletteNames = list(classes = "teal"),
    nrowsHeatmaps = 2, withSampleNames = F), "omic not found")

  res.pathwayReport <- plotMultiPathwayReport(list(twoCPathwayObj),
                                              MOcolors = c(exp = "red", met = "green", mut = "blue"),
                                              top = 40, fontsize = 10)

  expect_false(is.null(res.pathwayReport))
  expect_s4_class(res.pathwayReport, "HeatmapList")

  expect_error(
    plotMultiPathwayReport(
      list(twoCPathwayObj),
      MOcolors = c(exp = "red", met = "green", mut = "blue", cnv = "yellow"),
      top = 40, fontsize = 10), "Length of MOcolors differs")

  twoCModuleObj <- multiOmicsTwoClassesModuleTest(dummy_omics,
                                                  dummy_react[[1]],
                                                  dummy_annot,
                                                  baseFormula = "classes ~",
                                                  nullModel = "classes ~ 1")
  res.moduleHeat <- plotModuleHeat(
    twoCModuleObj, 2, paletteNames = c(
      exp = "red", met = "green", mut = "blue", cnv = "yellow"),
    additionalAnnotations = dummy_annot,
    additionalPaletteNames = list(classes = "teal"),
    sortBy = c("met2k", "expPC1", "classes"), withSampleNames = F)

  expect_false(is.null(res.moduleHeat))
  expect_identical(class(res.moduleHeat), c("gTree", "grob", "gDesc"))

  expect_error(plotModuleHeat(
    twoCModuleObj, 2, sortBy = c("expPC2", "mut", "classes"),
    paletteNames = c(exp = "red", met = "green"),
    additionalAnnotations = dummy_annot,
    additionalPaletteNames = list(classes = "teal"),
    nrowsHeatmaps = 2, withSampleNames = F), "omic not found")

  res.ModuleInGraph <- plotModuleInGraph(twoCModuleObj,
                                         moduleNumber = 3, legendLabels = c("expr", "methylation", "cnv"),
                                         paletteNames = c(exp = "red", met = "green", mut = "blue", cnv = "yellow"))

  expect_false(is.null(res.ModuleInGraph))
  expect_true(is.list(res.ModuleInGraph))
  expect_true(all(vapply(res.ModuleInGraph, function(x) !is.null(x),
                         logical(1))))

  expect_error(plotModuleInGraph(
    twoCModuleObj,
    moduleNumber = 3, legendLabels = c("expr", "methylation", "cnv"),
    paletteNames = c(
      exp = "orange", met = "green", mut = "blue", cnv = "yellow")),
    "paletteNames value is not allowed")
  expect_error(plotModuleInGraph(
    twoCModuleObj,
    moduleNumber = 3, legendLabels = c("expr", "methylation", "cnv"),
    paletteNames = c(exp = "red")), "Missing palette")

  res.moduleReport <- plotModuleReport(
    twoCModuleObj, MOcolors = c( exp = "red", met = "green", mut = "blue"))
  expect_false(is.null(res.moduleReport))
  expect_s4_class(res.moduleReport, "HeatmapList")
  expect_error(
    plotModuleReport(
      twoCModuleObj,
      MOcolors = c(exp = "red", met = "green", mut = "blue", cnv = "yellow")),
    "Length of MOcolors differs")
})

test_that("Kaplan Meier plots work", {
  dummy_react <- readRDS(test_path("fixtures", "graph.rds"))
  dummy_omicsSurv <- createFakeOmics("survival")

  res.survPathway <- multiOmicsSurvivalPathwayTest(
    dummy_omicsSurv, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = T)

  res.PathwayKM <- plotPathwayKM(
    res.survPathway,
    formula = "Surv(days, status) ~ met3k + expPC2",
    paletteNames = "Paired", envir=env)

  expect_false(is.null(res.PathwayKM))
  expect_identical(class(res.PathwayKM), c("ggsurvplot", "ggsurv", "list"))

  res.survModules <- multiOmicsSurvivalModuleTest(
    dummy_omicsSurv, dummy_react[[1]], survFormula = "Surv(days, status) ~",
    autoCompleteFormula = TRUE)

  res.ModuleKM <- plotModuleKM(
    res.survModules, 2,
    formula = "Surv(days, status) ~ met2k + expPC1",
    paletteNames = "Paired", risk.table = TRUE, size = 2, inYears = TRUE, envir=env)

  expect_false(is.null(res.ModuleKM))
  expect_identical(class(res.ModuleKM), c("ggsurvplot", "ggsurv", "list"))

  expect_error(plotModuleKM(
    res.survModules, 2,
    formula = "Surv(days, status) ~ met2k + expPC1",
    paletteNames = "Paired", risk.table = TRUE, size = 2, inYears = TRUE,
    additional_discrete = "recurrence_status", envir=env),
    "discrete variables were not found")

  expect_error(plotModuleKM(
    res.survModules, 2,
    formula = "Surv(days, status) ~ met2k + expPC1",
    paletteNames = "Paired", risk.table = TRUE, size = 2, inYears = TRUE,
    additional_continuous = "recurrence_days", envir=env),
    "continuous variables were not found")

})
