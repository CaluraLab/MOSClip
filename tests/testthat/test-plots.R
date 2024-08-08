env <- test_env()
graph <- readRDS(test_path("fixtures", "graph.rds"))
mo <- fake_mo(genes=nodes(graph), type="survival", omics="exp")
res <- multiOmicsSurvivalModuleTest(mo, graph)
assign("graph", graph, envir = env)
assign("mo", mo, envir = env)
assign("res", res, envir = env)


test_that("multiOmicsSurvialTest", {
  res <- multiOmicsSurvivalModuleTest(mo, graph)
  expect_s4_class(res, "MultiOmicsModules")
  p <- plotModuleKM(res, 1, formula="Surv(days, status) ~ expPC1", envir=env)
  expect_s3_class(p, "ggsurvplot")
  expect_false(is.null(p))
})






