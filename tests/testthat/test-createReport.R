test_that("multiPathwayReport", {
  res <- readRDS(test_path("fixtures", "resPathways.rds"))
  priority_to <- c("mut","exp")
  report <- multiPathwayReport(res, priority_to = priority_to)
  expect_false(is.null(report))
  covs <- unique(sapply(colnames(report), substr, 1, 3))[-1]
  expect_identical(covs[1:length(priority_to)], priority_to)
  expect_equal(nrow(report), 10)
  expect_equal(ncol(report), 9)
  expect_error(multiPathwayReport(res[[1]]),
               "A list of pathway results are expected.")
  res[[1]] <- list(1)
  expect_error(multiPathwayReport(res),
               "A list of pathway results are expected.")
})


test_that("multiPathwayModuleReport", {
  res <- readRDS(test_path("fixtures", "resModules.rds"))
  priority_to <- c("met","cnv")
  report <- multiPathwayModuleReport(res, priority_to = priority_to)
  expect_false(is.null(report))
  covs <- unique(sapply(colnames(report), substr, 1, 3))[-(1:3)]
  expect_identical(covs[1:length(priority_to)], priority_to)
  expect_s3_class(report, "data.frame")
  expect_equal(ncol(report), 11)
  expect_error(multiPathwayModuleReport(res[[1]]),
               "A list of pathway modules results are expected.")
  res[[1]] <- list(1)
  expect_error(multiPathwayModuleReport(res),
               "A list of pathway modules results are expected.")
})


