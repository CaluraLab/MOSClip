test_that("makeOmics", {
  fake_exp <- dummy_expression_like_dataset()
  fake_exp2 <- fake_exp[,-1]
  new_order <- sample(colnames(fake_exp))
  fake_ord <- fake_exp[, new_order]
  new_order_dup <- sample(row.names(fake_exp), replace = TRUE)
  fake_dup <- fake_exp[new_order_dup,]
  
  expect_message(makeOmics(experiments = ExperimentList(exp=fake_exp)), 
                 "Data and relative methods to analyze them must be equal in length.")
  expect_message(makeOmics(experiments = ExperimentList(exp=fake_exp, exp2=fake_exp),
                           modelInfo = c("summarizeWithPca", "test"),
                           specificArgs = "test"),
                 "data and specificArgs must be equal in length.")
  expect_message(makeOmics(experiments = ExperimentList(exp=fake_exp),
                           modelInfo = "test",
                           specificArgs = "test"),
                 "test modelInfo not found in method. Try availableOmicMethods.")
  expect_message(makeOmics(experiments = ExperimentList(exp=fake_exp, exp2=fake_exp2),
                           modelInfo = c("summarizeWithPca", "summarizeWithPca"),
                           specificArgs = list("pca", "pca")),
                 "Mismatch in sample numbers")
  expect_message(makeOmics(experiments = ExperimentList(exp=fake_exp, exp2=fake_dup),
                           modelInfo = c("summarizeWithPca", "summarizeWithPca"),
                           specificArgs = list("pca", "pca")),
                 "Duplicated row.names found in omics ")
  expect_message(mo <- makeOmics(experiments = ExperimentList(exp=fake_exp, exp2=fake_ord),
                                 modelInfo = c("summarizeWithPca", "summarizeWithPca"),
                                 specificArgs = list("pca", "pca")),
                 "Samples order mismatch")
  expect_null(mo)
  expect_s4_class(fake_mo(), "Omics")

  exp0 <- matrix(0, nrow=5, ncol=200)
  colnames(exp0) <- paste0("P_", seq_len(NCOL(exp0)))
  expect_error(makeOmics(experiments = ExperimentList(exp=exp0),
               modelInfo=c("summarizeWithPca"),
               specificArgs=list("pca")), 
               "Some experiment matrices contains only 0 or NA values: ")
})

