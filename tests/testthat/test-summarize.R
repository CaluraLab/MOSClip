test_that("summarizeToBinaryEvents", {
  mat <- dummy_mutation_like_dataset()
  res <- summarizeToBinaryEvents(mat, features = c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "logical")
  res <- summarizeToBinaryEvents(data=NULL)
  expect_null(res)
  res <- summarizeToBinaryEvents(mat, features = c("gene_x"))
  expect_null(res)
  })


test_that("summarizeToNumberOfEvents", {
  mat <- dummy_mutation_like_dataset()
  res <- summarizeToNumberOfEvents(mat, features = 
                                     c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "integer")
  res <- summarizeToNumberOfEvents(data=NULL)
  expect_null(res)
  res <- summarizeToNumberOfEvents(mat, features = c("gene_x"))
  expect_null(res)
  })


test_that("summarizeInCluster", {
  mat <- dummy_methylation_like_dataset()
  res <- summarizeInCluster(mat, features = c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "integer")
  res <- summarizeInCluster(data=NULL)
  expect_null(res)
  res <- summarizeInCluster(mat, features = c("gene_x"))
  expect_null(res)
})


test_that("summarizeWithPca", {
  mat <- dummy_expression_like_dataset()
  res <- summarizeWithPca(mat, features = c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "double")
  res <- summarizeInCluster(data=NULL)
  expect_null(res)
  res <- summarizeInCluster(mat, features = c("gene_x"))
  expect_null(res)
})


test_that("summarizeToNumberOfDirectionalEvents", {
  mat <- dummy_cnv_like_dataset()
  res <- summarizeToNumberOfDirectionalEvents(mat, features = 
                                                c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "integer")
  res <- summarizeToNumberOfDirectionalEvents(data=NULL)
  expect_null(res)
  res <- summarizeToNumberOfDirectionalEvents(mat, features = c("gene_x"))
  expect_null(res)
})


test_that("summarizeToBinaryDirectionalEvents", {
  mat <- dummy_cnv_like_dataset()
  res <- summarizeToBinaryDirectionalEvents(mat, features = 
                                              c("gene_1", "gene_2", "gene_3"))
  expect_type(res$x[[1]], "logical")
  res <- summarizeToBinaryDirectionalEvents(data=NULL)
  expect_null(res)
  res <- summarizeToBinaryDirectionalEvents(mat, features = c("gene_x"))
  expect_null(res)
})