

if (file.exists("reactome-entrez.RData")) {
  load("reactome-entrez.RData")
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
}


if (file.exists("tests/testthat/fixtures/reactSmallDummy.rds")) {
  graph <- readRDS("tests/testthat/fixtures/reactSmallDummy.rds")
} else {
  graph <- reactome["Intrinsic Pathway for Apoptosis"]
  saveRDS(graph, file="tests/testthat/fixtures/reactSmallDummy.rds")}


if (file.exists("tests/testthat/fixtures/graph_list.rds")) {
  graph_list <- readRDS("tests/testthat/fixtures/graph_list.rds")
} else {
  graph_list <- reactome[1:10]
  saveRDS(reactome[1:10], file="tests/testthat/fixtures/graph_list.rds")}





