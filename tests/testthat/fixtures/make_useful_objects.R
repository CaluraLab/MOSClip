

if (file.exists("reactome-entrez.RData")) {
  load("reactome-entrez.RData")
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
}


if (file.exists("tests/testthat/fixtures/graph.rds")) {
  graph <- readRDS("tests/testthat/fixtures/graph.rds")
} else {saveRDS(reactome[[1]], file="tests/testthat/fixtures/graph.rds")}



