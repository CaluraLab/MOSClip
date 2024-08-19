

if (file.exists("reactome-entrez.RData")) {
  load("reactome-entrez.RData")
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
}


if (file.exists("tests/testthat/fixtures/graph.rds")) {
  graph <- readRDS("tests/testthat/fixtures/graph.rds")
} else {
  graph <- reactome[[1]]
  saveRDS(reactome[[1]], file="tests/testthat/fixtures/graph.rds")}


if (file.exists("tests/testthat/fixtures/graph_list.rds")) {
  graph_list <- readRDS("tests/testthat/fixtures/graph_list.rds")
} else {
  graph_list <- reactome[1:10]
  saveRDS(reactome[1:10], file="tests/testthat/fixtures/graph_list.rds")}


if (file.exists("tests/testthat/fixtures/resPathways.rds")) {
  res <- readRDS("tests/testthat/fixtures/resPathways.rds")
} else {
  genes <- unique(unlist(unname(lapply(graph_list, nodes))))
  mo <- fake_mo(genes = genes)
  res <- lapply(graph_list, function(x) {multiOmicsSurvivalPathwayTest(mo, x)})
  saveRDS(res, file="tests/testthat/fixtures/resPathways.rds")
}


if (file.exists("tests/testthat/fixtures/resModules.rds")) {
  res <- readRDS("tests/testthat/fixtures/resModules.rds")
} else {
  genes <- unique(unlist(unname(lapply(graph_list, nodes))))
  mo <- fake_mo(genes = genes, type = "two-classes")
  res <- lapply(graph_list, function(x) {multiOmicsTwoClassesModuleTest(mo, x, mo@colData)})
  saveRDS(res, file="tests/testthat/fixtures/resModules.rds")
}



