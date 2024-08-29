# Creating reactome object for test files
load("inst/TCGA-OV-hg38-preprocessed-externalSurv.RData")
expression <- rnaSeq
row.names(expression) <- paste0("ENTREZID:", row.names(expression))

if (file.exists("reactome-entrez.RData")) {
  whole_reactome <- readRDS(test_path("fixtures", "reactome-entrez"))
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
  saveRDS(reactome, file = "test/testthat/fixtures/reactome-entrez.rds")
}

# Creating smaller version for smaller tests
nodesLength <- sapply(react_dummy, function(g) {
  length(intersect(graphite::nodes(g), row.names(expression)))})
reactSmallDummy <- react_dummy[nodesLength >= 20 & nodesLength <= 100]

saveRDS(reactSmallDummy[1], file = "reactSmallDummy.rds")

# Saving reactome object with 5 pathways for frequencies plot
saveRDS(reactSmallDummy[1:5],
        file = "tests/testthat/fixtures/reactBiggerDummy.rds")

