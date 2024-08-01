# Creating reactome object for test files

react_dummy <- graphite::pathways("hsapiens", "reactome")
react_dummy <- graphite::convertIdentifiers(reactome, "entrez")
nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(pseudoExpNorm)))})
reactSmallDummy <- reactome[nodesLength >= 20 & nodesLength <= 100]

saveRDS(reactSmallDummy[1], file = "reactSmallDummy.rds")