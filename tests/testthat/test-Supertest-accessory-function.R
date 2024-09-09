
library(graphite)
pathHierarchy <- downloadPathwayRelationFromReactome()
pathHierarchyGraph <- igraph::graph_from_data_frame(
  d = pathHierarchy, directed = TRUE)


pathwayName = c("Circadian Clock", "Signaling Pathways")
reactome <- pathways("hsapiens", "reactome")[pathwayName]


omicsClasses2pathways <- list(exp=pathwayName, 'cnv;exp'="Not a pathway")
test_that("generateWarn", {
  expect_warning(omicsClasses2fathers <- lapply(omicsClasses2pathways,
                                                annotePathwayToFather, 
                                                graphiteDB=reactome,
                                                hierarchy=pathHierarchyGraph))
})