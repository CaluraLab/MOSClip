test_that("generateWarn", {
  pathHierarchy <- downloadPathwayRelationFromReactome()
  pathHierarchyGraph <- igraph::graph_from_data_frame(
    d = pathHierarchy, directed = TRUE)
  
  pathwayName = c("Circadian Clock", "Signaling Pathways")
  reactome <- graphite::pathways("hsapiens", "reactome")[pathwayName]
  
  
  omicsClasses2pathways <- list(exp=pathwayName, 'cnv;exp'="Not a pathway")
  expect_warning(omicsClasses2fathers <- lapply(omicsClasses2pathways,
                                                annotePathwayToFather, 
                                                graphiteDB=reactome,
                                                hierarchy=pathHierarchyGraph))
  expect_warning(correspondence <- lapply(names(
    omicsClasses2pathways), match_pathway_to_fathers,
    omicsClasses2pathways=omicsClasses2pathways,
    omicsClasses2fathers=omicsClasses2fathers))
  
})