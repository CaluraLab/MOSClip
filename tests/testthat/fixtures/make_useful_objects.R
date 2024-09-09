
if (file.exists("tests/testthat/fixtures/reactBiggerDummy.rds")) {
    graph <- readRDS("tests/testthat/fixtures/reactBiggerDummy.rds")
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
  graphs <- c("Intrinsic Pathway for Apoptosis", 
              "PI3K Cascade",
              "Translesion synthesis by Y family DNA polymerases bypasses lesions on DNA template",
              "Recognition of DNA damage by PCNA-containing replication complex",
              "Recognition and association of DNA glycosylase with site containing an affected pyrimidine")
  dummy_reactome <- reactome[graphs]
  saveRDS(dummy_reactome, file = "tests/testthat/fixtures/reactBiggerDummy.rds")
}






