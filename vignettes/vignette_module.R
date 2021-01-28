library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EDASeq)
library(graphite)
library(houseOfClipUtility)
library(biocmosclip)
library(MultiAssayExperiment)

load("vignettes/TCGA-OV-hg38-preprocessed-externalSurv.RData")

survAnnotations <- fup$pfs

expression <- rnaSeq
row.names(expression) <- paste0("ENTREZID:", row.names(expression))
mutation <- mutations$data
row.names(mutation) <- paste0("ENTREZID:", row.names(mutation))
names(metClustValues$eMap) <- paste0("ENTREZID:", row.names(metClustValues$eMap))
row.names(cnv) <- paste0("ENTREZID:", row.names(cnv))

survAnnot <- na.omit(survAnnotations)
patients <- row.names(survAnnot)
patients <- intersect(patients, colnames(expression))
patients <- intersect(patients, colnames(metClustValues$MET_Cancer_Clustered))
patients <- intersect(patients, colnames(mutation))
patients <- intersect(patients, colnames(cnv))

survAnnot <- survAnnot[patients, , drop = F]
dirname = "vignettes"

if (!file.exists("vignettes/survAnnot.RData")) {
  file = paste0(dirname, "/survAnnot-", as.character(Sys.Date()), ".RData")
  link = "survAnnot.RData"
  save(survAnnot, file = file)
  file.symlink(file, link)
}

expression <- expression[, patients, drop = F]
keep = apply(expression >= 100, 1, any)
expNorm <- betweenLaneNormalization(expression[keep, , drop = F], which = "upper")
pseudoExpNorm <- log2(expNorm + 1)

methylation <- metClustValues
methylation$MET_Cancer_Clustered <- methylation$MET_Cancer_Clustered[, patients, 
                                                                     drop = F]

mutation <- mutation[, patients, drop = F]
cnv <- cnv[, patients, drop = F]

if (file.exists("vignettes/reactome-entrez.RData")) {
  load("vignettes/reactome-entrez.RData")
} else {
  reactome <- pathways("hsapiens", "reactome")
  reactome <- convertIdentifiers(reactome, "entrez")
  file = paste0(dirname, "/reactome-entrez-", as.character(Sys.Date()), ".RData")
  link = "vignettes/reactome-entrez.RData"
  save(reactome, file = file)
  file.symlink(file, link)
}

nodesLength <- sapply(reactome, function(g) {
  length(intersect(graphite::nodes(g), row.names(pseudoExpNorm)))
})
reactSmall <- reactome[nodesLength >= 20 & nodesLength <= 100]
reactHuge <- reactome[nodesLength >= 10]
reactUsed <- reactSmall[c(2,10,20,43)]

pathHierarchy <- houseOfClipUtility::downloadPathwayRelationFromReactome()
pathHierarchyGraph <- igraph::graph.data.frame(d = pathHierarchy, directed = TRUE)


multiOmics <- Omics(experiments = ExperimentList(expr = pseudoExpNorm, met = methylation$MET_Cancer_Clustered, 
                                                 mut = mutation, cnv = cnv),
                    colData = survAnnot,
                    modelInfo = c("summarizeWithPca", "summarizeInCluster", 
                                  "summarizeToNumberOfEvents", "summarizeToNumberOfDirectionalEvents"),
                    specificArgs = list(pcaArgs = list(name = "exp",
                                                       shrink = "FALSE", method = "sparse", maxPCs = 3),
                                        clusterArgs = list(name = "met",  dict = methylation$eMap, max_cluster_number = 3),
                                        countEvent = list(name = "mut", min_prop = 0.05), cnvAgv = list(name = "cnv", min_prop = 0.05)))

genesToConsider <- row.names(pseudoExpNorm)

if (file.exists("vignettes/multiOmicsReactome.RData")) {
  load("vignettes/multiOmicsReactome.RData")
} else {
  multiOmicsReactome <- lapply(reactUsed, function(g) {
    print(g@title)  # uncomment this to see the progression along pathways
    set.seed(1234)
    fcl = multiOmicsSurvivalModuleTest(multiOmics, g, useThisGenes = genesToConsider)
    fcl
  })
  file = paste0(dirname, "/multiOmicsReactome-", as.character(Sys.Date()), ".RData")
  link = "multiOmicsReactome.RData"
  save(multiOmicsReactome, file = file)
  file.symlink(file, link)
}


moduleSummary <- multiPathwayModuleReport(multiOmicsReactome)
moduleSummary


useThisPathways <- unique(moduleSummary$pathway[moduleSummary$pvalue <= 0.05])
sModule <- moduleSummary[moduleSummary$pathway %in% useThisPathways, , drop = T]


##############################################
####
#### REALLY SLOW
####
#############################################

# if (file.exists("perms.RData")) {
#   load("perms.RData")
# } else {
#   perms <- resampling(fullMultiOmics = multiOmics, reactSmall, nperm = 100, 
#                       pathwaySubset = useThisPathways)
#   file = paste0(dirname, "/perms-", as.character(Sys.Date()), ".RData")
#   link = "perms.RData"
#   save(perms, file = file)
#   file.symlink(file, link)
# }


plotModuleReport(multiOmicsReactome[["PI3K Cascade"]], 
                 fontsize_row = 14,
                 MOcolors = c(exp = "red", met = "green", cnv = "yellow", mut = "blue"))
