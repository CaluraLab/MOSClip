## Colors
MOSpalette <- c(
    "#E33402", "#04B880", "#3871CC",
    "#E5C002", "#8C04B6", "#01B7BF"
)
names(MOSpalette) <- c("red", "green", "blue", "yellow", "violet", "teal")

MOSpaletteSchema <- data.frame(
    dark = grDevices::adjustcolor(MOSpalette,
        offset = c(
            -0.4, -0.4, -0.4, 0
        )
    ),
    smart = MOSpalette,
    light = grDevices::adjustcolor(MOSpalette,
        offset = c(0.5, 0.5, 0.5, 0)
    ),
    transparent = grDevices::adjustcolor(MOSpalette, alpha.f = 0.6),
    white = grDevices::adjustcolor(MOSpalette,
        offset = c(0.96, 0.96, 0.96, 0)
    ),
    stringsAsFactors = FALSE, check.names = FALSE
)
rownames(MOSpaletteSchema) <- c(
    "red", "green", "blue",
    "yellow", "violet", "teal"
)

# The followings are functions, use it like redShades(100) where 100 is the
# number of shades. In the case of binary map we can use redShades(2).
redShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["red", "white"],
        MOSpaletteSchema["red", "smart"]
    ),
    bias = 1, space = "rgb", interpolate = "linear", alpha = FALSE
)

greenShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["green", "white"],
        MOSpaletteSchema["green", "smart"]
    ), bias = 1, space = "rgb", interpolate = "linear", alpha = FALSE
)

blueShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["blue", "white"],
        MOSpaletteSchema["blue", "smart"]
    ), bias = 1, space = "rgb", interpolate = "linear", alpha = FALSE
)

yellowShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["yellow", "white"],
        MOSpaletteSchema["yellow", "smart"]
    ),
    bias = 1, space = "rgb",
    interpolate = "linear", alpha = FALSE
)

violetShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["violet", "white"],
        MOSpaletteSchema["violet", "smart"]
    ),
    bias = 1, space = "rgb", interpolate = "linear", alpha = FALSE
)

tealShades <- grDevices::colorRampPalette(
    colors = c(
        MOSpaletteSchema["teal", "white"],
        MOSpaletteSchema["teal", "smart"]
    ),
    bias = 1, space = "rgb", interpolate = "linear", alpha = FALSE
)

pvalueShades <- (grDevices::colorRampPalette(
    colors = c("#edf7f5", "#2796bd"), bias = 20,
    space = "rgb", interpolate = "linear", alpha = FALSE
))(100)

omicsRegexp <- "(PC[0-9]+|[1-9]k[1-9]?|TRUE|FALSE|POS|NEG|POSTRUE|NEGTRUE)$"

#' Omics class object with TCGA ovarian data
#'
#' An `Omics` class object containing data from TCGA ovarian cancer. The TCGA
#' data was manually selected and preprocessed. It contains 4 omics: expression,
#' methylation, mutation, and copy number variation. Additionally, it contains
#' specific arguments to perform the dimensionality reduction.
#' The datasets were downloaded from TCGA using `TCGABiolink` R package,
#' selecting only patients with primary solid tumors.
#' Expression matrix was processed first, converting gene identifiers into
#' Entrez IDs. The profiles of genes present more than once were averaged.
#' Genes with at least 100 counts in at least one patients were selected,
#' to avoid data sparsity.
#' Mutation matrix was filtered, keeping only genes with expression data
#' available. We chose to consider only missense and nonsense mutations and
#' mutation impact was also considered following Mutect2 pipeline.
#' CNV values were transformed into numeric values.
#' Methylation data were processed with Methyl Mix R package. Patients that had
#' both normal and primary tumors samples were selected. With the help of a
#' dictionary array probes were connected to CpG clusters, and finally CpG
#' clusters were mapped to genes (Entrez ID).
#' Survival annotation curated by Liu et al. (2018) was used to extract PFS
#' information.
#' Only patients with matched data across the four omics were considered.
#' After the selection of patients and genes, we performed expression
#' normalization and log2 of the counts+1 transformation. This will ensure us
#' to work with expression data approximately close to a normal distribution,
#' the most suitable distribution for the subsequent `MOSClip` tests.
#' Genes and samples were manually selected to create this small example
#' dataset for demonstration purposes.
#'
#' @format ## `multiOmics`
#' An Omics with 4 omics:
#' \describe{
#'   \item{exp}{Matrix with 151 rows and 50 columns of RNA expression values}
#'   \item{met}{A matrix with 178 rows and 50 columns of methylation data with
#'   probes clustered}
#'   \item{mut}{A matrix with 107 rows and 50 columns of mutation counts}
#'   \item{cnv}{A matrix with matrix with 145 rows and 50 columns of copy
#'   number}
#'   ...
#' }
#'
#' @usage data('multiOmics')

"multiOmics"

#' Omics class object with TCGA ovarian data for topological analysis
#'
#' An `Omics` class object containing data from TCGA ovarian cancer.
#' The data are the same as in \code{\link{multiOmics}} object.
#' Arguments in `specificArgs` slot have been set to efficiently run a
#' topological pathway analysis, i.e., the topological method is used for PCA
#' and shrink parameter is set to TRUE.
#' This method can't be used for analyses on modules.
#'
#' @format ## `multiOmicsTopo`
#' An Omics with 4 omics:
#' \describe{
#'   \item{exp}{Matrix with 151 rows and 50 columns of RNA expression values}
#'   \item{met}{A matrix with 178 rows and 50 columns of methylation data with
#'   probes clustered}
#'   \item{mut}{A matrix with 107 rows and 50 columns of mutation counts}
#'   \item{cnv}{A matrix with matrix with 145 rows and 50 columns of copy
#'   number}
#'   ...
#' }
#'
#' @usage data('multiOmicsTopo')

"multiOmicsTopo"

#' PathwayList of pathways from Reactome
#'
#' A PathwayList with three pathways necessary for the analysis:
#' 'Activation of Matrix Metalloproteinases', 'FGFR1 mutant receptor
#' activation', and 'VEGFA-VEGFR2 Pathway'. Pathways were downloaded using
#' `graphite` package and the names of the nodes were converted into
#' Entrez IDs.
#'
#' @format ## `reactSmall`
#' A PathwayList with Reactome pathways for hsapiens
#' \describe{
#'   \item{entries}{Three Reactome pathways with their nodes}
#' }
#'
#' @usage data('reactSmall')

"reactSmall"

#' ExperimentList class object with TCGA ovarian data
#'
#' An ExperimentList class object containing data from TCGA ovarian cancer.
#' The TCGA data was manually selected and preprocessed. It contains 4 omics:
#' expression, methylation, mutation, and copy number variation.
#'
#' @format ## `ExperimentList`
#' An ExperimentList with 4 omics:
#' \describe{
#'   \item{exp}{Matrix with 101 rows and 50 columns of RNA expression values}
#'   \item{met}{A matrix with 97 rows and 50 columnsof methylation data with
#'   probes clustered}
#'   \item{mut}{A matrix with 55 rows and 50 columns of mutation counts}
#'   \item{cnv}{A matrix with matrix with 101 rows and 50 columns of copy
#'   number}
#'   ...
#' }
#'
#' @usage data('ovarianDataset')

"ovarianDataset"

#' MOSClip: Multi-Omics Survival Clip
#'
#' `MOSClip` R package implements a statistical approach able to integrate
#' multi-omic data and look for survival associated gene modules. It
#' integrates multiple omics - trascriptomics, methylomics, genomic
#' mutations, and genomic copy number variations - using various data
#' dimensionality reduction strategies and multivariate models. Exploiting
#' graph theory, pathways can be decomposed into their connected
#' components, that we call modules. The analysis can then be performed at
#' the level of entire pathways or pathway modules. `MOSClip` pathway
#' analysis serves two primary purposes: testing the survival association
#' of pathways or modules using the Cox proportional hazard model, and
#' conducting a two-class analysis with a generalized linear model.
#' Additionally, the package offers valuable graphical tools to visualize
#' and interpret the results.
#'
#' To conduct a multi-omic survival analysis on pathways or modules use:
#' - \code{\link{multiOmicsSurvivalPathwayTest}}
#' - \code{\link{multiOmicsSurvivalModuleTest}}
#'
#' To perform a two-class comparison enrichment analysis on pathways or modules
#' use:
#' - \code{\link{multiOmicsTwoClassPathwayTest}}
#' - \code{\link{multiOmicsTwoClassModuleTest}}
#'
#' @references Paolo Martini, Monica Chiogna, Enrica Calura,
#' and Chiara Romualdi. 2019. “MOSClip: Multi-Omic and Survival Pathway
#' Analysis for the Identification of Survival Associated Gene and Modules.”
#' Nucleic Acids Research 47 (14): e80. \url{https://doi.org/10.1093/nar/gkz324}
#'
#' @docType package
#' @name MOSClip

"_PACKAGE"
