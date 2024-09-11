
## Colors
MOSpalette <- c ("#E33402", "#04B880", "#3871CC", "#E5C002", "#8C04B6",
                 "#01B7BF")
names(MOSpalette) <- c("red","green","blue","yellow","violet","teal")

MOSpaletteSchema <- data.frame(
  dark=grDevices::adjustcolor(MOSpalette,offset = c(-0.4, -0.4, -0.4, 0)),
  smart=MOSpalette,
  light=grDevices::adjustcolor(MOSpalette,offset = c(0.5, 0.5, 0.5, 0)),
  transparent=grDevices::adjustcolor(MOSpalette, alpha.f = 0.6),
  white=grDevices::adjustcolor(MOSpalette,offset = c(0.96, 0.96, 0.96, 0)),
  stringsAsFactors = FALSE, check.names = FALSE)
rownames(MOSpaletteSchema) <- c("red","green","blue","yellow","violet","teal")

# The followings are functions, use it like redShades(100) where 100 is the 
# number of shades. In the case of binary map we can use redShades(2).
redShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["red","white"],
           MOSpaletteSchema["red","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

greenShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["green","white"],
           MOSpaletteSchema["green","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

blueShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["blue","white"],
           MOSpaletteSchema["blue","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

yellowShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["yellow","white"],
           MOSpaletteSchema["yellow","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

violetShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["violet","white"],
           MOSpaletteSchema["violet","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

tealShades <- grDevices::colorRampPalette(
  colors=c(MOSpaletteSchema["teal","white"],
           MOSpaletteSchema["teal","smart"]),
  bias=1, space="rgb",
  interpolate="linear", alpha=FALSE)

pvalueShades <- grDevices::colorRampPalette(
  colors=c("#edf7f5", "#2796bd"),
  bias=20, space="rgb",
  interpolate="linear", alpha=FALSE)(100)

omicsRegexp <- "(PC[0-9]+|[1-9]k[1-9]?|TRUE|FALSE|POS|NEG|POSTRUE|NEGTRUE)$"

#' Omics class object with TCGA ovarian data
#'
#' A Omics class object containing data from TCGA ovarian cancer. The TCGA data 
#' was manually selected and preprocessed. It contains 4 omics: expression,
#' methylation, mutation, and copy number variation. Additionally, it contains
#' specific arguments to perform the dimensionality
#' reduction
#'
#' @format ## `multiOmics`
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

"multiOmics"

#' Omics class object with TCGA ovarian data for two-class analysis
#'
#' A Omics class object containing data from TCGA ovarian cancer. The TCGA data 
#' was manually selected and preprocessed. It contains 4 omics: expression,
#' methylation, mutation, and copy number variation. Additionally, it contains
#' specific arguments to perform the dimensionality reduction and colData.
#'
#' @format ## `multiOmicsTwoClass`
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

"multiOmicsTwoClass"

#' PathwayList of pathways from Reactome
#'
#' A PathwayList with three pathways necessary for the analysis:
#' "Intrinsic Pathway for Apoptosis", "Opioid Signalling", and 
#' "Mitochondrial protein import"
#'
#' @format ## `reactSmall`
#' A PathwayList with Reactome pathways for hsapiens
#' \describe{
#'   \item{entries}{Three reactome pathways with their nodes}
#' }

"reactSmall"
