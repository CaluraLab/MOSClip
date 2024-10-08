% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{multiOmics}
\alias{multiOmics}
\title{Omics class object with TCGA ovarian data}
\format{
\subsection{\code{multiOmics}}{

An Omics with 4 omics:
\describe{
\item{exp}{Matrix with 151 rows and 50 columns of RNA expression values}
\item{met}{A matrix with 178 rows and 50 columns of methylation data with
probes clustered}
\item{mut}{A matrix with 107 rows and 50 columns of mutation counts}
\item{cnv}{A matrix with matrix with 145 rows and 50 columns of copy
number}
...
}
}
}
\usage{
data('multiOmics')
}
\description{
An \code{Omics} class object containing data from TCGA ovarian cancer. The TCGA
data was manually selected and preprocessed. It contains 4 omics: expression,
methylation, mutation, and copy number variation. Additionally, it contains
specific arguments to perform the dimensionality reduction.
The datasets were downloaded from TCGA using \code{TCGABiolink} R package,
selecting only patients with primary solid tumors.
Expression matrix was processed first, converting gene identifiers into
Entrez IDs. The profiles of genes present more than once were averaged.
Genes with at least 100 counts in at least one patients were selected,
to avoid data sparsity.
Mutation matrix was filtered, keeping only genes with expression data
available. We chose to consider only missense and nonsense mutations and
mutation impact was also considered following Mutect2 pipeline.
CNV values were transformed into numeric values.
Methylation data were processed with Methyl Mix R package. Patients that had
both normal and primary tumors samples were selected. With the help of a
dictionary array probes were connected to CpG clusters, and finally CpG
clusters were mapped to genes (Entrez ID).
Survival annotation curated by Liu et al. (2018) was used to extract PFS
information.
Only patients with matched data across the four omics were considered.
After the selection of patients and genes, we performed expression
normalization and log2 of the counts+1 transformation. This will ensure us
to work with expression data approximately close to a normal distribution,
the most suitable distribution for the subsequent \code{MOSClip} tests.
Genes and samples were manually selected to create this small example
dataset for demonstration purposes.
}
\keyword{datasets}
