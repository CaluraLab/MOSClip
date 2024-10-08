% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiOmicsMethod.R
\name{summarizeInCluster}
\alias{summarizeInCluster}
\title{Summarize Using Cluster Analysis}
\usage{
summarizeInCluster(
  data,
  features,
  name = "clu",
  dictionary = NULL,
  max_cluster_number = 3,
  cliques = NULL
)
}
\arguments{
\item{data}{a data matrix}

\item{features}{a vector with the features to analyze}

\item{name}{prefix of the covariates}

\item{dictionary}{translate features (genes) into sets
(row.names of the data)}

\item{max_cluster_number}{the maximum number of cluster to evaluate}

\item{cliques}{the features organized in cliques. Only use for topology}
}
\value{
a list with summary of the omic:
\item{x}{summary of the omic for each sample}
\item{usedGenes}{genes list of genes used to calculate the summary}
\item{namesCov}{names of the covariates}
\item{cls}{the genes in clusters}
\item{method}{method used for the analysis}
\item{omicName}{name of the omic}
}
\description{
Given a matrix it summarize in classes
}
\details{
The user can define a maximum of classes. The function
guess the optimal number of clusters using NbClust methods.
}
