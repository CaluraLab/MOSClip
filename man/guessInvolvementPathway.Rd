% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/guessInvolvement.R
\name{guessInvolvementPathway}
\alias{guessInvolvementPathway}
\title{Guess the most influent features from MultiOmics Survival or Two-class
results.}
\usage{
guessInvolvementPathway(
  pathway,
  loadThr = 0.6,
  n = 3,
  atleast = 1,
  min_prop_pca = 0.1,
  min_prop_events = 0.1,
  ...
)
}
\arguments{
\item{pathway}{\code{MultiOmicsModules} object from a pathway}

\item{loadThr}{the loading threshold to select genes (PCA only)}

\item{n}{the maximum number of genes to retrieve (cluster and binary only)}

\item{atleast}{the minimum number of features to select (PCA only)}

\item{min_prop_pca}{the minimal proportion to compute the PCA classes}

\item{min_prop_events}{the minimal proportion to compute the event classes}

\item{...}{additional arguments passed to \code{get} function}
}
\value{
a list. Each item of the list corresponds to an omic that is
summarized with the specific 'extractSummary' functions. Each item is the
summary for an omic summarized using the setted method: pvalues are present
only for cluster method.
}
\description{
Given a pathway analyzed by \code{multiOmicsSurvivalPathwayTest} or
\code{multiOmicsTwoClassPathwayTest}, it retrieves for each omic the most
influent features.
}
