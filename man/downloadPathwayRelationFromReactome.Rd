% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/houseClip-pathway-relations.R
\name{downloadPathwayRelationFromReactome}
\alias{downloadPathwayRelationFromReactome}
\title{Download Reactome Pathway Relations}
\usage{
downloadPathwayRelationFromReactome(url = NULL, speciesAbbr = "HSA")
}
\arguments{
\item{url}{the location of the file. Can be local. If NULL pick the package
reactome file.}

\item{speciesAbbr}{species acronim}
}
\value{
A data frame with 2 columns:
\item{parent}{The Reactome pathway ID of the parent pathway.}
\item{child}{The Reactome pathway ID of the child pathway.}
}
\description{
Download Pathway Relations from Reactome. The file is retrieved from the
\href{https://reactome.org/download/current/ReactomePathwaysRelation.txt}{url}
}
\examples{
downloadPathwayRelationFromReactome()

}
