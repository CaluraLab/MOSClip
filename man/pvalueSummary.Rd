% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supertest&frequencyplot.R
\name{pvalueSummary}
\alias{pvalueSummary}
\title{Compute pvalue Summary}
\usage{
pvalueSummary(multiPathwayReportData, excludeColumns = NULL, as.list = FALSE)
}
\arguments{
\item{multiPathwayReportData}{data.frame, the output of the
\code{\link{multiPathwayReport}} or \code{\link{multiPathwayModuleReport}}
functions.}

\item{excludeColumns}{a vector of characters listing the columns of
\code{multiPathwayReportData} object to be excluded by the analysis.
In the case \code{multiPathwayReportData} derives from
\code{\link{multiPathwayModuleReport}}
you should set \code{excludeColumns = c('pathway','module')}.}

\item{as.list}{return a list rather than a data.frame}
}
\value{
a list
}
\description{
Compute pvalue Summary
}
