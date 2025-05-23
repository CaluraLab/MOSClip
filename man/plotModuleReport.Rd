% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingMO.R
\name{plotModuleReport}
\alias{plotModuleReport}
\title{Plot a table of a \code{MultiOmicsModules} (MOM) object}
\usage{
plotModuleReport(
  modulesObj,
  MOcolors = NULL,
  priority_to = NULL,
  fontsize = 12,
  ...
)
}
\arguments{
\item{modulesObj}{\code{MultiOmicsModules} class object}

\item{MOcolors}{character vector with the omic colors.
The colors should be among the colors in \code{\link{showMOSpalette}}}

\item{priority_to}{a vector with the covariates (omic names) that should go
first}

\item{fontsize}{Size of the font to be used in the plot}

\item{...}{additional argument to be passed to pheatmap}
}
\value{
a Heatmap list object from ComplexHeatmap package of the results
contained in the \code{MultiOmicsModules} object provided
}
\description{
Given a \code{MultiOmicsModules} object, it plots its results in a
tabular fashion
}
\examples{
data(multiOmics)
data(reactSmall)

genesToUse <- row.names(multiOmics[[1]])

MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
    survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
    useTheseGenes = genesToUse
)

plotModuleReport(MOM_survival,
    MOcolors = c(
        exp = "red", met = "green", mut = "blue",
        cnv = "yellow"
    )
)

}
