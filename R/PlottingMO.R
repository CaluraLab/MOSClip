#' Plot heatmaps of the pathway by omics
#'
#' Given the pathway, it creates the heatmaps of the mostly involved genes for
#' each omic.
#'
#' @param pathway `MultiOmicsPathway` class object
#' @param sortBy one or more covariates to sort the samples
#' @param paletteNames name of the colors for each omic
#' @param additionalAnnotations optional additional sample annotations
#' (e.g. survival annotation)
#' @param additionalPaletteNames colors for additional annotations. The colors
#' available are the ones in \code{\link{showMOSpalette}}
#' @param discr_prop_pca the minimal proportion to compute the PCA classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param withSampleNames show the sample names in the plot
#' @param nrowsHeatmaps magnification respect to annotation of sample
#' (annotations take 1 row)
#' @param orgDbi a Dbi organism to be used. Default is `org.Hs.eg.db`
#' @param ... additional arguments passed to `guessInvolvementPathway` function
#' (internal use)
#'
#' @return An object of class `ggplot` plotted with ComplexHeatMap package.
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' survAnnot <- data.frame(
#'     status = multiOmics$status,
#'     days = multiOmics$days,
#'     row.names = colnames(multiOmics[[1]])
#' )
#'
#' # Creating the MultiOmicsPathway object
#' MOP_survival <- multiOmicsSurvivalPathwayTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' # Plotting
#' plotPathwayHeat(MOP_survival,
#'     sortBy = c("expPC2", "mut", "status", "days"),
#'     paletteNames = c(exp = "red", met = "green", 
#'                      mut = "blue", cnv = "yellow"),
#'     additionalAnnotations = survAnnot,
#'     additionalPaletteNames = list(status = "teal", days = "violet"),
#'     nrowsHeatmaps = 2, withSampleNames = F
#' )
#'
#' @importFrom checkmate assertClass
#' @importFrom grid gpar grid.newpage grid.draw rectGrob
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation '%v%' draw
#' @importFrom stats relevel
#' @importFrom ggplotify as.ggplot
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @export

plotPathwayHeat <- function(
    pathway, sortBy = NULL, paletteNames = NULL, additionalAnnotations = NULL,
    additionalPaletteNames = NULL, discr_prop_pca = 0.15,
    discr_prop_events = 0.05, withSampleNames = TRUE, nrowsHeatmaps = 3,
    orgDbi = "org.Hs.eg.db", ...) {
    checkmate::assertClass(pathway, "MultiOmicsPathway")

    involved <- guessInvolvementPathway(pathway,
        min_prop_pca = discr_prop_pca,
        min_prop_events = discr_prop_events,
        ...
    )

    if (length(paletteNames) != length(involved)) {
        repTimes <- ceiling(length(involved) / length(paletteNames))
        paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
    }

    annotationFull <- formatAnnotations(involved, sortBy = NULL)
    idx <- which(unlist(lapply(annotationFull, class)) == "numeric")
    if (length(idx) > 0) {
        annotationFull[, idx] <- lapply(annotationFull[, idx], function(col) {
            stats::relevel(as.factor(col), ref = "1")
        })
    }

    omics <- guessOmics(colnames(annotationFull))
    if (is.null(paletteNames)) {
        paletteNames <- names(paletteNames)[seq_along(unique(omics))]
        paletteNames <- guessOmicsColors(omics)
    }

    if (length(paletteNames) != length(unique(omics))) {
        stop(paste0(
            "Length of MOcolors differs from the number of omics: ",
            paste(unique(omics), collapse = " ,")
        ))
    }

    if (is.null(names(paletteNames))) {
        names(paletteNames) <- unique(omics)
    }

    annotationPalettes <- paletteNames

    ann_col <- lapply(colnames(annotationFull), function(name) {
        omic <- guessOmic(name)
        if (!omic %in% names(annotationPalettes)) {
            stop("Missing omic in annotationPalettes")
        }
        discreteColor <- annotationPalettes[[omic]]
        values <- sort(unique(annotationFull[, name]))
        if (!is.null(levels(values))) {
            values <- levels(values)
        }
        if (length(values) == 1) {
            annot <- as.character(MOSpaletteSchema[discreteColor, c("smart")])
            names(annot) <- values
        } else if (length(values) == 2) {
            annot <- as.character(MOSpaletteSchema[
                discreteColor,
                c("smart", "light")
            ])
            names(annot) <- values
        } else if (length(table(annotationFull[, name])) == 3) {
            annot <- as.character(MOSpaletteSchema[
                discreteColor,
                c("dark", "smart", "light")
            ])
            names(annot) <- values
        } else {
            stop("I'm puzzled! Too many values to map")
        }
        annot
    })
    names(ann_col) <- colnames(annotationFull)

    if (!is.null(additionalAnnotations)) {
        additionalAnnotations <- matchAnnotations(
            annotationFull,
            additionalAnnotations
        )
        annotationFull <- cbind(annotationFull, additionalAnnotations)

        if (!is.null(additionalPaletteNames)) {
            add_ann_col <- lapply(colnames(additionalAnnotations), 
                                  function(name) {
                values <- sort(unique(additionalAnnotations[[name]]))
                discreteColor <- additionalPaletteNames[[name]]

                if (length(values) == 2) {
                    annot <- as.character(MOSpaletteSchema[
                        discreteColor,
                        c("smart", "light")
                    ])
                    names(annot) <- values
                } else if (length(table(annotationFull[, name])) == 3) {
                    annot <- as.character(MOSpaletteSchema[
                        discreteColor,
                        c("dark", "smart", "light")
                    ])
                    names(annot) <- values
                } else {
                    annot <- colorRamp2(
                        c(min(values), max(values)),
                        getContinousPalette(discreteColor, 2)
                    )
                    names(annot) <- levels(values)
                }
                annot
            })
            names(add_ann_col) <- colnames(additionalAnnotations)
            ann_col <- c(ann_col, add_ann_col)
        }
    }

    annotationFull <- sortAnnotations(annotationFull, sortBy)
    annotationFull <- annotationFull[, rev(colnames(annotationFull))]

    ha <- HeatmapAnnotation(df = annotationFull, col = ann_col)
    ht_list <- ha
    for (n in seq_along(involved)) {
        heatMatrix <- involved[[n]]$sigModule
        heatMatrix <- heatMatrix[, row.names(annotationFull), drop = FALSE]
        row.names(heatMatrix) <- conversionToSymbols(row.names(heatMatrix), 
                                                     orgDbi)
        Ht <- Heatmap(heatMatrix,
            cluster_columns = FALSE, cluster_rows = FALSE,
            show_column_names = FALSE,
            heatmap_legend_param = list(title = NULL),
            col = c("#FFFFFF",
                    rev(as.vector(ann_col[involved[[n]]$covsConsidered][[1]]))
                    )
            )
        ht_list <- ht_list %v% Ht
    }

    gb <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,
        legend_grouping = "original"
    ))
    gb <- ggplotify::as.ggplot(gb)
    return(gb)
}

#' Plot Kaplan-Meier survival curves of a specific pathway
#'
#' Given a `MultiOmicsPathway` class object, it plots Kaplan-Meier curves,
#' in which the strata corresponds to the chosen omics
#'
#' @param pathway `MultiOmicsPathway` class object
#' @param formula a formula to compute the plot
#' @param fileName optional filenames to save the plot
#' @param paletteNames a palette containing three colors
#' @param h the height of the plot
#' @param w the width of the plot
#' @param risk.table logical value. If TRUE, shows the `risk.table`. Default
#' is TRUE.
#' @param pval logical value. Shows p-value of the curves
#' @param size line width of the KM curves
#' @param inYears logical value. If TRUE, converts days to years
#' @param discr_prop_pca the minimal proportion to compute the PCA classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to
#' include
#' @param additional_continuous names of the additional continuous variables to
#' include
#' @param ... additional arguments passed to `guessInvolvementPathway` and
#' `get` function (internal use)
#'
#' @return a ggsurvplot class object
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' # Creating the MultiOmicsPathway object
#' MOP_survival <- multiOmicsSurvivalPathwayTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' plotPathwayKM(MOP_survival,
#'     formula = "Surv(days, status) ~ mut + expPC2",
#'     paletteNames = "Paired", inYears = TRUE
#' )
#'
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#'
#' @export

plotPathwayKM <- function(
    pathway, formula = "Surv(days, status) ~ PC1", fileName = NULL,
    paletteNames = NULL, h = 9, w = 7, risk.table = TRUE, pval = TRUE, size = 1,
    inYears = FALSE, discr_prop_pca = 0.15, discr_prop_events = 0.05,
    additional_discrete = NULL, additional_continuous = NULL, ...) {
    checkmate::assertClass(pathway, "MultiOmicsPathway")

    involved <- guessInvolvementPathway(pathway,
        min_prop_pca = discr_prop_pca, min_prop_events = discr_prop_events, ...)

    annotationFull <- formatAnnotations(involved, sortBy = NULL)
    multiOmicObj <- get(pathway@multiOmicObj, ...)

    coxObj <- createCoxObj(multiOmicObj@colData, pathway@pathView)
    daysAndStatus <- coxObj[, c("status", "days"), drop = FALSE]

    if (inYears) {
        daysAndStatus$days <- daysAndStatus$days / 365.24
    }

    coxObj <- data.frame(
        daysAndStatus,
        annotationFull[row.names(daysAndStatus), , drop = FALSE]
    )

    fit <- survminer::surv_fit(formula(formula), data = coxObj)

    palette <- NULL
    if (!is.null(paletteNames)) {
        if (length(paletteNames) == 1) {
            palette <- paletteNames
        } else {
            classes <- names(fit$strata)
            if (length(classes) == length(paletteNames)) {
                palette <- paletteNames
            }
        }
    }

    p <- survminer::ggsurvplot(fit,
        data = coxObj, risk.table = risk.table, pval = pval,
        palette = palette, size = size
    )

    if (!is.null(fileName)) {
        ggplot2::ggsave(filename = fileName, p, height = h, width = w)
    } else {
        p
    }
}

#' Plot a Heatmap of a Module by Omics
#'
#' It creates a heatmap of the most involved genes of each omic of a specific
#' module from a `MultiOmicsModule` object.
#'
#' @param moduleobj `MultiOmicsModule` class object
#' @param moduleNumber module number of interest
#' @param sortBy a covariate (omic) to sort by
#' @param paletteNames a palette containing three colors
#' @param additionalAnnotations optional additional sample annotations
#' @param additionalPaletteNames optional additional colors for annotations
#' @param withSampleNames show sample names
#' @param fontsize_row font size for row labels
#' @param fontsize_col font size for column labels
#' @param nrowsHeatmaps magnification respect to annotation of sample
#' (annotations take 1 row)
#' @param orgDbi a Dbi organism to be used. Default is `org.Hs.eg.db`
#' @param discr_prop_pca the minimal proportion to compute the PCA classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param ... additional arguments passed to `guessInvolvement` function
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' survAnnot <- data.frame(
#'     status = multiOmics$status,
#'     days = multiOmics$days,
#'     row.names = colnames(multiOmics[[1]])
#' )
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' plotModuleHeat(MOM_survival, 1,
#'     sortBy = c("mut", "expPC1", "status", "days"),
#'     additionalAnnotations = survAnnot,
#'     additionalPaletteNames = list(status = "teal", days = "violet"),
#'     withSampleNames = F
#' )
#'
#' @return A heatmap of a pathway module (results of the module test)
#'
#' @importFrom checkmate assertClass
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation '%v%'
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid gpar grid.newpage grid.draw rectGrob grid.grabExpr
#' @importFrom graphics plot
#' @importFrom stats relevel
#' @importFrom ggplotify as.ggplot
#'
#' @export

plotModuleHeat <- function(
    moduleobj, moduleNumber, sortBy = NULL, paletteNames = NULL,
    additionalAnnotations = NULL, additionalPaletteNames = NULL,
    withSampleNames = TRUE, fontsize_row = 10, fontsize_col = 1,
    nrowsHeatmaps = 3, orgDbi = "org.Hs.eg.db",
    discr_prop_pca = 0.15, discr_prop_events = 0.05, ...) {
    checkmate::assertClass(moduleobj, "MultiOmicsModules")

    moduleGenes <- moduleobj@modules[[moduleNumber]]

    involved <- guessInvolvement(moduleobj,
        moduleNumber = moduleNumber, min_prop_pca = discr_prop_pca,
        min_prop_events = discr_prop_events, ...)

    if (length(paletteNames) != length(involved)) {
        repTimes <- ceiling(length(involved) / length(paletteNames))
        paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
    }

    # Create annotation and sort
    annotationFull <- formatAnnotations(involved, sortBy = NULL)
    idx <- which(unlist(lapply(annotationFull, class)) == "numeric")
    if (length(idx) > 0) {
        annotationFull[, idx] <- lapply(annotationFull[, idx], function(col) {
            stats::relevel(as.factor(col), ref = "1")
        })
    }

    omics <- guessOmics(colnames(annotationFull))
    if (is.null(paletteNames)) {
        paletteNames <- names(paletteNames)[seq_along(unique(omics))]
        paletteNames <- guessOmicsColors(omics)
    }

    if (length(paletteNames) != length(unique(omics))) {
        stop(paste0(
            "Length of MOcolors differs from the number of omics: ",
            paste(unique(omics),
                  collapse = ", ")
        ))
    }

    if (is.null(names(paletteNames))) {
        names(paletteNames) <- unique(omics)
    }

    annotationPalettes <- paletteNames

    ann_col <- lapply(colnames(annotationFull), function(name) {
        omic <- guessOmic(name)
        if (!omic %in% names(annotationPalettes)) {
            stop("Missing omic in annotationPalettes")
        }
        discreteColor <- annotationPalettes[[omic]]
        values <- sort(unique(annotationFull[, name]))
        if (length(values) == 1) {
            annot <- as.character(MOSpaletteSchema[discreteColor, c("smart")])
            names(annot) <- values
        } else if (length(values) == 2) {
            annot <- as.character(MOSpaletteSchema[
                discreteColor,
                c("smart", "light")
            ])
            names(annot) <- values
        } else if (length(table(annotationFull[, name])) == 3) {
            annot <- as.character(MOSpaletteSchema[
                discreteColor,
                c("dark", "smart", "light")
            ])
            names(annot) <- values
        } else {
            stop("I'm puzzled! Too many values to map")
        }
        annot
    })
    names(ann_col) <- colnames(annotationFull)

    if (!is.null(additionalAnnotations)) {
        additionalAnnotations <- matchAnnotations(
            annotationFull,
            additionalAnnotations
        )
        annotationFull <- cbind(annotationFull, additionalAnnotations)

        if (!is.null(additionalPaletteNames)) {
            add_ann_col <- lapply(colnames(additionalAnnotations), 
                                  function(name) {
                values <- sort(unique(additionalAnnotations[[name]]))
                discreteColor <- additionalPaletteNames[[name]]

                if (length(values) == 2) {
                    annot <- as.character(MOSpaletteSchema[
                        discreteColor,
                        c("smart", "light")
                    ])
                    names(annot) <- values
                } else if (length(table(annotationFull[, name])) == 3) {
                    annot <- as.character(MOSpaletteSchema[
                        discreteColor,
                        c("dark", "smart", "light")
                    ])
                    names(annot) <- values
                } else {
                    annot <- colorRamp2(
                        c(min(values), max(values)),
                        getContinousPalette(discreteColor, 2)
                    )
                    names(annot) <- levels(values)
                }
                annot
            })
            names(add_ann_col) <- colnames(additionalAnnotations)
            ann_col <- c(ann_col, add_ann_col)
        }
    }

    annotationFull <- sortAnnotations(annotationFull, sortBy)
    annotationFull <- annotationFull[, rev(colnames(annotationFull))]

    ha <- HeatmapAnnotation(df = annotationFull, col = ann_col)
    ht_list <- ha
    for (n in seq_along(involved)) {
        heatMatrix <- involved[[n]]$sigModule
        heatMatrix <- heatMatrix[, row.names(annotationFull), drop = FALSE]
        row.names(heatMatrix) <- conversionToSymbols(row.names(heatMatrix), 
                                                     orgDbi)

        Ht <- Heatmap(heatMatrix,
            cluster_columns = FALSE, cluster_rows = FALSE,
            show_column_names = FALSE, 
            heatmap_legend_param = list(title = NULL),
            col = c("#FFFFFF", 
                    as.vector(ann_col[[involved[[n]]$covsConsidered[1]]]))
        )
        ht_list <- ht_list %v% Ht
    }

    gb <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_list,
        legend_grouping = "original"
    ))
    
    ggplot_obj <- ggplotify::as.ggplot(gb)
    return(ggplot_obj)
}


#' Plot Kaplan-Meier survival curves of a specific module
#'
#' Given a `MultiOmicsModule` class object and a specific module number, it
#' plots Kaplan-Meier curves, in which the strata corresponds to the omics
#'
#' @param MOM a `MultiOmicsModule` class object
#' @param moduleNumber numeric value. The module number of interest
#' @param formula a formula for the survival analysis. It should be written as
#' 'Surv(days, status) ~ omic'. To plot more than one omic, write them
#' separated by a '+' character after the separator (~)
#' @param fileName optional filenames to save the plot
#' @param paletteNames a palette name to be used
#' @param h the height of the plot
#' @param w the width of the plot
#' @param risk.table logical value. If TRUE, shows the `risk.table`. Default
#' is TRUE.
#' @param pval logical value. If TRUE, shows the p-value of the curves. Default
#' is TRUE.
#' @param size line width of the KM curves
#' @param inYears set time in years
#' @param discr_prop_pca the minimal proportion to compute the PCA classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to
#' include
#' @param additional_continuous names of the additional continous variables to
#' include
#' @param ... additional arguments passed to `guessInvolvement` and `get`
#' function
#'
#' @return a ggsurvplot class object
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' plotModuleKM(MOM_survival, 1,
#'     formula = "Surv(days, status) ~ mut + expPC2",
#'     paletteNames = "Paired", inYears = TRUE
#' )
#'
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#'
#' @export

plotModuleKM <- function(
    MOM, moduleNumber, formula = "Surv(days, status) ~ PC1",
    fileName = NULL, paletteNames = NULL, h = 9, w = 7, risk.table = TRUE,
    pval = TRUE, size = 1, inYears = FALSE, discr_prop_pca = 0.15,
    discr_prop_events = 0.05, additional_discrete = NULL,
    additional_continuous = NULL, ...) {
    checkmate::assertClass(MOM, "MultiOmicsModules")

    involved <- guessInvolvement(MOM,
        moduleNumber = moduleNumber, min_prop_pca = discr_prop_pca,
        min_prop_events = discr_prop_events, ...)

    multiOmicObj <- get(MOM@multiOmicObj, ...)
    coxObj <- createCoxObj(multiOmicObj@colData, 
                           MOM@modulesView[[moduleNumber]])

    annotationFull <- formatAnnotations(involved, sortBy = NULL)

    days_status_names <- c("status", "days")

    daysAndStatus <- coxObj[, days_status_names, drop = FALSE]
    if (inYears) {
        daysAndStatus$days <- daysAndStatus$days / 365.24
    }

    additional_clinic <- NULL
    if (!is.null(additional_discrete)) {
        not_found <- setdiff(additional_discrete, colnames(coxObj))
        if (length(not_found) > 0) {
            stop(paste0("Some discrete variables were not found: ", 
                        paste(not_found, collapse = ", ", "\n"), 
                        "We found the following variables: ", 
                        paste(colnames(coxObj), collapse = ", ")))
        }

        fixed_covs <- setdiff(additional_discrete, colnames(annotationFull))
        additional_clinic <- coxObj[row.names(daysAndStatus), 
                                    fixed_covs, drop = FALSE]
    }

    if (!is.null(additional_continuous)) {
        not_found <- setdiff(additional_continuous, colnames(coxObj))
        if (length(not_found) > 0) {
            stop(paste0("Some continuous variables were not found: ", 
                        paste(not_found, collapse = ", ", "\n"), 
                        "We found the following variables: ", 
                        paste(colnames(coxObj), collapse = ", ")))
        }

        fixed_covs <- setdiff(additional_continuous, colnames(annotationFull))
        df <- cbind(daysAndStatus, coxObj[row.names(daysAndStatus), 
                                          fixed_covs, drop = FALSE])
        discretized_covs <- createDiscreteClasses(df, fixed_covs)

        if (!is.null(additional_clinic)) {
            additional_clinic <- cbind(additional_clinic, 
                                       discretized_covs[, fixed_covs,
                                                        drop = FALSE])
        } else {
            additional_clinic <- discretized_covs[, fixed_covs, drop = FALSE]
        }
    }

    if (is.null(additional_clinic)) {
        coxObj <- data.frame(daysAndStatus, 
                             annotationFull[row.names(daysAndStatus), ,
                                            drop = FALSE])
    } else {
        coxObj <- data.frame(daysAndStatus, 
                             annotationFull[row.names(daysAndStatus), ,
                                            drop = FALSE], 
                             additional_clinic)
    }

    fit <- survminer::surv_fit(formula(formula), data = coxObj)

    palette <- NULL
    if (!is.null(paletteNames)) {
        if (length(paletteNames) == 1) {
            palette <- paletteNames
        } else {
            classes <- names(fit$strata)
            if (length(classes) == length(paletteNames)) {
                palette <- paletteNames
            }
            if (!is.null(names(paletteNames))) {
                diff <- setdiff(names(paletteNames), classes)
                if (length(diff) != 0) {
                    stop("Names of paletteNames must be equal to classes")
                }
            }
        }
    }

    p <- survminer::ggsurvplot(fit,
        data = coxObj, risk.table = risk.table, pval = pval,
        palette = unname(palette), size = size
    )

    if (!is.null(fileName)) {
        ggplot2::ggsave(filename = fileName, p, height = h, width = w)
    } else {
        p
    }
}

#' Plot a Directed Graph of the MultiOmicsModules Object
#'
#' From a `MultiOmicsModules` object, it plots the position of a given module
#' in the pathway. The omics are also represented in the graph.
#'
#' @param modulesobj a `MultiOmicsModule` class object
#' @param pathList a `PathwayList` from `graphite` package that contains the
#' pathways to be used
#' @param moduleNumber a module number
#' @param orgDbi if needed, indicates an organism Dbi to translate the vectors
#' @param legendLabels set up your favourite names for the omics
#' @param paletteNames named vector of MOSpalettes, names replace makeLegend
#' arguments
#' @param fileName optional filenames to save the plot
#' @param discr_prop_pca the minimal proportion to compute the PCA classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param pathTitle title of the graph, to be searched in `pathList`
#' @param ... additional arguments passed to `guessInvolvement` function
#'
#' @return a MOSClip plot in form of a list class object
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' plotModuleInGraph(MOM_survival, reactSmall,
#'     moduleNumber = 1,
#'     paletteNames = c(exp = "red", met = "green", 
#'                      mut = "blue", cnv = "yellow")
#' )
#'
#' @importFrom checkmate assertClass
#' @importFrom igraph V V<- simplify graph_from_graphnel
#' @importFrom AnnotationDbi select
#' @importFrom graphics plot legend par
#' @importFrom grDevices dev.off pdf rainbow
#'
#' @export
plotModuleInGraph <- function(
    modulesobj, pathList, moduleNumber, orgDbi = "org.Hs.eg.db",
    paletteNames = NULL, legendLabels = NULL, fileName = NULL, 
    discr_prop_pca = 0.15,
    discr_prop_events = 0.05, pathTitle = NULL, ...) {
    checkmate::assertClass(modulesobj, "MultiOmicsModules")

    if (!is(pathList[[1]], "graphNEL")) {
        net <- igraph::graph_from_graphnel(convertPathway(
            pathList[[modulesobj@title]], NULL))
    } else {
        net <- igraph::graph_from_graphnel(pathList[[pathTitle]])
    }

    moduleGenes <- modulesobj@modules[[moduleNumber]]
    net <- igraph::simplify(net, remove.multiple = TRUE, remove.loops = TRUE)
    color <- rep("grey", length(V(net)))
    color[names(V(net)) %in% moduleGenes] <- "tomato"

    involved <- guessInvolvement(modulesobj,
        moduleNumber = moduleNumber, min_prop_pca = discr_prop_events,
        min_prop_events = discr_prop_events, ...)
    mark.groups <- lapply(involved, function(x) {
        row.names(x$subset)
    })

    group.names <- vapply(involved, function(x) {
        guessOmic(x$covsConsidered)}, character(1))

    colLength <- length(mark.groups)
    if (colLength < 3) {
        mark.col <- rainbow(3, alpha = 0.33)[seq_len(colLength)]
    } else {
        mark.col <- rainbow(colLength, alpha = 0.33)
    }

    mark.border <- NA

    if (!is.null(paletteNames)) {
        if (is.null(names(paletteNames))) {
            stop("paletteNames must be named vector")
        }

        if (!is.null(names(paletteNames))) {
            mismatch <- setdiff(group.names, names(paletteNames))
            if (length(mismatch) > 0) {
                stop("Missing palette for one or more omics")
            }
            paletteNames <- paletteNames[group.names]
        }

        err <- setdiff(paletteNames, row.names(MOSpaletteSchema))
        if (length(err) != 0) {
            stop("paletteNames value not allowed.")
        }

        mark.col <- MOSpaletteSchema[paletteNames, ]$transparent
    }

    if (is.null(legendLabels)) {
        legendLabels <- group.names
    } else {
        if (length(legendLabels) != length(group.names)) {
            warning("Some legendLabels were not found in data")
        }
        legendLabels <- legendLabels[seq_along(group.names)]
    }

    labels <- conversionToSymbols(names(V(net)), orgDbi)

    par(mar = c(1.5, 1, 1, 0.5))

    if (!is.null(fileName)) {
        pdf(fileName)
    }
    plot(net,
        edge.arrow.size = 0.5, edge.curved = 0.2, vertex.label = labels, 
        vertex.label.cex = 0.6, vertex.label.family = "sans", 
        vertex.label.font = 2, vertex.color = color, 
        vertex.frame.color = "gray", vertex.label.color = "black", 
        vertex.size = 15, mark.groups = mark.groups, mark.col = mark.col, 
        mark.border = NA
    )
    legend(
        x = -1, y = -1, legendLabels, pch = 21, horiz = TRUE, col = "#777777",
        pt.bg = mark.col, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1
    )
}

#' Summarize and plot pathways' info from a list of `MultiOmicsPathway` (MOP)
#'
#' Given the list of MOPs, it plots a table of its results.
#'
#' @param multiPathwayList a `list` of `MultiOmicsPathway` class objects
#' @param top numeric value. Plot only the top number of pathways
#' @param MOcolors character vector with the omic colors.
#' The colors should be among the colors in \code{\link{showMOSpalette}}
#' @param priority_to a vector with the covariates (omic names) that should go
#' first
#' @param fontsize the font size to be used. Default is 12.
#' @param \dots additional argument to be passed to pheatmap
#'
#' @return a Heatmap list object from ComplexHeatmap package of the results
#' contained in the `MultiOmicsPathway` object provided
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' MOP_list <- lapply(reactSmall, function(g) {
#'     fcl <- multiOmicsSurvivalPathwayTest(multiOmics, g,
#'         survFormula = "Surv(days, status) ~",
#'         autoCompleteFormula = TRUE,
#'         useTheseGenes = genesToUse
#'     )
#'     fcl
#' })
#'
#' plotMultiPathwayReport(MOP_list,
#'     MOcolors = c(
#'         exp = "red", met = "green", mut = "blue",
#'         cnv = "yellow"
#'     ),
#'     fontsize = 12
#' )
#'
#' @importFrom pheatmap pheatmap
#'
#' @export

plotMultiPathwayReport <- function(
    multiPathwayList, top = 25, MOcolors = NULL, priority_to = NULL,
    fontsize = 6, ...) {
    if (!is.list(multiPathwayList)) {
        stop("multiPathwayList must be a list.")
    }

    summary <- multiPathwayReport(multiPathwayList)
    top <- min(top, NROW(summary))

    annCol <- guessOmics(colnames(summary))
    omics <- annCol[2:length(annCol)]

    if (is.null(priority_to) & (!is.null(MOcolors))) {
        if (!is.null(names(MOcolors))) {
            priority_to <- names(MOcolors)
        }
    }

    if (is.null(MOcolors)) {
        MOcolors <- guessOmicsColors(omics)
    }

    if (length(MOcolors) != length(unique(omics))) {
        stop(paste0("Length of MOcolors differs from the number of omics: ", 
                    paste(unique(omics), collapse = ", ")))
    }

    if (is.null(names(MOcolors))) {
        names(MOcolors) <- unique(omics)
    }

    omics <- omics[order(match(omics, priority_to))]
    colors <- createColors(omics, MOcolors)
    names(colors) <- unique(omics)

    msummary <- as.matrix(summary[seq_len(top), 2:(ncol(summary))])

    msummary <- order_by_covariates(msummary, 0, priority_to)

    cell_text <- function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", msummary[i, j]), x, y, 
                  gp = gpar(fontsize = fontsize))
    }

    ta <- HeatmapAnnotation(Omics = omics, col = list(Omics = colors))

    ht1 <- Heatmap(
        matrix = summary$pvalue[seq_len(top)], 
        col = pvalueShades[seq_len(round((summary$pvalue[top] *100) + 1))], 
        cluster_rows = FALSE, show_heatmap_legend = FALSE, name = "Pvalue",
        column_names_gp = gpar(fontsize = fontsize), 
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", summary$pvalue[i]), x, y, 
                      gp = gpar(fontsize = fontsize))}
    )
    dots <- list(...)
    defargs <- list(
        matrix = msummary, name = "p-Value", col = pvalueShades, 
        cluster_rows = FALSE, cluster_columns = FALSE, cell_fun = cell_text, 
        top_annotation = ta, row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = fontsize)
    )
    args <- matchArguments(dots, defargs)

    ht2 <- do.call(Heatmap, args)
    ht1 + ht2
}

#' Plot a table of a `MultiOmicsModules` (MOM) object
#'
#' Given a `MultiOmicsModules` object, it plots its results in a
#' tabular fashion
#'
#' @inheritParams plotMultiPathwayReport
#' @param modulesObj `MultiOmicsModules` class object
#' @param fontsize Size of the font to be used in the plot
#'
#' @return a Heatmap list object from ComplexHeatmap package of the results
#' contained in the `MultiOmicsModules` object provided
#'
#' @examples
#' data(multiOmics)
#' data(reactSmall)
#'
#' genesToUse <- row.names(multiOmics[[1]])
#'
#' MOM_survival <- multiOmicsSurvivalModuleTest(multiOmics, reactSmall[[1]],
#'     survFormula = "Surv(days, status) ~", autoCompleteFormula = TRUE,
#'     useTheseGenes = genesToUse
#' )
#'
#' plotModuleReport(MOM_survival,
#'     MOcolors = c(
#'         exp = "red", met = "green", mut = "blue",
#'         cnv = "yellow"
#'     )
#' )
#'
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom grid grid.text
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export

plotModuleReport <- function(
    modulesObj, MOcolors = NULL, priority_to = NULL, fontsize = 12,
    ...) {
    checkmate::assertClass(modulesObj, "MultiOmicsModules")

    summary <- formatModuleReport(modulesObj)
    rownames(summary) <- paste0(rownames(summary), "\u00B0 module")

    annCol <- guessOmics(colnames(summary))
    omics <- annCol[2:length(annCol)]

    if (is.null(priority_to) & (!is.null(MOcolors))) {
        if (!is.null(names(MOcolors))) {
            priority_to <- names(MOcolors)
        }
    }
    if (is.null(MOcolors)) {
        MOcolors <- guessOmicsColors(omics)
    }
    if (length(MOcolors) != length(unique(omics))) {
        stop(paste0("Length of MOcolors differs from the number of omics: ", 
                    paste(unique(omics), collapse = ", ")))
    }
    if (is.null(names(MOcolors))) {
        names(MOcolors) <- unique(omics)
    }

    omics <- omics[order(match(omics, priority_to))]
    colors <- createColors(omics, MOcolors)
    names(colors) <- unique(omics)
    pvalcol <- colorRamp2(c(0, 1), c("#edf7f5", "#2796bd"))

    msummary <- as.matrix(summary[, 2:ncol(summary)])
    msummary <- order_by_covariates(msummary, 0, priority_to)

    cell_text <- function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", msummary[i, j]), x, y, 
                  gp = gpar(fontsize = fontsize))
    }

    ta <- HeatmapAnnotation(Omics = omics, col = list(Omics = colors))

    ht1 <- Heatmap(
        matrix = summary$pvalue, col = pvalueShades, cluster_rows = FALSE,
        show_heatmap_legend = FALSE, name = "Pvalue", 
        column_names_gp = gpar(fontsize = fontsize),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", summary$pvalue[i]), x, y, 
                      gp = gpar(fontsize = fontsize))}
    )
    dots <- list(...)
    defargs <- list(
        matrix = msummary, name = "p-Value", col = pvalueShades, 
        cluster_rows = FALSE, cluster_columns = FALSE, cell_fun = cell_text, 
        top_annotation = ta, row_names_gp = gpar(fontsize = fontsize),
        column_names_gp = gpar(fontsize = fontsize)
    )
    args <- matchArguments(dots, defargs)

    ht2 <- do.call(Heatmap, args)
    (ht1 + ht2)
}
