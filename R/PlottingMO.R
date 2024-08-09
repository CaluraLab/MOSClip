#' Plot heatmaps of the pathway by omics
#'
#' Given the pathway, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsPathway pathway object
#' @param sortBy a covariate to sort by
#' @param paletteNames three palettes
#' @param additionalAnnotations optional additional sample annotations
#' @param additionalPaletteNames optional additional colors for annotations
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param withSampleNames create also the samples names
#' @param nrowsHeatmaps magnification respect to annotation of sample (annotations take 1 row)
#'
#' @return NULL
#' 
#' @importFrom checkmate assertClass
#' @importFrom grid gpar grid.newpage grid.draw rectGrob
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation '%v%' draw
#' @importFrom stats relevel
#' @importFrom ggplotify as.ggplot
#' 
#' @export
plotPathwayHeat <- function(pathway, multiOmic, sortBy = NULL,
                            paletteNames = NULL,
                            additionalAnnotations = NULL,
                            additionalPaletteNames = NULL,
                            discr_prop_pca = 0.15, discr_prop_events = 0.05,
                            withSampleNames = TRUE, nrowsHeatmaps = 3,
                            orgDbi = "org.Hs.eg.db") {
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway, multiOmic,
                                      min_prop_pca = discr_prop_pca,
                                      min_prop_events = discr_prop_events)
  
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }
  
  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy = NULL)
  idx <- which(unlist(lapply(annotationFull, class))=="numeric")
  if (length(idx)>0) {
    for (i in idx) {
      annotationFull[,i] <- stats::relevel(
        as.factor(annotationFull[,i]), ref = "1")}
  }
  
  omics <- guessOmics(colnames(annotationFull))
  if(is.null(paletteNames)){
    paletteNames <- names(paletteNames)[1:length(unique(omics))]
    paletteNames <- guessOmicsColors(omics)
  }
  
  if(length(paletteNames) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ",
                paste(unique(omics), collapse = " ,")))
  }
  
  if (is.null(names(paletteNames))) {names(paletteNames) <- unique(omics)}
  
  annotationPalettes <- paletteNames
  
  ann_col <- lapply(colnames(annotationFull), function(name) {
    omic <- guessOmic(name)
    if (!omic %in% names(annotationPalettes)){
      stop(paste0(omic, " omic not found in annotationPalettes"))
    }
    discreteColor <- annotationPalettes[[omic]]
    values <- sort(unique(annotationFull[, name]))
    if (!is.null(levels(values))) {values <- levels(values)}
    if (length(values)==1) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("smart")])
      names(annot) <- values
    } else if (length(values)==2) {
      annot <- as.character(
        MOSpaletteSchema[discreteColor, c("smart", "light")])
      names(annot) <- values
    } else if (length(table(annotationFull[, name]))==3) {
      annot <- as.character(
        MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
      names(annot) <- values
    } else {
      stop("I'm puzzled! Too many values to map")
    }
    annot
  })
  names(ann_col) <- colnames(annotationFull)
  
  if (!is.null(additionalAnnotations)) {
    additionalAnnotations <- matchAnnotations(annotationFull,
                                              additionalAnnotations)
    annotationFull <- cbind(annotationFull, additionalAnnotations)
    
    if (!is.null(additionalPaletteNames)) {
      add_ann_col <- lapply(colnames(additionalAnnotations), function(name) {
        values <- sort(unique(additionalAnnotations[[name]]))
        discreteColor <- additionalPaletteNames[[name]]
        
        if (length(values)==2) {
          annot <- as.character(
            MOSpaletteSchema[discreteColor, c("smart", "light")])
          names(annot) <- values
        } else if (length(table(annotationFull[, name]))==3) {
          annot <- as.character(
            MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
          names(annot) <- values
        } else {
          annot <- getContinousPalette(discreteColor, length(values))
          # names(annot) <- levels(values)
        }
        annot
      })
      names(add_ann_col) <- colnames(additionalAnnotations)
      ann_col=c(ann_col, add_ann_col)
    }
  }
  
  annotationFull <- sortAnnotations(annotationFull, sortBy)
  annotationFull <- annotationFull[, rev(colnames(annotationFull))]
  
  ha <- HeatmapAnnotation(df = annotationFull, col = ann_col)
  ht_list = ha
  for(n in seq_along(involved)){
    heatMatrix <- involved[[n]]$sigModule
    heatMatrix <- heatMatrix[, row.names(annotationFull), drop=F]
    row.names(heatMatrix) <- conversionToSymbols(row.names(heatMatrix), orgDbi)
    Ht <- Heatmap(
      heatMatrix, cluster_columns = F, cluster_rows = F, show_column_names = F,
      heatmap_legend_param = list(title = NULL), col = c("#FFFFFF", 
              rev(as.vector(ann_col[involved[[n]]$covsConsidered][[1]]))))
    ht_list = ht_list %v% Ht}
  suppressMessages(ComplexHeatmap::draw(ht_list, legend_grouping = "original"))
  gb <- grid::grid.grabExpr(ComplexHeatmap::draw(
    ht_list, legend_grouping = "original"))
  gb <- ggplotify::as.ggplot(gb)
  invisible(gb)
}

#' Plot KM of the pathway by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param formula a formula to compute the plot
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#' @param risk.table logical to show risk.table
#' @param pval logical to show pvalue
#' @param size line width of the KM curves
#' @param inYears set time in years
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to include
#' @param additional_continuous names of the additional continous variables to include
#'
#' @return NULL
#'
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#' 
#' @export
plotPathwayKM <- function(pathway, multiOmic,
                          formula = "Surv(days, status) ~ PC1",
                          fileName=NULL, paletteNames = NULL, h = 9, w=7, 
                          risk.table=TRUE, pval=TRUE, size=1, inYears=FALSE,
                          discr_prop_pca=0.15, discr_prop_events=0.05,
                          additional_discrete=NULL, additional_continuous=NULL){
  
  checkmate::assertClass(pathway, "MultiOmicsPathway")
  
  involved <- guessInvolvementPathway(pathway, multiOmic, 
                                      min_prop_pca=discr_prop_pca,
                                      min_prop_events=discr_prop_events)
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  multiOmicObj <- multiOmic
  coxObj <- createCoxObj(multiOmicObj@colData, pathway@pathView)
  daysAndStatus <- coxObj[, c("status", "days"), drop=F]
  
  if (inYears)
    daysAndStatus$days <- daysAndStatus$days/365.24
  
  coxObj <- data.frame(daysAndStatus, 
                       annotationFull[row.names(daysAndStatus), , drop=F])
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  
  palette=NULL
  if (!is.null(paletteNames)) {
    if (length(paletteNames)==1) {
      palette=paletteNames
    } else {
      classes <- names(fit$strata)
      if (length(classes)==length(paletteNames))
        palette = paletteNames
    }
  }
  
  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = risk.table,
                             pval=pval, palette=palette, size=size)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
  # invisible(list(plot = p, fit = fit, coxObj = coxObj))
}

#' Plot heatmaps of the module by omics
#'
#' Given the pathway and the module, it creates the heatmaps of the mostly involved genes for each omic.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param sortBy a covariate to sort by
#' @param paletteNames three palettes
#' @param additionalAnnotations optional additional sample annotations
#' @param additionalPaletteNames optional additional colors for annotations
#' @param withSampleNames create also the samples names
#' @param fontsize_row size of the fonts for rows
#' @param fontsize_col like fontsize_row but for columns
#' @param nrowsHeatmaps magnification respect to annotation of sample (annotations take 1 row)
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation '%v%'
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid gpar grid.newpage grid.draw rectGrob grid.grabExpr
#' @importFrom graphics plot
#' @importFrom stats relevel
#' @importFrom ggplotify as.ggplot
#' 
#' @export
plotModuleHeat <- function(pathway, multiOmic, moduleNumber, sortBy = NULL,
                           fileName = NULL, paletteNames = NULL,
                           additionalAnnotations = NULL,
                           additionalPaletteNames = NULL,
                           withSampleNames = TRUE, fontsize_row = 10,
                           fontsize_col = 1, nrowsHeatmaps = 3,
                           orgDbi = "org.Hs.eg.db", discr_prop_pca = 0.15,
                           discr_prop_events = 0.05) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")
  
  moduleGenes <- pathway@modules[[moduleNumber]]
  
  involved <- guessInvolvement(pathway, multiOmic, moduleNumber = moduleNumber,
                               min_prop_pca = discr_prop_pca,
                               min_prop_events = discr_prop_events)
  
  if(length(paletteNames)!=length(involved)) {
    repTimes <- ceiling(length(involved)/length(paletteNames))
    paletteNames <- rep(paletteNames, repTimes)[seq_along(involved)]
  }
  
  # Create annotation and sort
  annotationFull <- formatAnnotations(involved, sortBy = NULL)
  idx <- which(unlist(lapply(annotationFull, class))=="numeric")
  if (length(idx)>0) {
    for (i in idx) {
      annotationFull[,i] <- stats::relevel(
        as.factor(annotationFull[,i]), ref="1")}
  }
  
  omics <- guessOmics(colnames(annotationFull))
  if(is.null(paletteNames)){
    paletteNames <- names(paletteNames)[1:length(unique(omics))]
    paletteNames <- guessOmicsColors(omics)
  }
  
  if(length(paletteNames) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ",
                paste(unique(omics), collapse = ", ")))
  }
  
  if (is.null(names(paletteNames))) {names(paletteNames) <- unique(omics)}
  
  annotationPalettes <- paletteNames
  
  ann_col <- lapply(colnames(annotationFull), function(name) {
    omic <- guessOmic(name)
    if (!omic %in% names(annotationPalettes)) {
      stop(paste0(omic, " omic not found in annotationPalettes"))
    }
    discreteColor <- annotationPalettes[[omic]]
    values <- sort(unique(annotationFull[, name]))
    if (length(values)==1) {
      annot <- as.character(MOSpaletteSchema[discreteColor, c("smart")])
      names(annot) <- values
    } else if (length(values)==2) {
      annot <- as.character(
        MOSpaletteSchema[discreteColor, c("smart", "light")])
      names(annot) <- values
    } else if (length(table(annotationFull[, name]))==3) {
      annot <- as.character(
        MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
      names(annot) <- values
    } else {
      stop("I'm puzzled! Too many values to map")
    }
    annot
  })
  names(ann_col) <- colnames(annotationFull)
  
  if (!is.null(additionalAnnotations)) {
    additionalAnnotations <- matchAnnotations(annotationFull,
                                              additionalAnnotations)
    annotationFull <- cbind(annotationFull, additionalAnnotations)
    
    if (!is.null(additionalPaletteNames)) {
      add_ann_col <- lapply(colnames(additionalAnnotations), function(name) {
        values <- sort(unique(additionalAnnotations[[name]]))
        discreteColor <- additionalPaletteNames[[name]]
        
        if (length(values)==2) {
          annot <- as.character(
            MOSpaletteSchema[discreteColor, c("smart", "light")])
          names(annot) <- values
        } else if (length(table(annotationFull[, name]))==3) {
          annot <- as.character(
            MOSpaletteSchema[discreteColor, c("dark", "smart", "light")])
          names(annot) <- values
        } else {
          # annot <- getContinousPalette(discreteColor, length(values))
          annot <- colorRamp2(c(min(values), max(values)),
                              getContinousPalette(discreteColor, 2))
          names(annot) <- levels(values)
        }
        annot
      })
      names(add_ann_col) <- colnames(additionalAnnotations)
      ann_col = c(ann_col, add_ann_col)
    }
  }
  
  annotationFull <- sortAnnotations(annotationFull, sortBy)
  annotationFull <- annotationFull[, rev(colnames(annotationFull))]
  
  ha <- HeatmapAnnotation(df = annotationFull, col = ann_col)
  ht_list = ha
  for(n in seq_along(involved)){
    heatMatrix <- involved[[n]]$sigModule
    heatMatrix <- heatMatrix[, row.names(annotationFull), drop=F]
    row.names(heatMatrix) <- conversionToSymbols(row.names(heatMatrix), orgDbi)
    Ht <- Heatmap(heatMatrix, cluster_columns = F, cluster_rows = F,
                 show_column_names = F, 
                 heatmap_legend_param = list(title = NULL),
                 col = c("#FFFFFF", 
                         as.vector(ann_col[[involved[[n]]$covsConsidered[1]]]))
                 )
    ht_list = ht_list %v% Ht}
  suppressMessages(ComplexHeatmap::draw(ht_list, legend_grouping = "original"))
  gb = grid::grid.grabExpr(ComplexHeatmap::draw(
    ht_list, legend_grouping = "original"))
  # invisible(ht_list)
  return(gb)
}


#' Plot KM of the module by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param formula a formula to compute the plot
#' @param fileName optional filenames to save the plot
#' @param paletteNames three palettes
#' @param h the height of the plot
#' @param w the width of the plot
#' @param risk.table logical to show risk.table
#' @param pval logical to show pvalue
#' @param size line width of the KM curves
#' @param inYears set time in years
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#' @param additional_discrete names of the additional discrete variables to include
#' @param additional_continuous names of the additional continous variables to include
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra arrangeGrob
#' @importFrom survminer ggsurvplot surv_fit
#' @importFrom ggplot2 ggsave
#' 
#' @export
plotModuleKM <- function(pathway, multiOmic, moduleNumber,
                         formula = "Surv(days, status) ~ PC1",
                         fileName = NULL, paletteNames = NULL, h = 9, w = 7,
                         risk.table = TRUE, pval = TRUE, size = 1,
                         inYears = FALSE,
                         discr_prop_pca = 0.15, discr_prop_events = 0.05,
                         additional_discrete = NULL,
                         additional_continuous = NULL) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")
  
  involved <- guessInvolvement(pathway, multiOmic, moduleNumber = moduleNumber,
                               min_prop_pca = discr_prop_pca,
                               min_prop_events = discr_prop_events)
  
  multiOmicObj <- multiOmic
  coxObj <- createCoxObj(
    multiOmicObj@colData, pathway@modulesView[[moduleNumber]])
  annotationFull <- formatAnnotations(involved, sortBy=NULL)
  
  days_status_names <- c("status", "days")
  
  daysAndStatus <- coxObj[, days_status_names, drop=F]
  if (inYears)
    daysAndStatus$days <- daysAndStatus$days/365.24
  
  additional_clinic = NULL
  if (!is.null(additional_discrete)){
    not_found <- setdiff(additional_discrete, colnames(coxObj))
    if (length(not_found) > 0)
      stop(paste0("Some discrete variables were not found: ", 
                  paste(not_found, collapse = ", ", "\n"),
                  "We found the following variables: ", 
                  paste(colnames(coxObj), collapse = ", ")))
    
    fixed_covs <- setdiff(additional_discrete, colnames(annotationFull))
    additional_clinic <- coxObj[row.names(daysAndStatus), fixed_covs, drop=F]
  }
  
  if (!is.null(additional_continuous)){
    not_found <- setdiff(additional_continuous, colnames(coxObj))
    if (length(not_found) > 0)
      stop(paste0("Some continuous variables were not found: ", 
                  paste(not_found, collapse = ", ", "\n"),
                  "We found the following variables: ", 
                  paste(colnames(coxObj), collapse = ", ")))
    
    fixed_covs <- setdiff(additional_continuous, colnames(annotationFull))
    df <- cbind(
      daysAndStatus, coxObj[row.names(daysAndStatus), fixed_covs, drop=F])
    discretized_covs <- createDiscreteClasses(df, fixed_covs)
    
    if (!is.null(additional_clinic)){
      additional_clinic <- cbind(
        additional_clinic, discretized_covs[, fixed_covs, drop=F])
    } else {
      additional_clinic <- discretized_covs[, fixed_covs, drop=F]
    }
  }
  
  if (is.null(additional_clinic)){
    coxObj <- data.frame(daysAndStatus, 
                         annotationFull[row.names(daysAndStatus), , drop=F])
  } else {
    coxObj <- data.frame(daysAndStatus, 
                         annotationFull[row.names(daysAndStatus), , drop=F],
                         additional_clinic)
  }
  
  fit <- survminer::surv_fit(formula(formula), data = coxObj)
  
  palette=NULL
  if (!is.null(paletteNames)) {
    if (length(paletteNames)==1) {
      palette=paletteNames
    } else {
      classes <- names(fit$strata)
      if (length(classes)==length(paletteNames))
        palette = paletteNames
      if (!is.null(names(paletteNames))) {
        diff = setdiff(names(paletteNames), classes)
        if (length(diff)!=0)
          stop("Names of paletteNames must be equal to classes")
      }
    }
  }

  p <- survminer::ggsurvplot(fit, data = coxObj, risk.table = risk.table,
                             pval=pval, palette=unname(palette), size=size)
  
  if(!is.null(fileName)) {
    ggplot2::ggsave(filename = fileName, p, height = h, width = w)
  } else {
    p
  }
}

#' Plot graph of the module by omics
#'
#' Given the pathway, it creates the Kaplan-meier curves following the formula.
#'
#' @param pathway MultiOmicsModule pathway object
#' @param moduleNumber a module number
#' @param orgDbi if needed, a organism Dbi to translate vectors
#' @param legendLabels set up your favourite names for the omics
#' @param paletteNames named vector of MOpalettes, names replace makeLegend arguments
#' @param fileName optional filenames to save the plot
#' @param discr_prop_pca the minimal proportion to compute the pca classes
#' @param discr_prop_events the minimal proportion to compute the event classes
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom igraph V V<- simplify graph_from_graphnel
#' @importFrom AnnotationDbi select
#' @importFrom graphics plot legend
#' @importFrom grDevices dev.off pdf rainbow
#' 
#' @export
plotModuleInGraph <- function(pathway, multiOmic, reactObj, moduleNumber,
                              orgDbi="org.Hs.eg.db",
                              paletteNames=NULL, legendLabels=NULL,
                              fileName=NULL, discr_prop_pca=0.15,
                              discr_prop_events=0.05) {
  
  checkmate::assertClass(pathway, "MultiOmicsModules")
  
  # dentro pathway dovrebbe esserci l'oggetto graphNEL
  net <- igraph::graph_from_graphnel(
    convertPathway(reactObj[[pathway@title]], NULL))
  moduleGenes <- pathway@modules[[moduleNumber]]
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)
  color <- rep("grey", length(V(net)))
  color[names(V(net)) %in% moduleGenes] <- "tomato"
  involved <- guessInvolvement(pathway, multiOmic, moduleNumber = moduleNumber,
                               min_prop_pca=discr_prop_events,
                               min_prop_events=discr_prop_events)
  mark.groups=lapply(involved, function(x) {
    row.names(x$subset)
  })
  
  group.names <- sapply(involved, function(x) {
    guessOmic(x$covsConsidered)
  })
  
  colLength <- length(mark.groups)
  if (colLength<3) {
    mark.col=rainbow(3, alpha=0.33)[seq_len(colLength)]
  } else {
    mark.col=rainbow(colLength, alpha=0.33)
  }
  mark.border=NA
  
  if (!is.null(paletteNames)) {
    # if (is.null(names(paletteNames)))
    #   stop("paletteNames must be named vector")
    
    if (!is.null(names(paletteNames))) {
      mismatch <- setdiff(group.names, names(paletteNames))
      if (length(mismatch)>0)
        stop(paste0("Missing palette for omics:" , paste(mismatch,
                                                       collapse = ", ")))
      paletteNames <- paletteNames[group.names]
    }
    
    # legendLabels <- names(paletteNames)
    err <- setdiff(paletteNames, row.names(MOSpaletteSchema))
    if (length(err)!=0)
      stop(paste0(err, " paletteNames value is not allowed."))
    
    mark.col <- MOSpaletteSchema[paletteNames, ]$transparent
  }
  
  if (is.null(legendLabels)) {
    # legendLabels <- c(paste("omic",seq_len(length(pathway@modulesView[[moduleNumber]]))))
    legendLabels <- group.names
  } else {
    if (length(legendLabels)!=length(group.names))
      warning("Some legendLabels were not found in data")
    legendLabels <- legendLabels[seq_along(group.names)]
  }
  
  labels <- conversionToSymbols(names(V(net)), orgDbi)
  
  if (!is.null(fileName)) {
    pdf(fileName)
  }
  plot(net, edge.arrow.size=.5, edge.curved=.2,
       vertex.label=labels, vertex.label.cex=.6, vertex.label.family="sans", 
       vertex.label.font=2, vertex.color=color, vertex.frame.color="gray", 
       vertex.label.color="black", vertex.size=15,
       mark.groups=mark.groups, mark.col=mark.col, mark.border=NA
  )
  legend(x=-1, y=-1, legendLabels, pch=21, horiz=TRUE,
         col="#777777", pt.bg=mark.col, pt.cex=2, cex=.8, bty="n", ncol=1)
  #if (!is.null(fileName)) {
   # dev.off()
  #}
}

#' Summarize and plot pathways' info from a list of MultiOmicsPathway (MOP)
#'
#' Given the list of MOPs, it plots the table.
#'
#' @param multiPathwayList MultiOmicsPathway list pathway object
#' @param top use top number of pathways
#' @param MOcolors character vector with the omic colors.
#' The colors should be among the colors in \code{showMOSpalette()}
#' @param priority_to a vector with the covariates (omic name) that should go first
#' @param \dots additional argument to be passed to pheatmap
#' 
#'
#' @return NULL
#'
#' @importFrom pheatmap pheatmap
#' 
#' @export
plotMultiPathwayReport <- function(multiPathwayList, top=25, MOcolors=NULL,
                                   priority_to=NULL, fontsize=6, ...){
  
  if(!is.list(multiPathwayList)) {stop("multiPathwayList must be a list.")}
  
  summary <- multiPathwayReport(multiPathwayList)
  top <- min(top, NROW(summary))
  
  annCol <- guessOmics(colnames(summary))
  omics <- annCol[2:length(annCol)]
  
  if (is.null(priority_to) & (!is.null(MOcolors))){
    if (!is.null(names(MOcolors)))
      priority_to = names(MOcolors)
  }
  
  if(is.null(MOcolors)) {MOcolors <- guessOmicsColors(omics)}
  
  if(length(MOcolors) != length(unique(omics))){
    stop(paste0("Length of MOcolors differs from the number of omics: ",
                paste(unique(omics), collapse = ", ")))
  }
  
  if (is.null(names(MOcolors))) {names(MOcolors) <- unique(omics)}
  
  omics <- omics[order(match(omics, priority_to))]
  colors <- createColors(omics, MOcolors)
  names(colors) <- unique(omics)
  
  msummary <- as.matrix(summary[seq_len(top),2:ncol(summary)])
  msummary <- order_by_covariates(msummary, 0, priority_to)
  
  cell_text <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", msummary[i, j]), x, y,
              gp = gpar(fontsize = fontsize))}
  
  ta <- HeatmapAnnotation(Omics = omics, col = list(Omics = colors))
  
  ht1 <- Heatmap(matrix = summary$pvalue[seq_len(top)],
                 col = pvalueShades[1:round((summary$pvalue[top]*100)+1)],
                 cluster_rows = F, show_heatmap_legend = F, name = "Pvalue",
                 column_names_gp = gpar(fontsize = fontsize),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", summary$pvalue[i]), x, y,
                             gp = gpar(fontsize = fontsize))
                 })
  dots <- list(...)
  defargs <- list(matrix = msummary, name = "p-Value", col = pvalueShades,
                  cluster_rows = F, cluster_columns = F, cell_fun = cell_text,
                  top_annotation = ta, row_names_gp = gpar(fontsize = fontsize),
                  column_names_gp = gpar(fontsize = fontsize))
  args <- matchArguments(dots, defargs)
  
  ht2 <- do.call(Heatmap, args)
  suppressWarnings(ht1 + ht2)
}

#' Summarize and plot pathways' info from a MultiOmicsModule (MOM) object
#'
#' Given a MOM, it plots the table.
#'
#' @inheritParams plotMultiPathwayReport
#' @param pathwayObj MultiOmicsModule of pathway object
#' 
#'
#' @return NULL
#' @importFrom checkmate assertClass
#' @importFrom pheatmap pheatmap
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom grid grid.text
#' @importFrom circlize colorRamp2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
plotModuleReport <- function(pathwayObj, MOcolors=NULL, priority_to=NULL,
                             fontsize = 12, ...) {
  
  checkmate::assertClass(pathwayObj, "MultiOmicsModules")
  
  summary <- formatModuleReport(pathwayObj)
  rownames(summary) <- paste0(rownames(summary), "° module")
  
  annCol <- guessOmics(colnames(summary))
  omics <- annCol[2:length(annCol)]
  
  if (is.null(priority_to) & (!is.null(MOcolors))) {
    if (!is.null(names(MOcolors)))
      priority_to = names(MOcolors)
  }
  if (is.null(MOcolors)) {MOcolors <- guessOmicsColors(omics)}
  if (length(MOcolors) != length(unique(omics))) {
    stop (paste0("Length of MOcolors differs from the number of omics: ",
                 paste(unique(omics), collapse = ", ")))
  }
  if (is.null(names(MOcolors))){names(MOcolors) <- unique(omics)}
  
  omics <- omics[order(match(omics,priority_to))]
  colors <- createColors(omics, MOcolors)
  names(colors) <- unique(omics)
  pvalcol <- colorRamp2(c(0,1), c("#edf7f5", "#2796bd"))
  
  msummary <- as.matrix(summary[,2:ncol(summary)])
  msummary <- order_by_covariates(msummary, 0, priority_to)
  
  cell_text <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", msummary[i, j]), x, y,
              gp = gpar(fontsize = fontsize))}
  
  ta <- HeatmapAnnotation(Omics = omics, col = list(Omics = colors))
  
  ht1 <- Heatmap(matrix = summary$pvalue, col = pvalueShades, cluster_rows = F,
                 show_heatmap_legend = F, name = "Pvalue",
                 column_names_gp = gpar(fontsize = fontsize),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", summary$pvalue[i]), x, y,
                             gp = gpar(fontsize = fontsize))
                 })
  dots <- list(...)
  defargs <- list(matrix = msummary, name = "p-Value", col = pvalueShades,
                  cluster_rows = F, cluster_columns = F, cell_fun = cell_text,
                  top_annotation = ta, row_names_gp = gpar(fontsize = fontsize),
                  column_names_gp = gpar(fontsize = fontsize))
  args <- matchArguments(dots, defargs)
  
  ht2 <- do.call(Heatmap, args)
  suppressWarnings(ht1 + ht2)
}
