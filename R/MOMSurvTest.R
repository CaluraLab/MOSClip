#' @importFrom methods new
#' @importFrom survival Surv

MOMSurvTest <- function(
    genes, omicsObj, survFormula = "Surv(days, status) ~",
    autoCompleteFormula = TRUE, robust = FALSE, include_from_annot = FALSE
) {

    # check if topological method has been used
    for (i in seq_along(omicsObj@ExperimentList@listData)) {
        if (omicsObj@modelInfo[i] == "summarizeWithPca") {
            if (!is.null(omicsObj@specificArgs[[i]]$method)) {
                if (omicsObj@specificArgs[[i]]$method == "topological") {
                  stop("Invalid method for module analysis: topological")
                }
            }
        }
    }

    moView <- createMOMView(omicsObj, genes)
    formula <- survFormula
    coxObj <- createCoxObj(omicsObj@colData, moView)

    add_covs <- unlist(
        lapply(
            moView, function(mo) {
                mo$namesCov
            }
        )
    )

    if (include_from_annot) {
        add_annot_covs <- colnames(coxObj)[!colnames(coxObj) %in%
            c("days", "status")]
        add_covs <- c(add_covs, add_annot_covs)
    }

    if (autoCompleteFormula)
        formula <- paste0(survFormula, paste(add_covs, collapse = "+"))

    if (is.null(coxObj)) {
        scox <- list()
    } else {
        if (robust) {
            scox <- survivalcoxr(coxObj, formula)
        } else {
            scox <- survivalcox(coxObj, formula)
        }
    }

    scox$moView <- moView

    scox
}
