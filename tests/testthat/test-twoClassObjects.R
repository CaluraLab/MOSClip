test_that(
    "Two Class Pathway Test works", {
        genes <- paste0(
            "ENTREZID:", c(
                10000, 90427, 5366, 23368, 637,
                3002, 79792, 840, 5414, 5599, 8655,
                140735, 708
            )
        )
        dummy_react <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
        dummy_Omics <- fake_mo(genes = genes, type = "two-classes")
        dummy_annot <- dummy_colData(type = "two-classes")

        twoCPathwayObj <- multiOmicsTwoClassPathwayTest(
            dummy_Omics, dummy_react[[1]], dummy_annot,
            baseFormula = "classes ~", nullModel = "classes ~ 1"
        )

        expect_s4_class(twoCPathwayObj, "MultiOmicsPathway")
        expect_true(
            !is.null(twoCPathwayObj@pvalue) &
                is.numeric(twoCPathwayObj@pvalue)
        )
        expect_true(
            !is.null(twoCPathwayObj@zlist) &
                is.numeric(twoCPathwayObj@zlist)
        )
        expect_true(is.list(twoCPathwayObj@pathView))
        expect_true(
            all(
                vapply(
                  twoCPathwayObj@pathView, function(x) !is.null(x),
                  logical(1)
              )
            )
        )

        test_annot <- rbind("A", dummy_annot)
        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Mismatch in the number of samples"
        )

        test_annot <- dummy_annot
        test_annot[1, ] <- "C"
        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            paste0(
                "Classes should be only two. ",
                "Check your dependent variables columns"
            )
        )

        test_annot <- dummy_annot
        test_order <- sample(
            row.names(test_annot),
            200, replace = FALSE
        )
        test_annot <- test_annot[test_order, , drop = FALSE]
        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            NA
        )

        test_annot <- dummy_annot[-1, , drop = FALSE]
        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Mismatch in the number of samples"
        )

        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "test ~ 1"
            ),
            "Invalid formula. Class column not found in classAnnot"
        )

        expect_error(
            multiOmicsTwoClassPathwayTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1",
                useTheseGenes = c("test1", "test2")
            ),
            "There is no nodes on the graph"
        )

        test_mo <- fake_mo(type = "two-classes")
        expect_error(
            multiOmicsTwoClassPathwayTest(
                test_mo, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Genes not found in some experiments of the omics object"
        )
    }
)



test_that(
    "Two Class Module Test works", {
        genes <- paste0(
            "ENTREZID:", c(
                10000, 90427, 5366, 23368, 637,
                3002, 79792, 840, 5414, 5599, 8655,
                140735, 708
            )
        )
        dummy_react <- readRDS(test_path("fixtures", "reactBiggerDummy.rds"))
        dummy_Omics <- fake_mo(genes = genes, type = "two-classes")
        dummy_annot <- dummy_colData(type = "two-classes")

        twoCModuleObj <- suppressWarnings(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            )
        )

        expect_s4_class(twoCModuleObj, "MultiOmicsModules")
        expect_true(
            !is.null(twoCModuleObj@alphas) &
                is.numeric(twoCModuleObj@alphas)
        )
        expect_true(
            !is.null(twoCModuleObj@zlists) &
                is.list(twoCModuleObj@zlists)
        )

        graph <- convertPathway(dummy_react[[1]], NULL)
        cliques <- extractCliquesFromDag(graph)

        expect_true(is.list(twoCModuleObj@modulesView))
        expect_identical(
            length(cliques),
            length(twoCModuleObj@modulesView)
        )
        expect_true(
            all(
                vapply(
                  twoCModuleObj@modulesView, function(x) !is.null(x),
                  logical(1)
              )
            )
        )

        test_annot <- rbind("A", dummy_annot)
        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Mismatch in the number of sample"
        )

        test_annot <- dummy_annot
        rownames(test_annot)[1] <- "test"
        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Mismatch in covariates and annotations row names."
        )

        test_annot <- dummy_annot
        test_annot[1, ] <- "C"
        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            paste0(
                "Classes should be only two. ",
                "Check your dependent variables columns"
            )
        )

        test_annot <- dummy_annot
        test_order <- sample(
            row.names(test_annot),
            200, replace = FALSE
        )
        test_annot <- test_annot[test_order, , drop = FALSE]
        expect_error(
            suppressWarnings(
                multiOmicsTwoClassModuleTest(
                  dummy_Omics, dummy_react[[1]],
                  test_annot, baseFormula = "classes ~",
                  nullModel = "classes ~ 1"
              )
            ),
            NA
        )

        test_annot <- dummy_annot[-1, , drop = FALSE]
        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], test_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Mismatch in the number of samples"
        )

        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "test ~ 1"
            ),
            "Invalid null formula. Class column not found in classAnnot"
        )

        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, genes, dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "graph argument should be a graphNEL object, not a gene list."
        )

        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1",
                useTheseGenes = c("test1", "test2")
            ),
            paste0(
                "There is no intersection between expression feature ",
                "names and the node names in the graph."
            )
        )

        dummy_Omics@specificArgs$exp$method = "topological"
        expect_error(
            multiOmicsTwoClassModuleTest(
                dummy_Omics, dummy_react[[1]], dummy_annot,
                baseFormula = "classes ~", nullModel = "classes ~ 1"
            ),
            "Invalid method for module analysis: topological"
        )

    }
)

