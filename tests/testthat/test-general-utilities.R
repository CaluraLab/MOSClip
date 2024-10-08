
test_that(
    "extractPvalues", {
        fakePvalues <- lapply(
            1:10, function(x) return(
                list(
                  pvalue = runif(
                    1, 0,
                    1
                )
              )
            )
        )
        naPvalue <- list(pvalue = NA)
        nullPvalue <- list(pvalue = NULL)

        res = sapply(fakePvalues, extractPvalues)
        expect_true(is.na(extractPvalues(naPvalue)))
        expect_true(is.na(extractPvalues(nullPvalue)))
        expect_equal(
            fakePvalues[[1]]$pvalue,
            extractPvalues(fakePvalues[[1]])
        )
        expect_equal(res, unname(unlist(fakePvalues)))
    }
)


test_that(
    "conversionToSymbols", {
        vec <- "SYMBOL:MMP15"
        vec1 <- c("SYMBOL:MMP15", "ENTREZID:10", "SYMBOL:CMA1")  #mixed ids
        vec2 <- c("SYMBOL:MMP15", "SYMBOL:MMP24", "SYMBOL:CMA1")
        vec3 <- c("ENTREZID:1", "ENTREZID:10", "ENTREZID:100")
        vec4 <- c("A1BG", "NAT2", "SYMBOL:ADA")  #not all are graphite style
        vec5 <- c(
            "UNIPROT:Q86VV6", "UNIPROT:Q9Y5R2", "UNIPROT:P23946",
            "UNIPROT:Q4FEB3"
        )  #many2one
        vec6 <- c(
            "SYMBOL:MMP15", "SYMBOL:ENRICA", "SYMBOL:CMA1",
            "SYMBOL:PAOLO"
        )  #symbols not in dictionary
        vec7 <- c(
            "ENTREZID:1", "ENTREZID:ENRICA", "ENTREZID:100",
            "ENTREZID:PAOLO"
        )  #entrez not in dictionary
        vec8 <- NULL

        expect_identical(
            conversionToSymbols(vec, "org.Hs.eg.db"),
            "MMP15"
        )
        expect_identical(
            conversionToSymbols(vec1, "org.Hs.eg.db"),
            c("SYMBOL:MMP15", "ENTREZID:10", "SYMBOL:CMA1")
        )
        expect_identical(
            conversionToSymbols(vec2, "org.Hs.eg.db"),
            c("MMP15", "MMP24", "CMA1")
        )
        expect_identical(
            conversionToSymbols(vec3, "org.Hs.eg.db"),
            c("A1BG", "NAT2", "ADA")
        )
        expect_identical(
            conversionToSymbols(vec4, "org.Hs.eg.db"),
            c("A1BG", "NAT2", "SYMBOL:ADA")
        )
        expect_identical(
            conversionToSymbols(vec5, "org.Hs.eg.db"),
            c("MMP24", "MMP24", "CMA1", "CMA1")
        )
        expect_identical(
            conversionToSymbols(vec6, "org.Hs.eg.db"),
            c("MMP15", "ENRICA", "CMA1", "PAOLO")
        )
        expect_identical(
            conversionToSymbols(vec7, "org.Hs.eg.db"),
            c("A1BG", "ENRICA", "ADA", "PAOLO")
        )
        expect_true(is.null(conversionToSymbols(vec8, "org.Hs.eg.db")))
    }
)

test_that(
    "extractPositivePortion", {
        expect_identical(
            extractPositivePortion(
                matrix(
                  c(1, 1, -1, 2),
                  2
              )
            ),
            matrix(
                c(1, 1, 0, 2),
                2
            )
        )
        expect_identical(
            extractPositivePortion(
                matrix(
                  c(1, 1, -1, 1),
                  2
              ),
                invert = T
            ),
            matrix(
                c(0, 0, 1, 0),
                2
            )
        )
    }
)

test_that(
    "check_minimal_proportion", {
        cls <- c(
            rep(0, 10),
            rep(1, 90)
        )
        expect_true(check_minimal_proportion(cls, min_prop = 0.1))
        cls[1] <- 1
        expect_false(check_minimal_proportion(cls, min_prop = 0.1))
        cls <- c(
            rep(0, 90),
            rep(20, 10)
        )
        expect_true(check_minimal_proportion(cls, min_prop = 0.1))
        cls[100] <- 0
        expect_false(check_minimal_proportion(cls, min_prop = 0.1))
    }
)




test_that(
    "filterMultiOmicsForSamples",
    {
        ex <- matrix(
            1:100, 10, 10, dimnames = list(
                letters[1:10],
                LETTERS[1:10]
            )
        )



        MO <- makeOmics(
            experiments = ExperimentList(e = ex, df = ex),
            modelInfo = c(
                "summarizeWithPca",
                "summarizeWithPca"
            ),
            specificArgs = list("a", "b")
        )

        smallMO <- filterMultiOmicsForSamples(MO, c("A"))

        expect_identical(
            smallMO@ExperimentList[[1]],
            as.matrix(smallMO@ExperimentList[[2]])
        )
        expect_error(filterMultiOmicsForSamples(MO, c("Z")))
    }
)


test_that(
    "createDataModule work", {
        mo <- fake_mo(type = "survival")
        mov <- createMOMView(mo, paste0("gene_", seq_len(200)))
        expect_length(mov, 4)
        mov[[1]]$omicName = "test"
        expect_error(
            createDataModule(mov[[1]], mo),
            paste0(
                "omicName not found in ExperimentList.\n",
                "Names of experiments in ExperimentList should match",
                "the name arguments given in specificArgs"
            )
        )
    }
)


test_that(
    "showMOSpalette work", {
        expect_no_error(showMOSpalette())
    }
)
