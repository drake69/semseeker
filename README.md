
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semseeker

<!-- badges: start -->

    #> ✓ Setting active project to '/Users/lcorsaro/Documents/Progetti_Sviluppo/
    #> semseeker'

[![](https://img.shields.io/badge/devel%20version-0.3.4-blue.svg)](https://github.com/drake69/semseeker)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://codecov.io/gh/drake69/semseeker/branch/main/graph/badge.svg)](https://codecov.io/gh/drake69/semseeker)
[![](https://img.shields.io/github/last-commit/drake69/semseeker.svg)](https://github.com/drake69/semseeker/commits/main)
[![R build
status](https://github.com/drake69/semseeker/workflows/R-CMD-check/badge.svg)](https://github.com/drake69/semseeker/actions)

<!-- [![R-CMD-check](https://github.com/drake69/semseeker/workflows/R-CMD-check/badge.svg)](https://github.com/drake69/semseeker/actions) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/drake69/semseeker/branch/main/graph/badge.svg)](https://app.codecov.io/gh/drake69/semseeker?branch=main) -->
<!-- badges: end -->

The goal of semseeker is to find all methylation localized and enriched
variants.

## Installation

You can install semseeker using devtools, future relese will be
available through CRAN.

Install the latest release:

    install.packages("devtools")
    library("devtools")
    install_github("drake69/semseeker")

## Quick Example

This is a basic example which shows how you can create the beta’s
methylation matrix to use for calculation using ChAMP:

    library(ChAMP)
    idat_folder <- "~/source_idat/"
    resultFolder = "~/result/"
    myLoadN <- champ.load(directory = idat_folder,
                          method = "minfi",
                          methValue="B",
                          autoimpute=TRUE,
                          filterDetP=TRUE,
                          ProbeCutoff=0,
                          SampleCutoff=0.1,
                          detPcut=0.01,
                          filterBeads=TRUE,
                          beadCutoff=0.05,
                          filterNoCG=TRUE,
                          filterSNPs=TRUE,
                          population=NULL,
                          filterMultiHit=TRUE,
                          filterXY=TRUE,
                          force=FALSE,
                          arraytype="450K")

    # normalize with ChAMP
    myNormN<-champ.norm(beta=myLoadN$beta,
                        rgSet=myLoadN$rgSet,
                        mset=myLoadN$mset,
                        resultsDir= resultFolder,
                        method="SWAN",
                        plotBMIQ=FALSE,
                        arraytype="450K",
                        cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
                        )

    saveRDS(myNormN,"~/normalizedData.rds")

This how to obtain the analyzed data:

    library(semseeker)

    normalizedData <- readRDS("~/normalizedData.rds")

    sampleSheet <- read.csv2("~/sample_sheet.csv")

    semseeker (sampleSheet = sampleSheet, 
            methylationData = normalizedData,
            resultFolder = "~/semseeker_result/")

# Complete Example

Look in to the example folder of the repository to seea complete and
working example with data from Gene Expression Omnibus (GEO)

# Input requirements

-   the samplesheet dataframe should contain a column called
    Sample_Group, the admitted values are: Case, Control, Reference. if
    you don0t have the Refernce population you can duplicate the Control
    population rows and use Reference in the Sample_Group column.
