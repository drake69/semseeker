
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semseeker

<!-- badges: start -->


[![](https://img.shields.io/badge/devel%20version-0.8.6-blue.svg)](https://github.com/drake69/semseeker)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://codecov.io/gh/drake69/semseeker/branch/main/graph/badge.svg)](https://codecov.io/gh/drake69/semseeker)
[![](https://img.shields.io/github/last-commit/drake69/semseeker.svg)](https://github.com/drake69/semseeker/commits/main)
[![R build
status](https://github.com/drake69/semseeker/workflows/R-CMD-check/badge.svg)](https://github.com/drake69/semseeker/actions)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Semseeker aims to identify all variants that are enriched and localized
in methylation.

## Installation

To install semseeker, you can use devtools; upcoming releases will be
accessible via CRAN.

Install the latest release:

    install.packages("devtools")
    library("devtools")
    install_github("drake69/semseeker")

## Quick Example

This basic example demonstrates the process of creating a methylation
matrix for beta values that can be utilized in calculations using ChAMP:

    library(ChAMP)
    idat_folder <- "~/source_idat/"
    result_folder = "~/result/"
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
                        resultsDir= result_folder,
                        method="SWAN",
                        plotBMIQ=FALSE,
                        arraytype="450K",
                        cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
                        )

    saveRDS(myNormN,"~/normalizedData.rds")

Here’s how you can obtain the analyzed data:

    library(semseeker)

    normalizedData <- readRDS("~/normalizedData.rds")

    sample_sheet <- read.csv2("~/sample_sheet.csv")

    semseeker (sample_sheet = sample_sheet, 
            methylation_data = normalizedData,
            result_folder = "~/semseeker_result/")

# Complete Example

You can find a complete and functional example, which includes data from
Gene Expression Omnibus (GEO), by examining the repository’s “example”
folder.

# Input requirements

- The “samplesheet” dataframe must include a column named “Sample_Group”
  with the following accepted values: “Case”, “Control”, and
  “Reference”. If you do not have a Reference population, you may
  duplicate the Control population rows and use “Reference” in the
  “Sample_Group” column.
