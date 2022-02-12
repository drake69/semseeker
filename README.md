<!-- README.md is generated from README.Rmd. Please edit that file -->

# semseeker
badge_devel("drake69/semseeker", "blue")

<!-- badges: start -->
`r badge_devel("drake69/semseeker", "blue")`
`r badge_lifecycle("stable")`
`r badge_repostatus("Active")`
`r badge_codecov("drake69/semseeker")`
`r badge_last_commit("drake69/semseeker")`
`r badge_github_actions("drake69/semseeker")`
`r badge_codefactor("drake69/semseeker")`
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

Look in to the example folder of the repository to seea complete and working example with data from Gene Expression Omnibus (GEO)

# Input requirements
- the samplesheet dataframe should contain a column called Sample_Group, the admitted values are: Case, Control, Reference. if you don0t have the Refernce population you can duplicate the Control population rows and use Reference in the Sample_Group column.
- the methylationData dataframe should have as columns name the same names in Sample_ID column of the sample sheet.

# Known limit

Actrually semseekwer works with EPIC data source, for data source as 450K and 27K some probes are missed due the changes of manifest.

# The outcomes are
<ul>
<li>
per each population
</li>
<li>
bed graph file with the delta methylation value above and under the
outline threshold
</li>
<li>
bed file of found MUTATIONS due to hyper methylation and hypomethylation
</li>
<li>
bed file of found LESIONS due to hyper methylation and hypomethylation
</li>
<li>
a cumulative bed file for lesions with a column identifying the sample
without annotations - a cumulative bed file for lesions with a column
identifying the sample annotated with genomic area, gene part, island
and DMR
</li>
<li>
chart: heatmaps to compare the burden difference cases vs. control per
genomic area
</li>
<li>
pivots: pivot table to compare the burden difference cases vs. control
per genomic area
</li>
</ul>
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
