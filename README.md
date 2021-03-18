
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semseeker

<!-- badges: start -->
<!-- badges: end -->

The goal of semseeker is to find all methylation localized and enriched
variants

## Installation

You can install semseeker using devtools, future relese will be
available through CRAN.

Install the latest release:

``` r
install.packages("devtools")
library("devtools")
install_github("drake69/semseeker")
```

## Example

This is a basic example which shows how you can create the beta’s
methylation matrix to use for calculation using ChAMP:

``` r
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
```

This how to obtain the analyzed data:

``` r
library(semseeker)

normalizedData <- readRDS("~/normalizedData.rds")

semseeker (sampleSheetPath = "~/source_idat/sample_sheet.csv", 
        methylationData = normalizedData,
        resultFolder = "~/semseeker_result/")
```
Requirements:
- the sample sheet should contain:
  - Sample_Name for the column of the sample name
  - Sample_Group for the dolumn of population
The outcomes are: 
- per each population 
- bed graph file with the delta methylation value above and under the outline threshold 
- bed file of found MUTATIONS due to hyper methylation and hypomethylation 
- bed file of found LESIONS due to hyper methylation and hypomethylation 
- a cumulative bed file for lesions with a column identifying the sample without annotations
- a cumulative bed file for lesions with a column identifying the sample annotated with genomic area, gene part, island and DMR - chart: heatmaps to compare the burden difference cases vs. control per genomic area

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
