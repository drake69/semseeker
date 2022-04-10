install.packages("devtools")
library("devtools")
install_github("drake69/semseeker")
library("semseeker")
packageVersion('semseeker')

devtools::install_github("hrbrmstr/dtupdate")
library(dtupdate)
github_update()

workingFolder <- file.path("~/Downloads/GSE139307")
dir.create(workingFolder)

sampleSheet <- buildDataSetFromGEO("GSE139307",workingFolder, 0)

# ChAMP need the sample name variable as first column
# so let's move the Sample_ID as first column
sampleSheet <- data.frame("Sample_ID"=sampleSheet$Sample_ID, sampleSheet[, colnames(sampleSheet)!="Sample_ID"])

write.table(
  sampleSheet,
  paste(workingFolder, "/", "final_samplesheet.csv", sep = ""),
  row.names = FALSE,
  sep=",",
  quote = FALSE
)

library(ChAMP)
idat_folder <- workingFolder
resultFolder = file.path( workingFolder,"/result/")

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
                      force=TRUE,
                      arraytype="450K")

# normalize with ChAMP
normalizedData<-champ.norm(beta=myLoadN$beta,
                           rgSet=myLoadN$rgSet,
                           mset=myLoadN$mset,
                           resultsDir= resultFolder,
                           method="SWAN",
                           plotBMIQ=FALSE,
                           arraytype="450K",
                           cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
)

sampleSheet <- read.csv("~/Downloads/GSE139307/final_samplesheet.csv")
normalizedData <- readRDS("~/Downloads/GSE139307/normalizeddata.rds")


# we need the Sample_Group, the sample sheet of the example has the column pathology we can use to distinguish Cases vs Controls
sampleSheet$Sample_Group <- sampleSheet$dioxin.group.ch1
sampleSheet[sampleSheet$Sample_Group=="CONTROL", "Sample_Group"] <- "Control"
sampleSheet[sampleSheet$Sample_Group!="Control", "Sample_Group"] <- "Case"

# we have not the reference Sample_Group so we reuse the Control sample group as reference
reference <- subset(sampleSheet,Sample_Group=="Control")
reference$Sample_Group <- "Reference"
sampleSheet <- rbind(sampleSheet, reference)



semseeker (sampleSheet = sampleSheet,
           methylationData = normalizedData,
           resultFolder = file.path(workingFolder,"/semseeker_result/"))
