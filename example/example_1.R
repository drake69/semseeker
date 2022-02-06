#
# Example of semseeker application with a GEO data set
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186766

install.packages("devtools")
library("devtools")
install_github("drake69/semseeker")

library("semseeker")

workingFolder <- file.path(getwd(),"/tmp")
sampleSheet <- buildDataSetFromGEO("GSE186766",workingFolder, TRUE)

# we need the Sample_Group, the sample sheet of the example has the column pathology we can use to distinguish Cases vs Controls

sampleSheet <- read.table("~/Documents/Progetti_Sviluppo/semseeker/tmp/final_samplesheet.csv", sep="\t", header = TRUE)

sampleSheet$Sample_Group <- sampleSheet$`pathology:ch1`
sampleSheet[sampleSheet$Sample_Group=="AN", "Sample_Group"] <- "Case"
sampleSheet[sampleSheet$Sample_Group=="control", "Sample_Group"] <- "Control"

reference <- sampleSheet[sampleSheet$Sample_Group=="Control",]
reference[, "Sample_Group"] <- "Reference"
sampleSheet <- rbind(sampleSheet, reference)

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
                      force=FALSE,
                      arraytype="EPIC")

# normalize with ChAMP
normalizedData<-champ.norm(beta=myLoadN$beta,
                    rgSet=myLoadN$rgSet,
                    mset=myLoadN$mset,
                    resultsDir= resultFolder,
                    method="SWAN",
                    plotBMIQ=FALSE,
                    arraytype="EPIC",
                    cores= detectCores(all.tests = FALSE, logical = TRUE) - 1
)

# saveRDS(normalizedData, file.path(workingFolder,"/normalizedData.rds"))
# normalizedData <- readRDS("~/normalizedData.rds")

semseeker (sampleSheet = sampleSheet,
           methylationData = normalizedData,
           resultFolder = file.path(workingFolder,"~/semseeker_result/"))
