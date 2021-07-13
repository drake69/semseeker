#' createPivotResultFromMultipleBedGeneric load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param resultFolder folder as root for bedfiles organized per anomaly
#' @param anomalyLabel anomaly definition used to lable folder and files eg
#' MUTATIONS, LESIONS
#' @param figureLable figure used to create the file HYPO HYPER
#' @param probeFeatures features of probe CHR and START and NAME
#' @param columnLabel
#'
#' @return list of pivot by column identified with columnLabel and by Sample

#'
createPivotResultFromMultipleBedGeneric <- function(resultFolder, genomicAreaMain, genomicAreaSub, sampleSheet) {

  POSITION <- NULL

  fileName <- paste(resultFolder, "/", genomicAreaMain , "._annotatedBed.bed", sep = "")
  annotatedBed <- read_csv(fileName)
  annotatedBed <- subset(annotatedBed, POPULATION !="Reference")
  annotatedBed$GROUP <- as.factor(annotatedBed$GROUP)
  # levels(annotatedBed$GROUP)

  pheno <- read_csv(sampleSheet)

  annotatedBed$ANOMALY <- as.factor(annotatedBed$ANOMALY)
  annotatedBed$FIGURE <- as.factor(annotatedBed$FIGURE)
  hypo <- subset(annotatedBed, FIGURE=="HYPER" & ANOMALY == "LESIONS" & GROUP ==genomicAreaSub )[, c(genomicAreaMain,"SAMPLENAME","POPULATION","freq")]
  hyper <- subset(annotatedBed,FIGURE=="HYPO"  & ANOMALY == "LESIONS"& GROUP ==genomicAreaSub)[, c(genomicAreaMain,"SAMPLENAME","POPULATION","freq")]

  colnames(hypo) <- c("GENE","SAMPLENAME","POPULATION","HYPO")
  colnames(hyper) <- c("GENE","SAMPLENAME","POPULATION","HYPER")

  hypo$freq <- -1 * hypo$HYPO
  merged <- merge(hypo,hyper,all=TRUE)[, c("GENE","SAMPLENAME","POPULATION","HYPO","HYPER")]
  merged$HYPO[is.na(merged$HYPO)] <-0
  merged$HYPER[is.na(merged$HYPER)] <-0
  merged$BALANCE <- merged$HYPO - merged$HYPER
  merged_downregulated <- subset(merged, BALANCE < 0)
  #pivot
  finalResult_down <- reshape2::dcast(data = merged_downregulated,POPULATION+SAMPLENAME ~ GENE, value.var = "BALANCE", sum)

  merged_upregulated <- subset(merged, BALANCE > 0)
  #pivot
  finalResult_up <- reshape2::dcast(data = merged_upregulated,POPULATION+SAMPLENAME ~ GENE, value.var = "BALANCE", sum)

  finalResultdim_up<-dim(finalResult_up)[2]
  pheno_final_up <- merge(pheno[,c("Sample_ID","Sample_Group")],finalResult_up[,2:finalResultdim_up], by.x = "Sample_ID", by.y = "SAMPLENAME")

  finalResultdim_down<-dim(finalResult_down)[2]
  pheno_final_down <- merge(pheno[,c("Sample_ID","Sample_Group")],finalResult_down[,2:finalResultdim_down], by.x = "Sample_ID", by.y = "SAMPLENAME")

  sheets <- list(
    SUMMARY = sampleSheet,
    HYPER_LESIONS = hyperLesions,
    HYPO_LESIONS = hypoLesions
  )
  openxlsx::write.xlsx(
    x = sheets,
    file = fileName,
    asTable = TRUE
  )


}

