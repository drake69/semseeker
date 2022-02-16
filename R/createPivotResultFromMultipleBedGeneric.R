#' createPivotResultFromMultipleBedGeneric load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param genomicAreaMain is the genomic area to investigate for example the gene
#' @param genomicAreaSub is the part of the genomic area main for example the Body of the gene
#' @param sampleSheet sample sheet to use
#'
#' @return list of pivot by column identified with columnLabel and by Sample
#'
createPivotResultFromMultipleBedGeneric <- function(genomicAreaMain, genomicAreaSub, sampleSheet) {


  fileName <- file_path_build( ssEnv$resultFolderData , c(genomicAreaMain , "annotatedBed"), "csv")
  annotatedBed <- utils::read.csv(fileName)
  annotatedBed <- subset(annotatedBed, annotatedBed$POPULATION !="Reference")
  annotatedBed$GROUP <- as.factor(annotatedBed$GROUP)
  # levels(annotatedBed$GROUP)

  pheno <- utils::read.csv(sampleSheet)

  annotatedBed$ANOMALY <- as.factor(annotatedBed$ANOMALY)
  annotatedBed$FIGURE <- as.factor(annotatedBed$FIGURE)
  hypo <- subset(annotatedBed, annotatedBed$FIGURE=="HYPER" & annotatedBed$ANOMALY == "LESIONS" & annotatedBed$GROUP ==genomicAreaSub )[, c(genomicAreaMain,"SAMPLEID","POPULATION","freq")]
  hyper <- subset(annotatedBed,annotatedBed$FIGURE=="HYPO"  & annotatedBed$ANOMALY == "LESIONS"& annotatedBed$GROUP ==genomicAreaSub)[, c(genomicAreaMain,"SAMPLEID","POPULATION","freq")]

  colnames(hypo) <- c("GENE","SAMPLEID","POPULATION","HYPO")
  colnames(hyper) <- c("GENE","SAMPLEID","POPULATION","HYPER")

  hypo$freq <- -1 * hypo$HYPO
  merged <- merge(hypo,hyper,all=TRUE)[, c("GENE","SAMPLEID","POPULATION","HYPO","HYPER")]
  merged$HYPO[is.na(merged$HYPO)] <-0
  merged$HYPER[is.na(merged$HYPER)] <-0
  merged$BALANCE <- merged$HYPO - merged$HYPER
  merged_downregulated <- subset(merged, merged$BALANCE < 0)
  #pivot
  finalResult_down <- reshape2::dcast(data = merged_downregulated,POPULATION+SAMPLEID ~ GENE, value.var = "BALANCE", sum)

  merged_upregulated <- subset(merged, merged$BALANCE > 0)
  #pivot
  finalResult_up <- reshape2::dcast(data = merged_upregulated,POPULATION+SAMPLEID ~ GENE, value.var = "BALANCE", sum)

  finalResultdim_up<-dim(finalResult_up)[2]
  hyperLesions <- merge(pheno[,c("Sample_ID","Sample_Group")],finalResult_up[,2:finalResultdim_up], by.x = "Sample_ID", by.y = "SAMPLEID")

  finalResultdim_down<-dim(finalResult_down)[2]
  hypoLesions <- merge(pheno[,c("Sample_ID","Sample_Group")],finalResult_down[,2:finalResultdim_down], by.x = "Sample_ID", by.y = "SAMPLEID")

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

