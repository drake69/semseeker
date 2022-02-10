
populationCheck <- function(sampleSheet, methylationData)
{
  sampleSheet <- as.data.frame(sampleSheet)
  sampleSheet <- sampleSheet[,!(colnames(sampleSheet) %in% c("Probes_Count", "MUTATIONS_HYPER", "LESIONS_HYPER", "MUTATIONS_HYPO", "LESIONS_HYPO", "MUTATIONS_BOTH", "LESIONS_BOTH"))]

  result <- NULL


  needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sampleSheet))]

  if (length(missedColumns) > 0) {
    result <- paste(result,  " Lost following columns ", missedColumns," ",Sys.time(), "Especting a column with name Sample_Group and possible values Reference,  Control and Case")
  }

  #Exist at least 3 samples per Sample_Group


  if(sum(is.na(sampleSheet$Sample_ID)))
  {
    result <- paste(result,  (" some samples have Sample_ID with NA !"))
  }

  if(sum(is.na(sampleSheet$Sample_Group)))
  {
    result <- paste(result,  (" some samples have Sample_Group with NA !"))
  }

  if (sum(colnames(methylationData) %in% sampleSheet$Sample_ID)!=ncol(methylationData))
  {
    result <- paste(result,  (" The methylation data has not the column's names as expected from the sample sheet in column Sample_ID!"))
  }
  populationGroups <- as.factor(sampleSheet$Sample_Group)
  # expected R as reference, S as Study, C as Control

  # reference population
  populations <-  c("Reference","Control","Case")
  sampleSheet$Sample_Group <- R.utils::toCamelCase(tolower(sampleSheet$Sample_Group), capitalize=TRUE)
  sampleSheet$Sample_Group <- as.factor(sampleSheet$Sample_Group)
  matchedPopulation <- levels(sampleSheet$Sample_Group) %in% populations
  if (is.element(FALSE, matchedPopulation)) {
    result <- paste(result,  " File:",sampleSheetPath, " Sample_Group should contain only: Reference, Control, Case" )
  }

  return(result)
}
