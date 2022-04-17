
population_check <- function(sample_sheet, methylation_data)
{
  sample_sheet <- as.data.frame(sample_sheet)
  sample_sheet <- sample_sheet[,!(colnames(sample_sheet) %in% c("Probes_Count", "MUTATIONS_HYPER", "LESIONS_HYPER", "MUTATIONS_HYPO", "LESIONS_HYPO", "MUTATIONS_BOTH", "LESIONS_BOTH"))]

  result <- NULL


  needColumns <- c("Sample_ID", "Sample_Group")
  missedColumns <- needColumns[!(needColumns %in% colnames(sample_sheet))]

  if (length(missedColumns) > 0) {
    result <- paste(result,  " Lost following columns ", missedColumns," ",Sys.time(), "Especting a column with name Sample_Group and possible values Reference,  Control and Case")
  }

  if(sum(is.na(sample_sheet$Sample_ID)))
  {
    result <- paste(result,  (" some samples have Sample_ID with NA !"))
  }

  if(sum(is.na(sample_sheet$Sample_Group)))
  {
    result <- paste(result,  (" some samples have Sample_Group with NA !"))
  }
  else
  {
    # Exist at least 3 samples per Sample_Group==Reference
    tabsample <- table(sample_sheet$Sample_Group)
    if(tabsample ["Reference"]<2)
    {
      result <- paste(result,  (" Required at least two samples for the Reference Sample Group"))
    }
  }


  if (sum(colnames(methylation_data) %in% sample_sheet$Sample_ID)!=ncol(methylation_data))
  {
    result <- paste(result,  (" The methylation data has not the column's names as expected from the sample sheet in column Sample_ID!"))
  }
  populationGroups <- as.factor(sample_sheet$Sample_Group)
  # expected R as reference, S as Study, C as Control

  # reference population
  populations <-  c("Reference","Control","Case")
  sample_sheet$Sample_Group <- R.utils::toCamelCase(tolower(sample_sheet$Sample_Group), capitalize=TRUE)
  sample_sheet$Sample_Group <- as.factor(sample_sheet$Sample_Group)
  matchedPopulation <- levels(sample_sheet$Sample_Group) %in% populations
  if (is.element(FALSE, matchedPopulation)) {
    result <- paste(result,  " The Sample_Group should contain only: Reference, Control, Case" )
  }

  return(result)
}
