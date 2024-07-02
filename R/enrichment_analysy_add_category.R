enrichment_analysy_add_category <- function(source, data)
{

  ssEnv <- get_session_info()
  #TO DO: normalize naming of differente data sources
  if(source=="phenolyzer")
    return(data)

  data$SS_CATEGORY <- NA
  if(any(source %in% c("WebGestalt","Phenolyzer_WebGestalt")))
  {
    #
    data$type <- toupper(data$type)
    data[data$type=='BP',"SS_CATEGORY"] <- "PROCESS"
    data[data$type=='CC',"SS_CATEGORY"] <- "COMPONENT"
    data[data$type=='MF',"SS_CATEGORY"] <- "FUNCTION"
  }

  if(source=="pathfindR")
  {
    #
    data$sources <- toupper(data$sources)
    data[data$source=='GO-BP',"SS_CATEGORY"] <- "PROCESS"
    data[data$source=='GO-MF',"SS_CATEGORY"] <- "FUNCTION"
    data[data$source=='KEGG',"SS_CATEGORY"] <- "PATHWAY"
    data[data$source=='REACTOME',"SS_CATEGORY"] <- "PATHWAY"
  }

  if(any(source %in% c("STRINGdb","Phenolyzer_STRINGdb")))
  {
    #
    data$category <- toupper(data$category)
    data$SS_CATEGORY <- data$category
    data[data$category=='PROCESS',"SS_CATEGORY"] <- "PROCESS"
    data[data$category=='FUNCTION',"SS_CATEGORY"] <- "FUNCTION"
    data[data$category=='KEGG',"SS_CATEGORY"] <- "PATHWAY"
  }

  data$SS_RANK <- 0
  with_category <- data[!is.na(data$SS_CATEGORY),]
  categories <- na.omit(unique(with_category$SS_CATEGORY))

  key_enrichment_format <- ssEnv$key_enrichment_format


  # Normalization functions
  normalize_minimize <- function(x) (max(x) - x) / (max(x) - min(x))
  normalize_maximize <- function(x) (x - min(x)) / (max(x) - min(x))

  # rank per SS_CATEGORY
  data <- data[is.na(data$SS_CATEGORY),]
  for (c in categories)
  {

    # fdr
    column_to_rank <- key_enrichment_format[key_enrichment_format$label==source,"column_of_pvalue"]
    if(max(with_category[,column_to_rank]) == min(with_category[,column_to_rank]))
      with_category$SCORE_FDR <- 1
    else
      with_category$SCORE_FDR <- normalize_minimize(with_category[,column_to_rank])

    column_to_rank <- key_enrichment_format[key_enrichment_format$label==source,"column_of_enrichment"]
    if(max(with_category[,column_to_rank]) == min(with_category[,column_to_rank]))
      with_category$SCORE_ENRICHMENT <- 1
    else
      with_category$SCORE_ENRICHMENT <- normalize_maximize(with_category[,column_to_rank])


    # do a total score
    with_category$SS_SCORE <- with_category$SCORE_FDR + with_category$SCORE_ENRICHMENT
    # orde by SCORE descending
    with_category <- with_category[order(with_category$SS_SCORE, decreasing = TRUE),]

    # remove columns SCORE_FDR and SCORE_ENRICHMENT
    with_category <- with_category[,!(names(with_category) %in% c("SCORE_FDR","SCORE_ENRICHMENT"))]

    end <- sum(with_category$SS_CATEGORY==c, na.rm = TRUE)
    to_rank <- with_category[with_category$SS_CATEGORY==c,]
    to_rank$SS_RANK <- 1:end
    data <- plyr::rbind.fill(data,to_rank)
  }

  data$SS_SCORE <- round(data$SS_SCORE,2)

  return(data)


}
