# compare inference associations of differente studies
association_cross_studies_overlaps <- function(inference_detail, studies,pvalue = 0.05, adjust_per_area = F,
  adjust_globally = F,pvalue_column="PVALUE_ADJ_ALL_BH",statistic_parameter, adjustment_method = "BH",
  result_folder, ...)
{

  browser()
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  color_palette <- ssEnv$color_palette
  localKeys <- unique(ssEnv$keys_areas_subareas_markers_figures[,c("MARKER","AREA")])

  aggregated_study_results <- data.frame()
  for (a in 1:nrow(localKeys))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each study in studies
    for (s in 1:nrow(studies))
    {
      #
      # get the inference details for the study
      result_folder_study <- studies[s,"RESULT_FOLDER"]
      ssEnv$result_folderInference <- dir_check_and_create(result_folder_study, "Inference")
      update_session_info(ssEnv)
      temp_res <- get_results_areas_inference(inference_details = inference_detail, marker = MARKER, area= AREA,
         adjust_per_area = adjust_per_area, adjust_globally = adjust_globally, pvalue_column= pvalue_column,
        adjustment_method = adjustment_method)
      if(nrow(temp_res) != 0)
        temp_res$STUDY <- studies[s,"STUDY"]
      aggregated_study_results <- plyr::rbind.fill(aggregated_study_results, temp_res)
    }
  }
  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y")," inference files aggregated!")
  aggregated_study_results <- unique(aggregated_study_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STUDY",statistic_parameter, pvalue_column)])

  #
  # change STUDY column to take int account the direction of the statistic
  # aggregated_study_results$STUDY <- as.character(paste0(aggregated_study_results$STUDY,"_", ifelse((aggregated_study_results[,statistic_parameter] > 0),"INCR","DECR") ))
  ssEnv$result_folderInference <- dir_check_and_create(result_folder, "Inference")
  update_session_info(ssEnv)
  ssEnv <- get_session_info(result_folder)
  # reshape to a table with study as columns, area as rows and values as cell without aggreagtion
  aggregated_study_results_table <- reshape2::dcast(aggregated_study_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ STUDY, value.var = pvalue_column)
  #
  studies_to_comb <- na.omit(unique(aggregated_study_results$STUDY))
  # calculate the mean of the pvalues
  aggregated_study_results_table[, pvalue_column] <- apply(aggregated_study_results_table[,studies_to_comb],1,function(x){
    if(any(is.na(x)))
      NA
    else
      mean(x)
  })
  # rename columns
  studies_to_comb_cols <- paste0(studies_to_comb,"_",pvalue_column)
  colnames(aggregated_study_results_table) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, pvalue_column)

  if(statistic_parameter != "")
  {
    aggregated_study_results_table_statistic_parameter <- reshape2::dcast(aggregated_study_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ STUDY, value.var = statistic_parameter)
    aggregated_study_results_table_statistic_parameter[, statistic_parameter] <- apply(aggregated_study_results_table_statistic_parameter[,studies_to_comb],1,function(x){
      if(any(is.na(x)))
        NA
      else
        mean(x)
    })
    studies_to_comb_cols <- paste0(studies_to_comb,"_",statistic_parameter)
    colnames(aggregated_study_results_table_statistic_parameter) <- c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST",studies_to_comb_cols, statistic_parameter)

    aggregated_study_results_table <- merge(aggregated_study_results_table, aggregated_study_results_table_statistic_parameter, by = c("AREA", "SUBAREA", "MARKER", "FIGURE", "AREA_OF_TEST"))
  }
  markers <- unique(aggregated_study_results_table[, c("MARKER")])
  for (i in 1:length(markers))
  {
    marker <- markers[i]
    filename <- inference_file_name(inference_detail, marker,ssEnv$result_folderInference,file_extension = "csv",suffix = "", prefix = "")
    aggregated_study_results_table_marker <- subset(aggregated_study_results_table, MARKER == marker)
    aggregated_study_results_table_marker$DEPTH <- 3
    write.csv2(aggregated_study_results_table_marker, filename, row.names = F)
  }
  filename <- inference_file_name(inference_detail, paste0(markers, collapse = "_") ,ssEnv$result_folderInference,file_extension = "csv",suffix = "AGGREGATED", prefix = "")
  write.csv2(aggregated_study_results_table, filename, row.names = F)
  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  aggregated files saved!")
  for( j in 2:length(studies_to_comb))
  {
    studies_comb <- combinat::combn(studies_to_comb, j)
    if (j==length(studies_to_comb))
      studies_comb <- data.frame("st"=studies_comb)
    # message(studies_comb)
    for(k in 1: ncol(studies_comb))
    {
      #
      results_inference_comb <- subset(aggregated_study_results, STUDY %in% studies_comb[,k])
      keys <- unique(results_inference_comb[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
      study_count <- length(na.omit(unique(results_inference_comb$STUDY)))

      # Load library
      # library(VennDiagram)
      for (i in 1:nrow(keys))
      {
        # i <- 1
        area_set <- results_inference_comb[
          results_inference_comb$AREA == keys[i, ]$AREA &
            results_inference_comb$SUBAREA == keys[i, ]$SUBAREA &
            results_inference_comb$MARKER == keys[i, ]$MARKER &
            results_inference_comb$FIGURE == keys[i, ]$FIGURE
          , c("AREA_OF_TEST", "STUDY"), ]
        # table(is.na(results_inference_comb))
        area_set <- na.omit(area_set)
        categories <- unique(na.omit(area_set$STUDY))
        # remove _
        categories <- gsub("_", " ", categories)
        if(length(categories)==1)
          next
        SPLIT <- split(area_set$AREA_OF_TEST, area_set$STUDY)
        filename <-
          paste(
            ssEnv$result_folderChart,  "/",
            keys[i, ]$AREA,
            "_",
            keys[i, ]$SUBAREA,
            "_",
            keys[i, ]$MARKER,
            "_",
            keys[i, ]$FIGURE,
            "_venn_diagramm.",
            ssEnv$plot_format,
            sep = ""
          )
        overlaps <- Reduce(intersect, SPLIT)
        if(length(overlaps)>0)
        {
          # Chart
          # Set up the Venn diagram parameters
          # color_palette <- color_palette[length(SPLIT)]
          #
          # Plot Venn diagram
          VennDiagram::venn.diagram(
            x = SPLIT,
            fill = color_palette[1:length(SPLIT)],
            alpha = 0.5,
            category.names = categories,
            cat.col = rep( x="black", length(SPLIT)),
            cat.cex = 1.2,
            cat.fontface = "bold",
            cex = 1.5,
            output = TRUE,
            individuals.in.intersections = TRUE,
            disable.logging = TRUE,
            filename = filename
          )
          log_event("DEBUG: ",format(Sys.time(), "%a %b %d %X %Y"),"  venn diagram completed !")

          filename <- inference_file_name(inference_detail, paste(keys[i, ]$AREA, keys[i, ]$SUBAREA, keys[i, ]$MARKER, keys[i, ]$FIGURE , sep="_"),ssEnv$result_folderInference,file_extension = "csv",suffix = "OVERLAPS", prefix = "")
          #
          overlaps <- data.frame("AREA_OF_TEST" = overlaps)
          write.csv2(overlaps, filename, append = TRUE)
        }
      }
    }
  }

  log_event("INFO: ",format(Sys.time(), "%a %b %d %X %Y"),"  job completed !")
  # remove all files ending with .log
  files <- list.files(ssEnv$result_folderInference, pattern = ".log")
  for (f in 1:length(files))
  {
    file <- files[f]
    file.remove(paste(ssEnv$result_folderInference, "/", file, sep = ""))
  }
}
