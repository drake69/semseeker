# compare inference associations of differente studies
compare_inference_association_cross_studies <- function(inference_detail, studies,pvalue = 0.05, adjust_per_area = F,
  adjust_globally = F,pvalue_column="PVALUE_ADJ_ALL_BH", adjustment_method = "BH",
  result_folder, ...)
{

  # browser()
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  color_palette <- ssEnv$color_palette
  localKeys <- ssEnv$keys_areas_subareas_markers_figures

  aggregated_study_results <- data.frame()
  for (a in 1:nrow(localKeys))
  {
    MARKER <- localKeys[a, "MARKER"]
    AREA <- localKeys[a, "AREA"]
    # for each study in studies
    for (s in 1:nrow(studies))
    {
      # get the inference details for the study
      result_folder_study <- studies[s,"RESULT_FOLDER"]
      ssEnv$result_folderInference <- semseeker:::dir_check_and_create(result_folder_study, "Inference")
      update_session_info(ssEnv)
      temp_res <- get_results_areas_inference(inference_details = inference_detail, marker = MARKER, area= AREA,
        pvalue = pvalue, adjust_per_area = adjust_per_area, adjust_globally = adjust_globally, pvalue_column= pvalue_column,adjustment_method = adjustment_method)
      if(nrow(temp_res) != 0)
        temp_res$STUDY <- studies[s,"STUDY"]
      aggregated_study_results <- plyr::rbind.fill(aggregated_study_results, temp_res)
    }
  }

  aggregated_study_results <- unique(aggregated_study_results[,c("AREA","SUBAREA","MARKER","FIGURE","AREA_OF_TEST","STUDY", pvalue_column)])

  # change STUDY column to take int account the direction of the statistic
  aggregated_study_results$STUDY <- as.character(paste0(aggregated_study_results$STUDY,"_", ifelse((aggregated_study_results$STATISTIC_PARAMETER > 0),"INCR","DECR") ))
  # browser()
  ssEnv$result_folderInference <- semseeker:::dir_check_and_create(result_folder, "Inference")
  update_session_info(ssEnv)
  ssEnv <- semseeker:::get_session_info(result_folder)
  # reshape to a teble with study as columns, area as rows and values as cell without aggreagtion
  aggregated_study_results_table <- reshape2::dcast(aggregated_study_results, AREA + SUBAREA + MARKER + FIGURE + AREA_OF_TEST ~ STUDY, value.var = pvalue_column)
  # aggregated_study_results_table <- na.omit(aggregated_study_results_table)
  write.csv2(aggregated_study_results_table, semseeker:::file_path_build(ssEnv$result_folderInference, "aggregated_study_results_table", "csv"), row.names = F)
  studies_to_comb <- na.omit(unique(aggregated_study_results$STUDY))
  for( j in 2:length(studies_to_comb))
  {
    studies_comb <- combinat::combn(studies_to_comb, j)
    if (j==length(studies_to_comb))
      studies_comb <- data.frame("st"=studies_comb)
    # message(studies_comb)
    for(k in 1: ncol(studies_comb))
    {
      results_inference_comb <- subset(aggregated_study_results, STUDY %in% studies_comb[,k])
      keys <- unique(results_inference_comb[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
      study_count <- length(na.omit(unique(results_inference_comb$STUDY)))

      # Load library
      # library(VennDiagram)
      for (i in 1:nrow(keys))
      {
        # i <- 1
        area_set <- results_inference_comb[results_inference_comb$SUBAREA == keys[i, ]$SUBAREA &
            results_inference_comb$AREA == keys[i, ]$AREA &
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
            "_venn_diagramm.png",
            sep = ""
          )
        overlaps <- Reduce(intersect, SPLIT)
        if(length(overlaps)>0)
        {
          # Chart
          # Set up the Venn diagram parameters
          # color_palette <- color_palette[length(SPLIT)]
          # browser()
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

          # # Chart
          # venn.plot <- venn.diagram(
          #   x = SPLIT,
          #   category.names = categories,
          #   filename = filename,
          #   output = TRUE,
          #   individuals.in.intersections = TRUE,
          #   disable.logging = T
          # )
          # filename <-
          #   paste(
          #     ssEnv$result_folderInference,  "/",
          #     keys[i, ]$AREA,
          #     "_",
          #     keys[i, ]$SUBAREA,
          #     "_",
          #     keys[i, ]$MARKER,
          #     "_",
          #     keys[i, ]$FIGURE,
          #     "_venn_diagramm.eps",
          #     sep = ""
          #   )
          # grDevices::setEPS()
          # grDevices::postscript(file = filename, fonts = "serif")
          # grid::grid.draw(venn.plot)
          # grDevices::dev.off()

          # overlaps$STUDY <- area_set$STUDY
          filename <-
            paste(
              ssEnv$result_folderInference,  "/",
              keys[i, ]$AREA,
              "_",
              keys[i, ]$SUBAREA,
              "_",
              keys[i, ]$MARKER,
              "_",
              keys[i, ]$FIGURE,
              ".csv",
              sep = ""
            )
          write.csv2(overlaps, filename)
        }
      }
    }
  }

  # remove all files ending with .log
  files <- list.files(ssEnv$result_folderInference, pattern = ".log")
  for (file in files)
  {
    file.remove(paste(ssEnv$result_folderInference, "/", file, sep = ""))
  }
}
