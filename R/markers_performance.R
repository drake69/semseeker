markers_performace <- function(inference_detail, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",aggr_fun ="mean", ...)
{
  if (nrow(inference_detail) > 1)
    stop("Only one inference detail allowed.")

  # inference_detail <- expand.grid("independent_variable"= c("Age"),
  #   "family_test"=c("exp_1"),
  #   "covariates"="",
  #   "transformation"="scale",
  #   "depth_analysis"=1,
  #   "filter_p_value" = FALSE)
  # result_folder <- "~/Documents/Dati_Lavoro/age/ewas_data_hub_base20/"
  # ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE)

  depth_analysis <- unique(inference_detail$depth_analysis)
  family_test <- inference_detail$family_test

  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  # for each marker read the file and extract the metrics and bind
  # them in a single dataframe
  markers <- unique(ssEnv$keys_markers_figures$MARKER)
  model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))


  for (marker in markers)
  {
    # marker <- "DELTAQ"
    # message(marker)
    # inference_detail <- inference_details[1,]
    file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
    # browser()
    if (!file.exists(file_name))
      next
    # read the file
    file <- read.csv2(file_name)
    file <- file[file$DEPTH == depth_analysis,]
    file$COUNT <- 1
    file$KEY <- paste(file$FIGURE,file$AREA,file$SUBAREA, sep="_")
    # bind the metrics
    if (!exists("final"))
      final <- file
    else
      final <- plyr::rbind.fill(final, file)
  }

  if (!exists("final"))
    return(NULL)

  # remove entries with pvalue > alpha
  final <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]

  # remove rows where p_value_colimn has NA
  final <- final[!is.na(final[,pvalue_column]),]

  # remove columns where any NA
  final <- final[,colSums(is.na(final))==0]


  # check which exists in columns dataframe
  model_metrics <- unique(model_metrics[model_metrics %in% colnames(final)])
  if(depth_analysis > 2)
    model_metrics <- c(model_metrics, "COUNT")
  model_metrics <- sort(model_metrics)

  if(any("EFFECT_SIZE_MAGNITUDE" %in% model_metrics)){
    # replace large with 4 medium with 3 small with 2 and neglibible with 1
    final$EFFECT_SIZE_MAGNITUDE <- as.character(final$EFFECT_SIZE_MAGNITUDE)
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="large"] <- 4
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="medium"] <- 3
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="small"] <- 2
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="negligible"] <- 1
    final$EFFECT_SIZE_MAGNITUDE <- as.numeric(final$EFFECT_SIZE_MAGNITUDE)
  }

  # browser()
  row_labels <- c()
  # create a list of figures
  figures <- unique(final$FIGURE)
  # creat empty list to store plots
  plot_list <- list()
  # sort the figures
  figures <- sort(figures)
  # sort the markers
  final$MARKER <- factor(final$MARKER, levels = sort(unique(final$MARKER)))
  for (fig in figures)
  {
    # fig <- "HYPO"
    for(metric in model_metrics){
      # metric <- "MAE"
      final_temp <- final[final$FIGURE==fig,]
      # group by MARKER and do mean of values of each group for metric
      final_temp <- final_temp[,colnames(final_temp) %in% (c("MARKER",metric))]
      if(nrow(final_temp) == 0)
        next
      final_temp <- aggregate(final_temp[,metric], by = list(final_temp[,"MARKER"]), FUN = ifelse(metric=="COUNT",sum, aggr_fun))
      colnames(final_temp) <- c("MARKER", "REBASED")
      bar_colors <- ssEnv$color_palette[1:length(unique(final_temp$MARKER))]

      min_rebased <- min(final_temp$REBASED, na.rm = TRUE)
      max_rebased <- max(final_temp$REBASED, na.rm = TRUE)

      # Create the plot
      ggp <- ggplot2::ggplot(final_temp, ggplot2::aes(x = MARKER, y = REBASED, fill = MARKER)) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
        ggplot2::scale_fill_manual(values = bar_colors) +  # Assign colors to bars
        ggplot2::labs(title = paste(metric, " ", fig), x = "Marker", y = "Value") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::coord_cartesian(ylim = c(min_rebased, max_rebased))  # Set y-axis to start from min_rebased

      plot_list[[length(plot_list) + 1]] <- ggp
    }
    row_labels <- c( row_labels, fig)
  }

  if(length(plot_list) == 0)
    return(NULL)

  # build  a panel from plot list
  gge <- gridExtra::grid.arrange(grobs = lapply(plot_list, ggplot2::ggplotGrob), ncol = length(model_metrics))

  # browser()
  # column_labels <- as.data.frame(model_metrics)
  #
  # # Use grid.arrange to arrange the plots with labels
  # gge <- gridExtra::grid.arrange(
  #   grobs = lapply(plot_list, ggplot2::ggplotGrob),
  #   ncol = length(model_metrics),
  #   top = gridExtra::tableGrob(tibble::tibble(column_labels), theme = tibble::tibble(padding = ggplot2::unit(1, "lines"))),
  #   left = gridExtra::textGrob(row_labels, rot = 90, vjust = 1)
  # )

  path <- dir_check_and_create(ssEnv$result_folderChart, "MARKERS_PERFORMANCE")
  file_name <- inference_file_name(inference_detail,paste0(markers,collapse="_"),path,file_extension="png", prefix = "", suffix=aggr_fun)
  # save the panel
  ggplot2::ggsave(file = file.path(file_name), gge, width = 4960, height = 2480, units = "px")
  rm(final)

}
