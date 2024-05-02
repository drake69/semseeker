family_test_performance <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH", ...)
{
  # check different families.test
  families.test <- unique(inference_details$family_test)
  if (length(families.test)!=nrow(inference_details))
    stop("All families.test must be different.")

  # all families must have the same scale
  transformations <- unique(inference_details$transformation)
  if (length(transformations) > 1)
    stop("All families.test must have the same transformation")
  # all families must have the same depth
  depths <- unique(inference_details$depth_analysis)
  if (length(depths) > 1)
    stop("All families.test must have the same depth")

  # ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE)
  ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
  # for each family test, read the file and extract the metrics and bind
  # them in a single dataframe
  browser()
  markers <- unique(ssEnv$keys_markers_figures$MARKER)
  model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))

  for (marker in markers)
  {
    # marker <- "DELTAQ"
    for (i in 1:nrow(inference_details))
    {
      # i <- 1
      inference_detail <- inference_details[i,]
      family_test <- inference_details[i,"family_test"]
      file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
      if(!file.exists(file_name))
        next
      # read the file
      file <- read.csv2(file_name)
      file <- file[file$DEPTH == depths,]
      # file$FAMILY.TEST <- paste(family_test, inference_detail$independent_variable, sep="_")
      file$KEY <- paste(file$FIGURE,file$AREA,file$SUBAREA, sep="_")
      # bind the metrics
      if (!exists("final"))
        final <- file
      else
        final <- plyr::rbind.fill(final, file)
    }
  }

  # remove entries with pvalue > alpha
  final <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]

  # remove rows where p_value_colimn has NA
  final <- final[!is.na(final[,pvalue_column]),]

  # remove columns where any NA
  final <- final[,colSums(is.na(final))==0]



  # check which exists in columns dataframe
  model_metrics <- unique(model_metrics[model_metrics %in% colnames(final)])

  if(any("EFFECT_SIZE_MAGNITUDE" %in% model_metrics)){
    # replace large with 4 medium with 3 small with 2 and neglibible with 1
    final$EFFECT_SIZE_MAGNITUDE <- as.character(final$EFFECT_SIZE_MAGNITUDE)
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="large"] <- 4
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="medium"] <- 3
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="small"] <- 2
    final$EFFECT_SIZE_MAGNITUDE[final$EFFECT_SIZE_MAGNITUDE=="negligible"] <- 1
    final$EFFECT_SIZE_MAGNITUDE <- as.numeric(final$EFFECT_SIZE_MAGNITUDE)
  }

  for (marker in markers)
  {
    if (nrow(final[final$MARKER==marker,]) == 0)
      next
    # marker <- "DELTAQ"
    # creat empty list to store plots
    plot_list <- list()
    # loop over metrics
    for (m in model_metrics)
    {
      # m <- "MAE"
      # create a plot
      p <- ggplot2::ggplot(final[final$MARKER==marker,c("FAMILY.TEST",m)], ggplot2::aes_string(x="FAMILY.TEST", y=m, fill="FAMILY.TEST")) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "none") +
        ggplot2::labs(title = "") +
        ggplot2::ylab(m) +
        ggplot2::xlab("") +
        ggplot2::scale_fill_manual(values = ssEnv$color_palette)
      # store the plot in the list
      plot_list[[m]] <- p
    }

    # build  a panel from plot list
    gge <- gridExtra::grid.arrange(grobs = lapply(plot_list, ggplot2::ggplotGrob), ncol = length(model_metrics))

    path <- dir_check_and_create(ssEnv$result_folderChart, "FAMILY.TEST_PERFORMANCE")
    # file_name <- inference_file_name(inference_detail, marker,path ,file_extension="png", prefix = "", suffix="")
    file_name <- file.path(path, paste0(marker,"_",paste0(families.test, collapse="_"),"_FAMILY_TEST_PERFORMANCE.png"))
    # save the panel
    ggplot2::ggsave(file = file.path(file_name), gge, width = 2048, height = 1024, units = "px")
  }
}
