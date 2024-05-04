markers_performace <- function(inference_details, result_folder, pvalue_column="PVALUE_ADJ_ALL_BH",aggr_fun ="mean", ...)
{
  for (id in 1:length(inference_details))
  {
    inference_detail <- inference_details[id,]
    depth_analysis <- inference_detail$depth_analysis
    family_test <- inference_detail$family_test
    transformation <- inference_detail$transformation

    ssEnv <- init_env( result_folder =  result_folder, start_fresh = FALSE, ...)
    # for each marker read the file and extract the metrics and bind
    # them in a single dataframe
    markers <- unique(ssEnv$keys_markers_figures$MARKER)
    model_metrics <- sort(c(ssEnv$model_metrics, pvalue_column))

    # browser()

    for (marker in markers)
    {
      # marker <- "DELTAQ"
      # message(marker)
      # inference_detail <- inference_details[1,]
      file_name <- inference_file_name(inference_detail, marker, ssEnv$result_folderInference,file_extension="csv", prefix = "", suffix="")
      # browser()
      if (!file.exists(file_name))
      {
        log_event("DEBUG: file ", file_name, " is missed !")
        next
      }
      # read the file
      file <- read.csv2(file_name)
      file <- file[file$DEPTH == depth_analysis,]
      # file$COUNT <- 1
      file$KEY <- paste(file$FIGURE,file$AREA,file$SUBAREA, sep="_")
      # bind the metrics
      if (!exists("final"))
        final <- file
      else
        final <- plyr::rbind.fill(final, file)
    }

    if (!exists("final"))
      return(NULL)

    # count the number of rows for each marker, figure
    count <- aggregate(!(final$MARKER==""), by=list(final$MARKER, final$FIGURE), FUN=sum)
    # assign to each row the count of rows matching by marker and figure
    final$COUNT <- count$x[match(paste(final$MARKER, final$FIGURE), paste(count$Group.1, count$Group.2))]

    # remove entries with pvalue > alpha
    final <- final[final[,pvalue_column] < as.numeric(ssEnv$alpha),]

    # remove rows where p_value_colimn has NA
    final <- final[!is.na(final[,pvalue_column]),]

    # count NA in the pvalue_column for each marker, figure
    count_na <- aggregate(final[,pvalue_column], by=list(final$MARKER, final$FIGURE), FUN=function(x) sum(is.na(x)))
    # assign to each row the count of NA in the pvalue_column matching by marker and figure
    final$COUNT_MISSED <- 0
    if(nrow(count_na) != 0)
    {
      count_na <- count_na[count_na$x >= 0,]
      # assign to each row the count of NA in the pvalue_column matching by marker and figure
      final$COUNT_MISSED <- count_na$x[match(paste(final$MARKER, final$FIGURE), paste(count_na$Group.1, count_na$Group.2))]
    }

    final$COUNT_MISSED <- (final$COUNT_MISSED / final$COUNT)

    # remove columns where any NA
    # final <- final[,colSums(is.na(final))==0]

    # count number of rows for each marker and figure
    count <- aggregate(final$COUNT, by=list(final$MARKER, final$FIGURE), FUN=sum)
    # assign to each row the count of rows matching by marker and figure
    final$COUNT_SIGN <- count$x[match(paste(final$MARKER, final$FIGURE), paste(count$Group.1, count$Group.2))]

    final$COUNT_SIGN <- (final$COUNT_SIGN / final$COUNT)


    # check which exists in columns dataframe
    model_metrics <- unique(model_metrics[model_metrics %in% colnames(final)])
    if(depth_analysis > 2)
      model_metrics <- c(model_metrics, c("COUNT","COUNT_MISSED"))
    model_metrics <- unique(sort(model_metrics))

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
    row <- 0
    # remove  NA from figures and from metrics
    figures <- figures[!is.na(figures)]
    model_metrics <- model_metrics[!is.na(model_metrics)]
    for (fig in figures)
    {
      row <- row  + 1
      col <- 0
      # fig <- "HYPO"
      for(metric in model_metrics){
        col <- col + 1
        # metric <- "MAE"
        final_temp <- final[final$FIGURE==fig,]
        # group by MARKER and do mean of values of each group for metric
        final_temp <- final_temp[,colnames(final_temp) %in% (c("MARKER",metric))]
        # browser()
        if(nrow(final_temp) == 0)
          next
        final_temp <- aggregate(final_temp[,metric], by = list(final_temp[,"MARKER"]), FUN = ifelse(metric %in% c("COUNT_MISSED","COUNT_SIGN"),mean, aggr_fun))
        colnames(final_temp) <- c("MARKER", "REBASED")
        bar_colors <- ssEnv$color_palette[1:length(unique(final_temp$MARKER))]

        min_rebased <- min(final_temp$REBASED, na.rm = TRUE)
        max_rebased <- max(final_temp$REBASED, na.rm = TRUE)

        main_title <- ""
        if (row==1)
          main_title <- paste(metric)

        y_title <- ""
        if (col==1)
          y_title <- fig

        # Create the plot
        ggp <- ggplot2::ggplot(final_temp, ggplot2::aes(x = MARKER, y = REBASED, fill = MARKER)) +
          ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
          ggplot2::scale_fill_manual(values = bar_colors) +  # Assign colors to bars
          ggplot2::labs(title = main_title, x = "", y = y_title) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::coord_cartesian(ylim = c(min_rebased, max_rebased))  # Set y-axis to start from min_rebased

        plot_list[[length(plot_list) + 1]] <- ggp
      }
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
    prfx <- paste(transformation,"_", paste0(markers,collapse="_"))
    file_name <- inference_file_name(inference_detail, prfx,path,file_extension="png", prefix = "", suffix=aggr_fun)
    # save the panel
    ggplot2::ggsave(file = file.path(file_name), gge, width = 4960, height = 2480, units = "px")
    rm(final)
  }
}
