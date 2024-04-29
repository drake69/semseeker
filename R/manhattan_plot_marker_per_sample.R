manhattan_plot_marker_per_sample <- function( sample_name = "NAME", probes_range = 1000:2000,
  hyper_color = ssEnv$color_palette[1],   hypo_color = "orange",   non_outlier_color = "grey", limit_label_color = ssEnv$color_palette[1],
  result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)
{

  limit_line_color <- ssEnv$color_palette[2]
  limit_line_color_q <- "cyan"
  limit_line_color_median <- "black"

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
    figures=c("BOTH","HYPER","HYPO"), areas="PROBE", ...)
  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
  #   figures=c("BOTH","HYPER","HYPO"), areas="PROBE", markers = "DELTAQ")
  chart_folder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart, "MARKER_PER_SAMPLE")

  create_excel_pivot()
  result_folderPivot <- semseeker:::dir_check_and_create(ssEnv$result_folderData,"Pivots")

  localKeys <- unique(ssEnv$keys_markers_figures$MARKER)
  tempKeys <- localKeys
  for(k in 1:length(localKeys)){
    marker <- as.character(localKeys[k])
    pivot_subfolder <- semseeker:::dir_check_and_create(result_folderPivot,marker)
    fname_both <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv")
    fname_hypo <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "HYPO", "PROBE","WHOLE"),"csv")
    fname_hyper <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "HYPER", "PROBE","WHOLE"),"csv")
    
    if (!file.exists(fname_both) | !file.exists(fname_hypo) | !file.exists(fname_hyper))
      tempKeys <- tempKeys[-k]
  }

  # remove SIGNAL from tempKeys
  tempKeys <- tempKeys[!grepl("SIGNAL", tempKeys)]

  for(j in 1:length(tempKeys))
  {
    
    marker <- as.character(tempKeys[j])
    pivot_subfolder <- semseeker:::dir_check_and_create(result_folderPivot,marker)
    fname_both <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv")
    if(!file.exists(fname_both))
      next
    marker_data_both <- utils::read.csv2(fname_both, sep  =  ";")
    marker_data_both <- marker_data_both[-1,]
    # if sample_name is not in the columns, then skip
    if(!(sample_name %in% colnames(marker_data_both)))
      next
    sample_marker_both <- marker_data_both[,sample_name]
    sample_marker_both <- as.numeric(as.vector(t(sample_marker_both)))
    probes <- marker_data_both$SAMPLEID

    fname_hypo <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "HYPO", "PROBE","WHOLE"),"csv")
    if(!file.exists(fname_hypo))
      next
    marker_data_hypo <- utils::read.csv2(fname_hypo, sep  =  ";")
    marker_data_hypo <- marker_data_hypo[-1,]
    # add to marker data hypo missed rows from both
    missed <- setdiff(marker_data_both$SAMPLEID,marker_data_hypo$SAMPLEID)
    missed <- data.frame(SAMPLEID = missed, matrix(0, nrow = length(missed), ncol = ncol(marker_data_hypo)-1))
    colnames(missed) <- colnames(marker_data_hypo)
    marker_data_hypo <- rbind(marker_data_hypo,missed)
    # add to marker missed column from both
    missed <- setdiff(colnames(marker_data_both),colnames(marker_data_hypo))
    missed <- data.frame(SAMPLEID = probes, matrix(0, nrow = length(probes), ncol = length(missed)))
    marker_data_hypo <- cbind(marker_data_hypo,missed)

    # give the same order as in both
    marker_data_hypo <- marker_data_hypo[match(marker_data_both$SAMPLEID,marker_data_hypo$SAMPLEID),]

    sample_marker_hypo <- marker_data_hypo[,sample_name]
    sample_marker_hypo <- as.numeric(as.vector(t(sample_marker_hypo)))
    sample_marker_hypo <- ifelse(sample_marker_hypo!=0,"Hypo","")

    fname_hyper <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "HYPER", "PROBE","WHOLE"),"csv")
    if(!file.exists(fname_hyper))
      next
    marker_data_hyper <- utils::read.csv2(fname_hyper, sep  =  ";")
    marker_data_hyper <- marker_data_hyper[-1,]
    # add to marker data hypo missed rows from both
    missed <- setdiff(marker_data_both$SAMPLEID,marker_data_hyper$SAMPLEID)
    missed <- data.frame(SAMPLEID = missed, matrix(0, nrow = length(missed), ncol = ncol(marker_data_hyper)-1))
    colnames(missed) <- colnames(marker_data_hyper)
    marker_data_hyper <- rbind(marker_data_hyper,missed)

    # add to marker missed column from both
    missed <- setdiff(colnames(marker_data_both),colnames(marker_data_hyper))
    missed <- data.frame(SAMPLEID = probes, matrix(0, nrow = length(probes), ncol = length(missed)))
    marker_data_hyper <- cbind(marker_data_hyper,missed)

    # give the same order as in both
    marker_data_hyper <- marker_data_hyper[match(marker_data_both$SAMPLEID,marker_data_hyper$SAMPLEID),]

    sample_marker_hyper <- marker_data_hyper[,sample_name]
    sample_marker_hyper <- as.numeric(as.vector(t(sample_marker_hyper)))
    sample_marker_hyper <- ifelse(sample_marker_hyper!=0,"Hyper","")

    outlier <- paste(sample_marker_hyper,sample_marker_hypo, sep ="")
    outlier <- ifelse(outlier=="","Non-outlier",outlier)

    marker_segment_color <-  ifelse(outlier =="Hyper", hyper_color, ifelse(outlier=="Hypo", hypo_color ,"white"))

    sample_marker_both[outlier=="Non-outlier"] <- 0
    data_to_plot <- data.frame(probes = probes, signal_to_plot = as.numeric(sample_marker_both), outlier)
    marker_bar_end <- 0

    data_to_plot <- data_to_plot[probes_range,]
    marker_segment_color <- marker_segment_color[probes_range]
    outlier <- outlier[probes_range]
    probes <- probes[probes_range]

    # Assuming data_to_plot is your data frame and it's already filtered or adjusted
    # to exclude the factor level corresponding to x=0

    pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(probes)), y = signal_to_plot, color = outlier)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(probes)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, non_outlier_color)) +
      ggplot2::labs(x = "Probe", y = marker, title ="") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0


    # pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(probes)), y = signal_to_plot, color = outlier)) +
    #   ggplot2::geom_point(alpha = 0.6) +
    #   ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(probes)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
    #   ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, non_outlier_color)) +
    #   # ggplot2::labs(x = "Probe", y = marker, title = paste0(marker," of all probes for sample ", sample_name)) +
    #   ggplot2::labs(x = "Probe", y = marker, title ="") +
    #   ggplot2::theme_classic() +
    #   ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))

    plot_filename <- semseeker:::file_path_build(chart_folder,c(sample_name,marker),"png")
    ggplot2::ggsave(filename =plot_filename,plot = pp,units = "in", width = 6, height = 2, dpi=300)
  }

}
