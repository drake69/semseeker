manhattan_plot_marker_per_probe <- function(probe_name = "cg11680158", max_sample=0, min_sample=0 , min_signal_probe=0, label_font_size=3,
  hyper_color = "grey", hypo_color = "orange",  non_outlier_color = "grey",  limit_label_color = "blue",  limit_line_color = "red", limit_line_color_median = "black",
  result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)
{

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
    figures="BOTH", areas="PROBE", ...)
  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, tech = tech)
  chart_folder <- semseeker:::dir_check_and_create(ssEnv$result_folderChart, "MARKER_PER_PROBE")

  create_excel_pivot()
  result_folderPivot <- semseeker:::dir_check_and_create(ssEnv$result_folderData,"Pivots")


  localKeys <- unique(ssEnv$keys_markers_figures$MARKER)
  tempKeys <- localKeys
  for(k in 1:length(localKeys)){
    marker <- as.character(localKeys[k])
    pivot_subfolder <- semseeker:::dir_check_and_create(result_folderPivot,marker)
    fname <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv")
    if (!file.exists(fname))
      tempKeys <- tempKeys[-k]
  }

  # remove SIGNAL from tempKeys
  tempKeys <- tempKeys[!grepl("SIGNAL", tempKeys)]

  signal_data <- readRDS(semseeker:::file_path_build( ssEnv$result_folderData, "1_signal_data","rds"))
  threshold_data <- fst::read_fst(semseeker:::file_path_build( ssEnv$result_folderData,"1_signal_thresholds","fst"))

  colnames(signal_data) <- gsub("-","_", colnames(signal_data))
  if(max_sample!=0)
    signal_data <- signal_data[,min_sample:max_sample]
  samples <- as.vector(colnames(signal_data))
  samples <- samples[-1]

  probe_signal <- as.vector(t(signal_data[signal_data$SAMPLEID==probe_name,]))
  probe_signal <- probe_signal[-1]
  probe_threshold <- threshold_data[threshold_data$PROBE==probe_name,]

  q1 <- probe_threshold$q1
  q3 <- probe_threshold$q3
  y_med <- probe_threshold$signal_median_values
  iqr <- probe_threshold$iqr
  y_sup <- probe_threshold$signal_superior_thresholds
  y_inf <- probe_threshold$signal_inferior_thresholds

  probe_signal <- as.numeric(probe_signal)
  outlier <-  as.vector(ifelse(probe_signal > y_sup, "Hyper", ifelse(probe_signal < y_inf, "Hypo" ,"Non-outlier")))
  marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))


  ########################################################################################################################################################################
  ############################ PLOT MARKERS FOR A GIVEN PROBE ################################################################################################
  ########################################################################################################################################################################

  for(j in 1:length(tempKeys))
  {
    marker <- as.character(tempKeys[j])
    pivot_subfolder <- semseeker:::dir_check_and_create(result_folderPivot,marker)
    fname <- semseeker:::file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv")
    message("Reading file ", fname)
    marker_data <- utils::read.csv2(fname, sep  =  ";")
    marker_data <- marker_data[-1,colnames(marker_data) %in% colnames(signal_data)]
    probe_marker <- marker_data[marker_data$SAMPLEID==probe_name,]
    probe_marker <- as.numeric(as.vector(t(probe_marker[-1])))

    probe_marker[outlier=="Non-outlier"] <- 0
    data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_marker), outlier)
    marker_bar_end <- 0

    pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, color = outlier)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
      ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, "Non-outlier" = non_outlier_color)) +
      # ggplot2::labs(x = "Sample", y = marker, title = paste0(marker, " found among all samples for probe ", probe_name)) +
      ggplot2::labs(x = "Sample", y = marker, title = "") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

    plot_filename <- semseeker:::file_path_build(chart_folder,c(probe_name,marker),".png")
    ggplot2::ggsave(filename =plot_filename,plot = pp,units = "in", width = 6, height = 2, dpi=300)
  }

  ########################################################################################################################################################################
  ############################ PLOT OUTLIERS FOR A GIVEN PROBE ################################################################################################
  ########################################################################################################################################################################


  data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_signal), outlier)
  delta_label_value <-  ifelse(probe_signal > y_sup, paste("\u03B4 ", round(probe_signal - y_sup,3),sep=""), ifelse(probe_signal < y_inf,  paste("\u03B4:", round(y_inf - probe_signal,3),sep=""),""))
  y_shift <- 0
  delta_label_position <-  ifelse(probe_signal > y_sup, (probe_signal + y_shift), ifelse(probe_signal < y_inf, (probe_signal - y_shift) ,0))
  marker_bar_end <-  ifelse(probe_signal > y_sup, y_sup, ifelse(probe_signal < y_inf, y_inf ,probe_signal))
  # marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))

  pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, color = outlier)) +
    ggplot2::geom_point(alpha = 0.6) +
    # ggplot2::geom_hline(yintercept = q1, linetype = "dashed", color = limit_line_color_q) +
    # ggplot2::geom_hline(yintercept = q3, linetype = "dashed", color = limit_line_color_q) +
    ggplot2::geom_hline(yintercept = y_med, linetype = "dashed", color = limit_line_color_median) +
    ggplot2::geom_hline(yintercept = y_sup, linetype = "dashed", color = limit_line_color) +
    ggplot2::geom_hline(yintercept = y_inf, linetype = "dashed", color = limit_line_color) +
    # ggplot2::geom_text(ggplot2::aes(x = 0, y = q1, label = "Q1"), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
    # ggplot2::geom_text(ggplot2::aes(x = 0, y = q3, label = "Q3"), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
    # ggplot2::geom_text(ggplot2::aes(x = 0, y = y_med, label = "Median"), color = limit_line_color_median, hjust = -0.2, vjust = -0.3) +
    # ggplot2::geom_text(ggplot2::aes(x = 0, y = y_sup, label = "Upper Limit"), color = limit_label_color, hjust = -0.2, vjust = 1.5) +
    # ggplot2::geom_text(ggplot2::aes(x = 0, y = y_inf, label = "Lower Limit"), color = limit_label_color, hjust = -0.2, vjust = -0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes( angle = 90, x = as.numeric(factor(samples)), y = delta_label_position, label = delta_label_value),
      color = marker_segment_color,
      check_overlap = TRUE,
      size = label_font_size,
      hjust = +1,
      vjust = -0.45) +
    ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, non_outlier_color)) +
    # ggplot2::labs(x = "Sample", y = "Signal", title = paste0("Signal outliers of all samples for probe ", probe_name)) +
    ggplot2::labs(x = "Sample", y = "Signal", title = "") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

  plot_filename <- semseeker:::file_path_build(chart_folder,c(probe_name,"SIGNAL_OUTLIER"),".png")
  ggplot2::ggsave(filename = plot_filename,plot = pp,units = "in", width = 6, height = 4, dpi=300)


}


