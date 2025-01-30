#' Title manhattan_plot_marker_per_probe
#'
#' @param probe_name
#' @param max_sample
#' @param min_sample
#' @param min_signal_probe
#' @param label_font_size
#' @param hyper_color
#' @param hypo_color
#' @param non_outlier_color
#' @param limit_label_color
#' @param limit_line_color
#' @param limit_line_color_median
#' @param reference_samples_color
#' @param result_folder
#' @param maxResources
manhattan_plot_marker_per_probe <- function(probe_name = "cg11680158", max_sample=0, min_sample=0 , min_signal_probe=0, label_font_size=3,
  hyper_color = "blue", hypo_color = "orange",  non_outlier_color = "grey",  limit_label_color = "blue",
  limit_line_color = "red", limit_line_color_median = "black",
  reference_samples_color = "cyan",
  result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)
{


  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE,
    figures="BOTH", areas="PROBE", ...)
  # ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, tech = tech)
  chart_folder <- dir_check_and_create(ssEnv$result_folderChart, "MARKER_PER_PROBE")

  sample_sheet <- read.csv(file_path_build(ssEnv$result_folderData, "sample_sheet_original","csv"), sep=";", header = TRUE, stringsAsFactors = FALSE)
  references <- sample_sheet[sample_sheet$Sample_Group == "Reference", "Sample_ID"]

  annotate_bed()
  create_excel_pivot()
  result_folderPivot <- dir_check_and_create(ssEnv$result_folderData,"Pivots")


  localKeys <- unique(ssEnv$keys_markers_figures$MARKER)
  tempKeys <- localKeys

  high_signal_probes <- 0
  if(probe_name=="")
    for(k in 1:length(localKeys)){
      marker <- as.character(localKeys[k])
      pivot_subfolder <- dir_check_and_create(result_folderPivot,marker)
      fname <- file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv", add_gz = TRUE)
      if (!file.exists(fname))
        tempKeys <- tempKeys[-k]
      else
      {
        if(high_signal_probes==0)
        {
          # get the row woth the max count of values
          count_m <- read.csv(fname, sep=";", skip = 2, row.names = 1)
          # count per each row the number of no zero values
          count_r <- apply(count_m, 1, function(x) length(x[x!=0]))
          if(max(count_r)>high_signal_probes)
          {
            high_signal_probes <- max(count_r)
            probe_name <- names(which.max(count_r))
          }
        }
      }
    }

  probe_name <- toupper(probe_name)
  # remove SIGNAL from tempKeys
  tempKeys <- tempKeys[!grepl("SIGNAL", tempKeys)]
  pivot_subfolder <- dir_check_and_create(result_folderPivot,"SIGNAL")
  # signal_data <- readRDS(file_path_build( ssEnv$result_folderData, "1_signal_data","rds"))
  browser()

  fname <- file_path_build( pivot_subfolder ,c("SIGNAL", "MEAN", "PROBE","PROBE"),"csv", add_gz=TRUE)
  data_name <- read.csv2(fname, sep  =  ";", nrows = 1)
  signal_data <- utils::read.csv2(fname, sep  =  ";", skip = 1)
  # remove column X1
  signal_data <- signal_data[, !grepl("X1", colnames(signal_data))]
  # rename column SAMPLE_GROUP as SAMPLEID
  colnames(signal_data)[colnames(signal_data)=="SAMPLE_GROUP"] <- "SAMPLEID"
  colnames(signal_data) <- colnames(data_name)


  threshold_data <- fst::read_fst(file_path_build( ssEnv$result_folderData,"1_signal_thresholds","fst"))

  colnames(signal_data) <- gsub("-","_", colnames(signal_data))
  if(max_sample!=0)
    signal_data <- signal_data[,min_sample:max_sample]
  samples <- as.vector(colnames(signal_data))
  samples <- samples[-1]

  signal_data$SAMPLEID <- toupper(signal_data$SAMPLEID)
  probe_signal <- as.vector(t(signal_data[signal_data$SAMPLEID==probe_name,]))
  probe_signal <- as.numeric(probe_signal[-1])

  threshold_data$PROBE <- toupper(threshold_data$PROBE)
  probe_threshold <- threshold_data[threshold_data$PROBE==probe_name,]

  q1 <- probe_threshold$q1
  q3 <- probe_threshold$q3
  y_med <- probe_threshold$signal_median_values
  iqr <- probe_threshold$iqr
  y_sup <- probe_threshold$signal_superior_thresholds
  y_inf <- probe_threshold$signal_inferior_thresholds

  probe_signal <- as.numeric(probe_signal)
  outlier <-  as.vector(ifelse(probe_signal > y_sup, "Hyper", ifelse(probe_signal < y_inf, "Hypo" ,"Non-outlier")))
  reference_col_id <- which(colnames(signal_data) %in% references)
  outlier[reference_col_id] <- "Reference"
  marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))


  ########################################################################################################################################################################
  ############################ PLOT MARKERS FOR A GIVEN PROBE ################################################################################################
  ########################################################################################################################################################################

  for(j in 1:length(tempKeys))
  {
    marker <- as.character(tempKeys[j])
    pivot_subfolder <- dir_check_and_create(result_folderPivot,marker)
    fname <- file_path_build( pivot_subfolder ,c(marker, "BOTH", "PROBE","WHOLE"),"csv", add_gz=TRUE)
    message("Reading file ", fname)
    if(!file.exists(fname))
    {
      log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y") ,"File does not exist", fname)
      next
    }
    data_name <- read.csv2(fname, sep  =  ";", nrows = 1)
    marker_data <- utils::read.csv2(fname, sep  =  ";", skip = 1)



    colnames(marker_data) <- colnames(data_name)
    if(any(!(colnames(signal_data) %in% colnames(marker_data))))
    {
      lost_columns <- colnames(signal_data)[!(colnames(signal_data) %in% colnames(marker_data))]
    #  add column with zeros to marker data
    for(l in 1:length(lost_columns))
      marker_data[[lost_columns[l]]] <- 0
    }
    # marker_data <- marker_data[,colnames(marker_data) %in% colnames(signal_data)]
    marker_data <- marker_data[,colnames(signal_data)]
    marker_data$SAMPLEID <- toupper(marker_data$SAMPLEID)
    probe_marker <- marker_data[marker_data$SAMPLEID==probe_name,-1]
    probe_marker <- as.numeric(probe_marker)

    outlier[probe_marker==0] <- "Non-outlier"
    outlier[reference_col_id] <- "Reference"

    data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_marker), outlier)
    marker_bar_end <- 0

    pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, color = outlier)) +
      ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, "Non-outlier" = non_outlier_color,"Reference" = reference_samples_color)) +
      # ggplot2::labs(x = "Sample", y = marker, title = paste0(marker, " found among all samples for probe ", probe_name)) +
      ggplot2::labs(x = "Sample", y = marker, title = "") +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))
      # ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

    plot_filename <- file_path_build(chart_folder,c(probe_name,marker),ssEnv$plot_format)
    ggplot2::ggsave(filename =plot_filename,plot = pp,units = "in", width = 6, height = 2, dpi=as.numeric(ssEnv$plot_resolution_ppi))
  }

  ########################################################################################################################################################################
  ############################ PLOT OUTLIERS FOR A GIVEN PROBE ################################################################################################
  ########################################################################################################################################################################

  #
  data_to_plot <- data.frame(samples = samples, signal_to_plot = as.numeric(probe_signal), outlier)
  delta_label_value <-  ifelse(probe_signal > y_sup, paste("\u03B4 ",
    format(probe_signal - y_sup, scientific=TRUE , digits=2),sep=""),
    ifelse(probe_signal < y_inf,  paste("\u03B4:", format(y_inf - probe_signal,scientific=TRUE, digits=2),sep=""),""))
  delta_label_value[as.numeric(probe_marker)==0] <- ""

  y_shift <- 0
  delta_label_position <-  ifelse(probe_signal > y_sup, (probe_signal + y_shift), ifelse(probe_signal < y_inf, (probe_signal - y_shift) ,0))
  marker_bar_end <-  ifelse(probe_signal > y_sup, y_sup, ifelse(probe_signal < y_inf, y_inf ,probe_signal))
  # marker_segment_color <-  ifelse(probe_signal > y_sup, hyper_color, ifelse(probe_signal < y_inf, hypo_color ,"white"))
  limit_line_color_q <- "red"

  pp <- ggplot2::ggplot(data_to_plot, ggplot2::aes(x = as.numeric(as.factor(samples)), y = signal_to_plot, color = outlier)) +
    ggplot2::geom_segment(ggplot2::aes(xend = as.numeric(factor(samples)), yend = marker_bar_end), color = marker_segment_color, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = y_med, linetype = "dashed", color = limit_line_color_median, size = 0.5) +
    ggplot2::geom_hline(yintercept = y_sup, linetype = "dashed", color = limit_line_color, size = 0.5) +
    ggplot2::geom_hline(yintercept = y_inf, linetype = "dashed", color = limit_line_color, size = 0.5) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = q1, linetype = "dashed", color = limit_line_color_q) +
    ggplot2::geom_hline(yintercept = q3, linetype = "dashed", color = limit_line_color_q) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = q1, label = "Q1"), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = q3, label = "Q3"), color = limit_line_color_q, hjust = -0.2, vjust = -0.3) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = y_med, label = "Median"), color = limit_line_color_median, hjust = -0.2, vjust = -0.3) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = y_sup, label = "Upper Limit"), color = limit_label_color, hjust = -0.2, vjust = 1.5) +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = y_inf, label = "Lower Limit"), color = limit_label_color, hjust = -0.2, vjust = -0.3) +
    ggplot2::geom_text(
      ggplot2::aes( angle = 90, x = as.numeric(factor(samples)), y = delta_label_position, label = delta_label_value),
      color = marker_segment_color,
      check_overlap = TRUE,
      size = label_font_size,
      hjust = +1,
      vjust = -0.45) +
    ggplot2::scale_color_manual(values = c("Hyper" = hyper_color, "Hypo"= hypo_color, "Non-outlier" = non_outlier_color,"Reference" = reference_samples_color)) +
    # ggplot2::labs(x = "Sample", y = "Signal", title = paste0("Signal outliers of all samples for probe ", probe_name)) +
    ggplot2::labs(x = "Sample", y = "Signal", title = "") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), plot.title = ggplot2::element_text(hjust = 0.5))
    # ggplot2::scale_x_continuous(limits = c(1, NA)) # Adjust the limits to exclude 0

  plot_filename <- file_path_build(chart_folder,c(probe_name,"SIGNAL_OUTLIER"),ssEnv$plot_format)
  ggplot2::ggsave(filename = plot_filename,plot = pp,units = "in", width = 6, height = 4, dpi=as.numeric(ssEnv$plot_resolution_ppi))



}


