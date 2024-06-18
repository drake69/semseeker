marker_performance_pathway_plot <- function(data, rules,file_prfx,path, disease)
{

  ssEnv <- get_session_info()
  # 
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(ggstance) # for position_dodgev

  column_of_id <- rules["column_of_id"]
  column_of_pvalue <- rules["column_of_pvalue"]
  column_of_description <- rules["column_of_description"]
  column_of_enrichment <- rules["column_of_enrichment"]

  colnames(data)[which(colnames(data)==column_of_id[1,1])] <- "ID"
  colnames(data)[which(colnames(data)==column_of_pvalue[1,1])] <- "P_Value"
  colnames(data)[which(colnames(data)==column_of_description[1,1])] <- "Description"
  colnames(data)[which(colnames(data)==column_of_enrichment[1,1])] <- "Enrichment"

  # Transform highest_p to -log10(highest_p)
  data <- data %>%
    mutate(log_fdr = -log10(P_Value))

  figures <- unique(data$FIGURE)

  for (f in figures)
  {
    data_to_plot <- subset(data, FIGURE == f)
    data_to_plot <- subset(data_to_plot, by_keyword == TRUE)

    if(nrow(data_to_plot) == 0)
      next

    data_to_plot$key <- data_to_plot$MARKER

    # divide in 5 bins the enrichment
    # data_to_plot$Fold_Enrichment_Shape <- cut(data_to_plot$Enrichment, breaks = 5, labels = c("21", "22", "23", "24", "25"))

    # divide in 5 quantile
    # more_than_one <- length(unique(data_to_plot$Enrichment))
    # if (more_than_one > 1)
    #   data_to_plot$Fold_Enrichment_Shape <- cut(data_to_plot$Enrichment, quantile(data_to_plot$Enrichment, probs = seq(0, 1, 0.2)), labels = c("21", "22", "23", "24", "25"))
    # else
    #   data_to_plot$Fold_Enrichment_Shape <- "21"


    # # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key, shape = factor(Fold_Enrichment_Shape))) +
    #   geom_point(aes(size = Enrichment), fill = NA, position = position_dodgev(height = 0.5), stroke = 1.5) +
    #   scale_shape_manual(values = c(21, 22, 23, 24,25)) + # Customize shape values as needed
    #   scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   guides(shape = FALSE) + # Hide the shape legend
    #   theme_minimal() +
    #   theme(
    #     axis.title.x = element_text(size = 12),
    #     axis.title.y = element_text(size = 12),
    #     axis.text.x = element_text(size = 10),
    #     axis.text.y = element_text(size = 10),
    #     plot.title = element_text(size = 14, face = 'bold')
    #   )

    # Define your custom colors
    custom_colors <- ssEnv$color_palette

    library(ggplot2)
    library(ggstance)  # Assuming position_dodgev is from ggstance

    # sort data by SS_RANK ascending
    data_to_plot <- data_to_plot[order(data_to_plot$SS_RANK),]

    # Create the lollipop plot with vertical dodge and shapes based on key
    ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = MARKER)) +
      geom_point(aes(shape = MARKER),
        fill = NA,
        position = position_dodgev(height = 0.5),
        stroke = 1.5) +
      scale_shape_manual(values = c(21, 22, 23, 24, 25)) + # Customize shape values as needed
      scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
      scale_color_manual(values = custom_colors) + # Add custom colors for points
      geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = ssEnv$color_palette[1], size = 1) + # Add vertical line at FDR = 1
      labs(
        title = "",
        x = '-log10(P_Value)',
        y = 'Description'
      ) +
      guides(shape = guide_legend(override.aes = list(size = 5))) + # Show the shape legend
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = 'bold')
      ) +
      xlim(c( min(data_to_plot$log_fdr) - 1 , max(data_to_plot$log_fdr) + 0.5)) # Extend the x-axis range

    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key)) +
    #   geom_point(aes(shape = factor(cut(Enrichment, breaks = 5, labels = FALSE))),
    #     fill = NA,
    #     position = position_dodgev(height = 0.5),
    #     stroke = 1.5) +
    #   scale_shape_manual(values = c(21, 22, 23, 24, 25)) + # Customize shape values as needed
    #   scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   scale_color_manual(values = custom_colors) + # Add custom colors for points
    #   geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = ssEnv$color_palette[1], size = 1) + # Add vertical line at FDR = 1
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   guides(shape = FALSE) + # Hide the shape legend
    #   theme_minimal() +
    #   theme(
    #     axis.title.x = element_text(size = 12),
    #     axis.title.y = element_text(size = 12),
    #     axis.text.x = element_text(size = 10),
    #     axis.text.y = element_text(size = 10),
    #     plot.title = element_text(size = 14, face = 'bold')
    #   )

    # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key)) +
    #   geom_point(aes(shape = factor(cut(Enrichment, breaks = 5, labels = FALSE))),
    #     fill = NA,
    #     position = position_dodgev(height = 0.5),
    #     stroke = 1.5) +
    #   scale_shape_manual(values = c(21, 22, 23, 24, 25)) + # Customize shape values as needed
    #   scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) + # Add vertical line at FDR = 1
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   guides(shape = FALSE) + # Hide the shape legend
    #   theme_minimal() +
    #   theme(
    #     axis.title.x = element_text(size = 12),
    #     axis.title.y = element_text(size = 12),
    #     axis.text.x = element_text(size = 10),
    #     axis.text.y = element_text(size = 10),
    #     plot.title = element_text(size = 14, face = 'bold')
    #   )

    # # Create the lollipop plot with vertical dodge and shapes based on Enrichment
    # ggplot(data_to_plot, aes(x = log_fdr, y = Description, size = Enrichment, color = key,)) +
    #   geom_point(aes(size = Enrichment), fill = NA, position = position_dodgev(height = 0.5), stroke = 1.5) +
    #   scale_shape_manual(values = c(21, 22, 23, 24,25)) + # Customize shape values as needed
    #   scale_size_continuous(range = c(3, 10)) + # Adjust the range for point sizes
    #   labs(
    #     title = "",
    #     x = '-log10(P_Value)',
    #     y = 'Description'
    #   ) +
    #   guides(shape = FALSE) + # Hide the shape legend
    #   theme_minimal() +
    #   theme(
    #     axis.title.x = element_text(size = 12),
    #     axis.title.y = element_text(size = 12),
    #     axis.text.x = element_text(size = 10),
    #     axis.text.y = element_text(size = 10),
    #     plot.title = element_text(size = 14, face = 'bold')
    #   )

    dest_folder <- dir_check_and_create(ssEnv$result_folderChart,"PATHWAYS")
    filename <- paste(dest_folder, "/",file_prfx,"_lollipop_plot_",f,"_",rules["label"],ifelse(disease=="","", paste("_", disease, sep="")), ".png",sep = "")
    ggplot2::ggsave(filename,  width = 16, height = 9, dpi = 300)
  }
}
