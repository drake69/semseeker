box.plot <- function(dataFrameToPlot, independent_variable,dependent_variable, transformation, family_test, samples_sql_condition="", area, subarea)
{
    if (!is.family_dicotomic(family_test))
      return()

    ssEnv <- get_session_info()
    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("COMPARISON",name_cleaning(as.character(samples_sql_condition))))
    filename  =  file_path_build(chartFolder,toupper(c(family_test,as.character(transformation), independent_variable,"Vs", dependent_variable,area, subarea)),ssEnv$plot_format)
    if(!file.exists(filename))
    {
      # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = ssEnv$plot_resolution)
      # formula <- as.formula(paste(dependent_variable,"~",independent_variable, sep=""))
      # dataFrameToPlot[, independent_variable] <- as.factor(dataFrameToPlot[, independent_variable])
      # # dataFrameToPlot$Sample_Group  <- stats::relevel(as.factor(dataFrameToPlot$Sample_Group), "Control")
      # # Number of boxplots
      num_boxplots <- length(unique(dataFrameToPlot[, independent_variable]))
      # # Define the labels for the independent variables
      # # labels <- levels(independent_variable)
      # # replace underscore with space
      # # labels <- gsub("_", " ", labels)
      # # Create the boxplot with specified colors and labels
      # if (num_boxplots > length(ssEnv$color_palette))
      #   graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2)
      # else
      #   graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2, col = ssEnv$color_palette[1:num_boxplots])
      # grDevices::dev.off()

      p <- ggplot2::ggplot(dataFrameToPlot, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = !!ggplot2::sym(independent_variable)), outlier.size = 1.5) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::theme(legend.position = "bottom")
      if (num_boxplots <= length(ssEnv$color_palette)) {
        p <- p + ggplot2::scale_fill_manual(values = ssEnv$color_palette[1:num_boxplots])
      }

      ## add pvalue
      if (family_test == "t.test") {
        p <- p + ggpubr::stat_compare_means (method = "t.test", label = "p.format", label.y = 0.5, size = 3)
      } else if (family_test == "wilcox.test") {
        p <- p + ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.5, size = 3)
      } else if (family_test == "kruskal.test") {

        # p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3)
        # Create the plot
        groups <- unique(dataFrameToPlot[,independent_variable])
        pairwise_comparisons <- combn(groups, 2, simplify = FALSE)

        p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.5, size = 3) +
          ggpubr::stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", size = 3)

      } else if (family_test == "anova") {
        p <- p + ggpubr::stat_compare_means(method = "anova", label = "p.format", label.y = 0.5, size = 3)
      }
      # ggplot2::ggsave(filename = filename, plot = p, width = 2480/ssEnv$plot_resolution, height = 2480/ssEnv$plot_resolution, dpi = ssEnv$plot_resolution)
      ggplot2::ggsave(filename = filename, plot = p, dpi = ssEnv$plot_resolution)

    }
}
