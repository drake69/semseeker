box.plot <- function(dataFrameToPlot, independent_variable,dependent_variable, transformation, family_test)
{
    if (!is.family_dicotomic(family_test))
      return()

    ssEnv <- get_session_info()
    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("COMPARISON"))
    filename  =  file_path_build(chartFolder,toupper(c(as.character(transformation), independent_variable,"Vs", dependent_variable)),"png")
    if(!file.exists(filename))
    {
      grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
      formula <- as.formula(paste(dependent_variable,"~",independent_variable, sep=""))
      dataFrameToPlot[, independent_variable] <- as.factor(dataFrameToPlot[, independent_variable])
      # dataFrameToPlot$Sample_Group  <- stats::relevel(as.factor(dataFrameToPlot$Sample_Group), "Control")
      # Number of boxplots
      num_boxplots <- length(unique(dataFrameToPlot[, independent_variable]))
      # Define the labels for the independent variables
      # labels <- levels(independent_variable)
      # replace underscore with space
      # labels <- gsub("_", " ", labels)
      # Create the boxplot with specified colors and labels
      if (num_boxplots > length(ssEnv$color_palette))
        graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2)
      else
        graphics::boxplot(formula, data= dataFrameToPlot,   cex = 2, col = ssEnv$color_palette[1:num_boxplots])
      grDevices::dev.off()
    }
}
