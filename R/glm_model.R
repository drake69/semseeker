#' Title
#'
#' @param family_test regression model to apply
#' @param tempDataFrame data frame to use for the model
#' @param sig.formula formula to apply the model
#'
glm_model <- function(family_test, tempDataFrame, sig.formula)
{
  result_glm  <- stats::glm( sig.formula, family = as.character(family_test), data = as.data.frame(tempDataFrame))
  res <- data.frame("pvalue" = summary(result_glm)$coeff[-1, 4][1])
  # 
  res$statistic_parameter <- (summary(result_glm)$coeff[-1, 1][1])
  res$aic_value <- (result_glm$aic)
  res$std.error <- summary(result_glm)$coeff[-1, 2][1]

  dep_var <- sig.formula_vars(sig.formula)
  dependent_variable <- dep_var$dependent_variable
  independent_variable <- dep_var$independent_variable

  if(family_test=="gaussian")
  {
    # 
    residuals <-  result_glm$resid
    #calculate shapiro of working residuals
    res$shapiro_pvalue_residuals <- if(length(residuals)>3 & length(unique(residuals))>3) (stats::shapiro.test(residuals)$p.value) else NA
    # calculate residuals homoschedadasticity
    res$residuals_breusch_pagan_pvalue_compl <- 1- lmtest::bptest(result_glm)$p.value

    res <- cbind(res,model_performance(fitted_values = result_glm$fitted.values, expected_values = tempDataFrame[,independent_variable],c(),c()))
    if(plot)
    {
      # 
      ggp <- ggplot2::ggplot(tempDataFrame, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, level = 0.95, color = ssEnv$color_palette_darker[3]) +
        ggplot2::ggtitle("") +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)

      # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (start_a * exp(start_b * x))))
      # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")
    }

  }

  if(family_test=="binomial")
  {
    if(plot)
    {
      # ggp <- box.plot(dataFrameToPlot, independent_variable,burdenValue, transformation, family_test)
      # 
      ggp <- ggplot2::ggplot(train.data, ggplot2::aes_string(x = independent_variable, y = dependent_variable)) +
        ggplot2::geom_point( color = ssEnv$color_palette[1] ) +
        ggplot2::geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, level = 0.95, color = ssEnv$color_palette_darker[3]) +
        ggplot2::ggtitle("") +
        ggplot2::xlab(independent_variable) +
        ggplot2::ylab(dependent_variable)

      # data_to_save <- cbind(train.data, predicted = apply(train.data[,independent_variable],1, function(x) (start_a * exp(start_b * x))))
      # colnames(data_to_save) <- c("Independent_Variable","Dependent_Variable","Predicted")
    }

  }

  if (plot){

    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,c("FITTED_MODEL"))
    filename  =  file_path_build(chartFolder,c(as.character(family_test), independent_variable,"Vs",as.character(transformation), dependent_variable),ssEnv$plot_format)

    ggplot2::ggsave(
      filename,
      plot = ggp,
      scale = 1,
      width = 1240,
      height = 1240,
      units = c("px"),
      dpi = 144
    )

  }

  #  calculate r-squared
  # res$r_squared <- summary(result_glm)$r.squared

  res$r_model <- "stats::glm"
  return (res)
}
