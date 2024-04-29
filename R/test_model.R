#' Title
#'
#' @param family_test which family test to apply
#' @param tempDataFrame data frame to use with the test
#' @param sig.formula formula to apply
#' @param burdenValue burden colon name
#' @param independent_variable independent variable for regressor
#'
test_model <- function (family_test, tempDataFrame, sig.formula,burdenValue,independent_variable , transformation, plot )
{

  ssEnv <- get_session_info()
  res <- data.frame(pvalue=1)
  if(family_test=="chisq.test")
  {
    if (plot)
      box.plot(tempDataFrame, independent_variable,burdenValue, transformation, family_test)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_chisq <- suppressWarnings(stats::chisq.test(as.matrix(contingency_table)))
    res$pvalue <- result_chisq$p.value
    res$r_model <- "stats_chisq.test"
    res$statistic_parameter <- result_chisq$statistic
    degrees_of_freedom <- result_chisq$parameter
    effect_size <- sqrt(result_chisq$statistic/nrow(tempDataFrame))
    rea$effect_size <- effect_size
    power_result <- pwr::pwr.chisq.test(w = effect_size, N = nrow(tempDataFrame) , df = degrees_of_freedom, sig.level = ssEnv$alpha, power = )
    res$power <- power_result$power
  }

  if (family_test=="fisher.test")
  {
    if (plot)
      box.plot(tempDataFrame, independent_variable,burdenValue, transformation, family_test)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_fisher <- suppressWarnings(stats::fisher.test(as.matrix(contingency_table)))
    res$pvalue <- result_fisher$p.value
    res$r_model <- "stats_fisher.test"
    res$statistic_parameter <- result_fisher$estimate
  }

  if (family_test=="kruskal.test")
  {
    if (plot)
      box.plot(tempDataFrame, independent_variable,burdenValue, transformation, family_test)

    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    dependent_variable <- round(as.numeric(tempDataFrame[,dependent_variable]),3)
    group <- as.factor(tempDataFrame[,independent_variable])
    result_fisher <- suppressWarnings(stats::kruskal.test(x = dependent_variable, g = group))
    res$pvalue <- result_fisher$p.value
    res$r_model <- "stats_kruskal.test"
    statistic_parameter <- result_fisher$statistic
    res$statistic_parameter <- statistic_parameter
    kw_result <- stats::pairwise.wilcox.test(dependent_variable, group)
    PVALUE_KW <- kw_result$p.value

    pvalue <- 0
    significative <- TRUE
    # for each group combination extract the p-value
    for (i in 1:nrow(PVALUE_KW)) {
      # i <-1
      for (j in 1:ncol(PVALUE_KW)) {
        # j <- 1
        p_value <- PVALUE_KW[i,j][1]
        row <- as.character(rownames(PVALUE_KW)[i])
        col <- as.character(colnames(PVALUE_KW)[j])
        pval_name <- paste0("PVALUE_KW_",as.character(row),"_",as.character(col),sep="")
        significative <- significative & p_value < ssEnv$alpha
        pvalue <- max(pvalue, p_value)
        p_value <- data.frame(p_value)
        colnames(p_value) <- pval_name
        if (exists("res"))
          res <- cbind(res, p_value)
        else
          res <- data.frame(p_value)
      }
    }

    res$significative <- significative
    res$pvalue <- pvalue
    # remove rowname from res
    rownames(res) <- NULL


  }

  if(family_test=="jsd")
  {
    # Sample observations for the first and second sample
    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(round(tempDataFrame[,dep_var[[2]]],3), tempDataFrame[,dep_var[[3]]])
    sample1 <- SPLIT[[1]]
    sample2 <- SPLIT[[2]]

    # Combine the unique elements from both samples to create a common event space
    common_events <- unique(c(sample1, sample2))

    # Create adjusted frequency tables for both samples
    frequency_table1_adjusted <- tabulate(match(sample1, common_events), nbins = length(common_events))
    frequency_table2_adjusted <- tabulate(match(sample2, common_events), nbins = length(common_events))

    # Convert adjusted frequencies to probabilities
    probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
    probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)

    # Calculate the Jensen-Shannon distance
    res$statistic_parameter <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
    res$r_model <- "philentropy.JSD"
  }

  if(family_test=="wilcoxon")
  {
    if (plot)
      box.plot(tempDataFrame, independent_variable,burdenValue, transformation, family_test)

    result_w  <- suppressWarnings(stats::wilcox.test(formula= sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE))
    res$pvalue <- result_w$p.value
    res$r_model <- "stats_wilcox.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])

    # Calculate effect size
    es_res <- effsize::VD.A(SPLIT[[1]], SPLIT[[2]])
    res$effect_size_estimate <- es_res$estimate
    res$effect_size_magnitude <- es_res$magnitude

    # Calculate the statistic parameter
    res$statistic_parameter <- median(SPLIT[[1]]) - median(SPLIT[[2]])
    # Calculate rank-biserial correlation as effect size
    res$rbc <- result_w$statistic / (length(SPLIT[[1]]) * length(SPLIT[[2]]))
    # Calculate power
    power_result = pwr::pwr.t2n.test(d = res$statistic_parameter, n1 = length(SPLIT[[1]]), n2=length(SPLIT[[2]]), sig.level = ssEnv$alpha, power = NULL)
    res$power <- power_result$power
  }

  if(family_test=="t.test")
  {
    if (plot)
      box.plot(tempDataFrame, independent_variable,burdenValue, transformation, family_test)

    result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame))
    res$pvalue <- result_w$p.value
    res$r_model <- "stats_t.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    res$statistic_parameter <- mean(SPLIT[[1]]) - mean(SPLIT[[2]])
    power_result = pwr::pwr.t2n.test(d = statistic_parameter, n1 = length(SPLIT[[1]]), n2=length(SPLIT[[2]]), sig.level = ssEnv$alpha, power = NULL)
    res$power <- power_result$power
  }

  if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
  {
    result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
    pvalue <- result_cor$p.value
    r_model <- "stats_cor.test"
    statistic_parameter <- result_cor$estimate
    power_result <- pwr::pwr.r.test(n = nrow(tempDataFrame) , r = statistic_parameter , sig.level = ssEnv$alpha , power = NULL)
    res$power <- power_result$power
  }

  return (res)

}
