#' Title
#'
#' @param family_test which family test to apply
#' @param tempDataFrame data frame to use with the test
#' @param sig.formula formula to apply
#' @param burdenValue burden colon name
#' @param independent_variable independent variable for regressor
#'
test_model <- function (family_test, tempDataFrame, sig.formula,burdenValue,independent_variable )
{

  if(family_test=="chisq.test")
  {
    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_chisq <- suppressWarnings(stats::chisq.test(as.matrix(contingency_table)))
    pvalue <- result_chisq$p.value
    r_model <- "stats_chisq.test"
    statistic_parameter <- result_chisq$statistic
  }

  if (family_test=="fisher.test")
  {
    tempDataFrame <- as.data.frame(tempDataFrame)
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    dependent_variable <- dep_var[[2]]
    independent_variable <- dep_var[[3]] # sample_group
    sample1 <- round(tempDataFrame[,dependent_variable],3)
    sample2 <- tempDataFrame[,independent_variable]
    # create a contingency table
    contingency_table <- table( dependent_variable= sample1, independent_variable= sample2)

    result_fisher <- suppressWarnings(stats::fisher.test(as.matrix(contingency_table)))
    pvalue <- result_fisher$p.value
    r_model <- "stats_fisher.test"
    statistic_parameter <- result_fisher$estimate
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
    statistic_parameter <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
    pvalue <- NA
    r_model <- "philentropy.JSD"
  }

  if(family_test=="wilcoxon")
  {
    result_w  <- suppressWarnings(stats::wilcox.test(formula= sig.formula, data = as.data.frame(tempDataFrame), exact=TRUE))
    pvalue <- result_w$p.value
    r_model <- "stats_wilcox.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    statistic_parameter <- mean(SPLIT[[1]]) - mean(SPLIT[[2]])
  }

  if(family_test=="t.test")
  {
    result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame))
    pvalue <- result_w$p.value
    r_model <- "stats_t.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    statistic_parameter <- mean(SPLIT[[1]]) - mean(SPLIT[[2]])
  }

  if( family_test=="pearson" | family_test=="kendall" | family_test=="spearman")
  {
    result_cor <- stats::cor.test(as.numeric(tempDataFrame[,burdenValue]), as.numeric(tempDataFrame[,independent_variable]), method = as.character(family_test))
    pvalue <- result_cor$p.value
    r_model <- "stats_cor.test"
    statistic_parameter <- result_cor$estimate
  }

  ci.lower <- NA
  ci.upper <- NA
  aic_value <- NA
  residuals <- NA
  shapiro_pvalue <- NA
  std.error <- NA
  n_permutations <- NA



  return (data.frame(ci.lower,ci.upper, pvalue, statistic_parameter,aic_value,residuals,shapiro_pvalue, r_model,std.error,n_permutations))

}
