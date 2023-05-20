#' Title
#'
#' @param anomaly investigated anomaly eg. MUTATIONS, DELTAR, DELTAQ
#' @param figure HYPO, HYPER
#' @param group genomic area (eg. GENE, ISLAND, DMR)
#' @param subgroup sub genomic area (TSS1550), depending on the genomic area
#' @param family fullname of the family used for the association analysis
#' @param adjust_method colnames of the pvalue adjusted to use
#' @param phenotype variable to select from the sample_sheet to use for coloring point
#' @param only_significant_areas
#'
#' @export
#'
manhattan_plot_per_area <- function(anomaly,figure,group,subgroup,family, adjust_method,phenotype, only_significant_areas=FALSE){

  if(only_significant_areas)
  {
    # get inference file

    # select for family, anomaly and figure, adjust_method, group and subgroup
  }

  pivot_data <- pivot_to_long_format(anomaly, figure, group,subgroup)

  group_variable <- pivot_data[,group]

  ggplot2::ggplot(pivot_data, aes(x = factor(group_variable), y = ifelse(VALUE < 0, 0, VALUE), fill = phenotype)) +
    geom_point(shape = 21, position = position_jitter(), alpha = 0.7, color = as.numeric(data_frame[,group])) +
    labs(x = group, y = "VALUE", title = paste("Outlier per ",group)) +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = paste(anomaly, " per ", group), x = group, y = anomaly) +
    theme(plot.title = element_text(hjust = 0.5))

}
