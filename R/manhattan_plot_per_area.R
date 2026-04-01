#' Title
#'
#' @param marker investigated marker eg. MUTATIONS, DELTAR, DELTAQ
#' @param figure HYPO, HYPER
#' @param area genomic area (eg. GENE, ISLAND, DMR)
#' @param subarea sub genomic area (TSS1550), depending on the genomic area
#' @param family fullname of the family used for the association analysis
#' @param adjust_method colnames of the pvalue adjusted to use
#' @param phenotype variable to select from the sample_sheet to use for coloring point
#' @param only_significant_areas TRUE if filter for pvalue < 0.05
#'
#' @examples
#' \dontrun{
#' manhattan_plot_per_area(
#'   marker        = "DELTARP",
#'   figure        = "HYPO",
#'   area          = "GENE",
#'   subarea       = "WHOLE",
#'   family        = "wilcoxon",
#'   adjust_method = "BH",
#'   phenotype     = "Sample_Group"
#' )
#' }
#' @export
#'
manhattan_plot_per_area <- function(marker,figure,area,subarea,family, adjust_method,phenotype, only_significant_areas=FALSE){

  if(only_significant_areas)
  {
    # get inference file

    # select for family, marker and figure, adjust_method, area and subarea
  }

  pivot_data <- pivot_to_long_format(marker, figure, area,subarea)

  group_variable <- pivot_data[,area]

  ggplot2::ggplot(pivot_data, ggplot2::aes(x = factor(group_variable), y = ifelse(pivot_data$VALUE < 0, 0, pivot_data$VALUE), fill = phenotype)) +
    ggplot2::geom_point(shape = 21, position = ggplot2::position_jitter(), alpha = 0.7, color = as.numeric(pivot_data[,area])) +
    ggplot2::labs(x = area, y = "VALUE", title = paste("Outlier per ",area)) +
    ggplot2::scale_fill_gradient(low = "white", high = ssEnv$color_palette[1]) +
    ggplot2::labs(title = paste(marker, " per ", area), x = area, y = marker) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

}
