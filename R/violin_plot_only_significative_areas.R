#' @importFrom doRNG %dorng%
violin_plot_only_significative_areas <- function(fileNameResults, inference_detail, figure, anomaly,metaarea, subgroup, independent_variable, transformation)
{
  ssEnv <- .pkgglobalenv$ssEnv
  inference_file_name <- inference_inference_file(inference_detail)

  # violin plot only significative areas

  metaareas <- utils::read.csv2(file.path(ssEnv$result_folderData, "/Pivots/", figure,"/", paste(figure,"_",anomaly,"_",metaarea,"_",subgroup, ".csv", sep="")))
  results_inference <- utils::read.csv2(file.path(ssEnv$result_folderInference,inference_file_name))
  results_inference <- subset(results_inference, "ANOMALY"==anomaly & "FIGURE"==figure & "GROUP" == group & "SUBGROUP" == subgroup
    & "INDIPENDENT.VARIABLE"==independent_variable)

  #pivot hase SAMPLEID over the genomic area of interest
  metaareas <- metaareas[metaareas$SAMPLEID %in% results_inference[results_inference$PVALUEADJ_ALL_BH<0.05, "AREA_OF_TEST"], ]

  sample_name <- colnames(metaareas)
  colnames(metaareas) <- metaareas[1,]

  metaareas <- metaareas[-1,]

  sample_sheet <- utils::read.csv2( file.path(ssEnv$result_folderData,"/sample_sheet_result.csv"), sep=";", dec=",")
  metaareas_f <- foreach::foreach(s = 2: ncol(metaareas), .combine = rbind) %dorng%
  # for( s in 2: ncol(metaareas) )
  {
    temp <- metaareas[,c(1,s)]
    POPULATION <- sample_sheet[sample_sheet$Sample_ID==sample_name[s],independent_variable]
    colnames(temp) <- c("AREA","VALUE")
    temp$POPULATION <- POPULATION
    temp$VALUE <- as.numeric(temp$VALUE)
    temp
    # if(exists("metaareas_f"))
    #   metaareas_f <- rbind(metaareas_f, temp)
    # else
    #   metaareas_f <- temp
  }
  metaareas_f$POPULATION <- sprintf("%02d",metaareas_f[,"POPULATION"])
  if(transformation=="log10")
    metaareas_f$VALUE <- log10(metaareas_f$VALUE)
  if (transformation=="scale")
    metaareas_f$VALUE <- scale(metaareas_f$VALUE)
  if (transformation=="none")
    metaareas_f$VALUE <- metaareas_f$VALUE

  p1 <- ggplot2::ggplot(metaareas_f, ggplot2::aes(x="POPULATION", y="VALUE")) + ggplot2::geom_violin()
  p1 <- p1 +  ggplot2::xlab("Dataset") + ggplot2::ylab("Epimutation Score")
  p1 <- p1 + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  # p1 <- p1 + ggplot2::stat_summary(fun.y=median, geom="point", shape=13, size=2)
  # p1 <- p1 + ylim(0,0.5)
  return (p1)
}
