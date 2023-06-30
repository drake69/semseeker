#' @importFrom doRNG %dorng%
violin_plot_only_significative_areas <- function(fileNameResults, inference_detail, figure, marker,metaarea, subgroup, independent_variable, transformation)
{
  inference_inference_file <- ""
  group <- ""
  s <- ""

  ssEnv <- get_session_info()
  inference_file_name <- inference_inference_file(inference_detail)

  # violin plot only significative areas

  areas <- utils::read.csv2(file.path(ssEnv$result_folderData, "/Pivots/", figure,"/", paste(figure,"_",marker,"_",metaarea,"_",subgroup, ".csv", sep="")))
  results_inference <- utils::read.csv2(file.path(ssEnv$result_folderInference,inference_file_name))
  results_inference <- subset(results_inference, "MARKER"==marker & "FIGURE"==figure & "AREA" == group & "SUBAREA" == subgroup
    & "INDIPENDENT.VARIABLE"==independent_variable)

  #pivot hase SAMPLEID over the genomic area of interest
  areas <- areas[areas$SAMPLEID %in% results_inference[results_inference$PVALUEADJ_ALL_BH<0.05, "AREA_OF_TEST"], ]

  sample_name <- colnames(areas)
  colnames(areas) <- areas[1,]

  areas <- areas[-1,]

  sample_sheet <- utils::read.csv2( file.path(ssEnv$result_folderData,"/sample_sheet_result.csv"), sep=";", dec=",")
  metaareas_f <- foreach::foreach(s = 2: ncol(areas), .combine = rbind) %dorng%
  # for( s in 2: ncol(areas) )
  {
    temp <- areas[,c(1,s)]
    SAMPLE_GROUP <- sample_sheet[sample_sheet$Sample_ID==sample_name[s],independent_variable]
    colnames(temp) <- c("AREA","VALUE")
    temp$SAMPLE_GROUP <- SAMPLE_GROUP
    temp$VALUE <- as.numeric(temp$VALUE)
    temp
    # if(exists("metaareas_f"))
    #   metaareas_f <- rbind(metaareas_f, temp)
    # else
    #   metaareas_f <- temp
  }
  metaareas_f$SAMPLE_GROUP <- sprintf("%02d",metaareas_f[,"SAMPLE_GROUP"])
  if(transformation=="log10")
    metaareas_f$VALUE <- log10(metaareas_f$VALUE)
  if (transformation=="scale")
    metaareas_f$VALUE <- scale(metaareas_f$VALUE)
  if (transformation=="none")
    metaareas_f$VALUE <- metaareas_f$VALUE

  p1 <- ggplot2::ggplot(metaareas_f, ggplot2::aes(x="SAMPLE_GROUP", y="VALUE")) + ggplot2::geom_violin()
  p1 <- p1 +  ggplot2::xlab("Dataset") + ggplot2::ylab("Epimutation Score")
  p1 <- p1 + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  # p1 <- p1 + ggplot2::stat_summary(fun.y=median, geom="point", shape=13, size=2)
  # p1 <- p1 + ylim(0,0.5)
  return (p1)
}
