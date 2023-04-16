violin_plot_only_significative_areas <- function(fileNameResults)
{
  ssEnv <- .pkgglobalenv$ssEnv
  inference_file <- inference_inference_file(inference_detail)

  # violin plot only significative areas
  sample_sheet <- read.csv( file.path(result_folder,"/Data/sample_sheet_result.csv"), sep=";", dec=",")

  metaareas <- read.csv2(file.path(result_folder, "/Data/Pivots/", figure,"/", paste(figure,"_",anomaly,"_",metaarea,"_",subgroup, ".csv", sep="")))
  results_inference <- read.csv2(file.path(result_folder,"/Inference/",inference_file))
  metaareas <- metaareas[metaareas$SAMPLEID %in% results_inference[results_inference$PVALUEADJ_ALL_BH<0.05, "AREA_OF_TEST"], ]

  sample_name <- colnames(metaareas)
  colnames(metaareas) <- metaareas[1,]

  metaareas <- metaareas[-1,]
  if(exists("res"))
    rm(res)
  for( s in 2: ncol(metaareas) )
  {
    temp <- metaareas[,c(1,s)]
    POPULATION <- sample_sheet[sample_sheet$Sample_ID==sample_name[s],independent_variable]
    colnames(temp) <- c("AREA","VALUE")
    temp$POPULATION <- POPULATION
    temp$VALUE <- as.numeric(temp$VALUE)
    if(exists("res"))
      res <- rbind(res, temp)
    else
      res <- temp
  }
  metaareas_f <- res
  metaareas_f$POPULATION <- sprintf("%02d",metaareas_f[,"POPULATION"])
  if(transformation=="log10")
    metaareas_f$VALUE <- log10(metaareas_f$VALUE)
  if (transformation=="scale")
    metaareas_f$VALUE <- scale(metaareas_f$VALUE)
  if (transformation=="none")
    metaareas_f$VALUE <- metaareas_f$VALUE

  p1 <- ggplot2::ggplot(metaareas_f, ggplot2::aes(x=POPULATION, y=VALUE)) + ggplot2::geom_violin()
  p1 <- p1 +  ggplot2::xlab("Dataset") + ggplot2::ylab("Epimutation Score")
  p1 <- p1 + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  # p1 <- p1 + ggplot2::stat_summary(fun.y=median, geom="point", shape=13, size=2)
  # p1 <- p1 + ylim(0,0.5)
  return (p1)

}
