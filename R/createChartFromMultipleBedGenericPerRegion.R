#' createChartFromMultipleBedGenericPerRegion load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param populations vectors of populations name used to group populations
#' for differential analysis
#' @param figures vectors of figure to identify lesions as HYPER or HYPO
#' @param anomalies vector of anomalies MUTATIONS or LESIONS
#' @param subGroups other grouping column over main group
#' @param probesPrefix definition of prefix to use probes to load
#' @param mainGroupLabel name of column to group the data set
#' @param subGroupLabel  name of column to sub group the data set
#' @param resultFolder folder as root for bedfiles organized per anomaly
#'
#' @return list of pivot by column identified with mainGroupLabel and by Sample

#'
createChartFromMultipleBedGenericPerRegion <-
  function(resultFolder, populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel ) {

  HYPO <- NULL
  HYPER <- NULL
  POPULATION <- NULL

  chartFolder <- paste(resultFolder, "/Charts/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel,"/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  chartFolder <- paste(resultFolder, "/Charts/", mainGroupLabel,"/Grouped", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  finalBed <-  annotateBed(  populations,figures ,anomalies,subGroups ,probesPrefix ,mainGroupLabel,subGroupLabel, resultFolder)

  if (is.null(finalBed))
    return()

  finalBed <- data.frame(finalBed, "KEY" = paste(finalBed$POPULATION,"_", finalBed[,mainGroupLabel], sep=""))
  finalBed[,mainGroupLabel] <- as.factor(finalBed[,mainGroupLabel])
  finalBed[,"KEY"] <- as.factor(finalBed[,"KEY"])
  finalBed[,"FIGURE"] <- as.factor(finalBed[,"FIGURE"])
  finalBed[,"POPULATION"] <- as.factor(finalBed[,"POPULATION"])

  numberOfCase <- length(unique(subset(finalBed, POPULATION == "Case" )$SAMPLENAME))
  numberOfControl <- length(unique(subset(finalBed, POPULATION == "Control" )$SAMPLENAME))

  # for (pop in populations)
  {
    pop <- "Reference"
    tempPopData <- subset(finalBed, finalBed[,"POPULATION"] != pop)
    for (grp in unique(tempPopData[,subGroupLabel]))
    {
      temp <- subset(tempPopData, tempPopData[,subGroupLabel]==grp)
      temp <- reshape2::dcast(data = temp, KEY + POPULATION  ~ FIGURE, value.var = "freq", sum)
      temp <- reshape2::dcast(data = temp, HYPO + HYPER  ~ POPULATION, value.var="HYPER", length)

      # temp$Case <- round(temp$Case / numberOfCase)
      # temp$Control <- round(temp$Control / numberOfControl)

      # temp$Case <- (temp$Case / numberOfCase)
      # temp$Control <- (temp$Control / numberOfControl)

      temp$SIZE <- temp$Case + temp$Control
      temp$ALPHA <- temp$Case / temp$Control

      # temp <- log10(temp)
      myplot <- ggplot2::ggplot(temp, ggplot2::aes(HYPO, HYPER))
      myplot <- myplot  + ggplot2::geom_point(ggplot2::aes(alpha = ALPHA, size = SIZE))

      # myplot <- myplot  + ggplot2::coord_cartesian(
      #   xlim = max(temp$HYPO),
      #   ylim = max(temp$HYPER),
      #   expand = TRUE,
      #   default = FALSE,
      #   clip = "on"
      # )
      # leg <- "Each symbol represent the number of locus having X hypometilated probes and Y hypermethylated probes\n" +
      #   " the size of the circle represents the count of locus having these figures \n"
      # " the transparency is the case/control ratio black means only case white means only control"

      myplot <- myplot  + ggplot2::ggtitle (paste("Count Hyper methylated probes \n Vs. Hypo Methylated Probes per ", mainGroupLabel, "  in the region: ", grp ," (", dim(temp)[1],")" ,  sep=""))
      # myplot

      try(
        ggplot2::ggsave(
          filename = paste( chartFolder,"/",paste( unique(tempPopData$POPULATION), collapse ="_Vs_"),"_", mainGroupLabel, "_", grp, ".jpg",sep=""),
          plot = myplot,
          device = NULL,
          path = NULL,
          scale = 1,
          width = NA,
          height = NA,
          units = c("in", "cm", "mm"),
          dpi = 300,
          limitsize = TRUE
        )
      )
    }
  }
  # return(tempResult)
  }


