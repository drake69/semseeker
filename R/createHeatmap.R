#' createChartFromMultipleBedGenericPerRegion load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param populations vectors of populations name used to group populations
#' for differential analysis
#' @param figures vectors of figure to identify lesions as HYPER or HYPO
#' @param anomalies vector of anomalies MUTATIONS or LESIONS
#' @param subGroups other grouping column over main group
#' @param probesPrefix definition of prefix to use probes to load
#' @param groupLabel name of column to group the data set
#' @param subGroupLabel  name of column to sub group the data set
#' @param resultFolder folder as root for bedfiles organized per anomaly
#'
#' @return list of pivot by column identified with groupLabel and by Sample

#'
createHeatmap <-
  function(inputBedDataFrame, anomalies, groupLabel, groupColumnIDs ,resultFolder) {

    chartFolder <- paste(resultFolder, "/Charts/", sep="")
    if (chartFolder != "" && !dir.exists(chartFolder)) {
      dir.create(chartFolder)
    }

    chartFolder <- paste(resultFolder, "/Charts/", groupLabel,"/", sep="")
    if (chartFolder != "" && !dir.exists(chartFolder)) {
      dir.create(chartFolder)
    }

    if (is.null(inputBedDataFrame))
      return()

    colnames(inputBedDataFrame) <- c("MAINGROUP","SAMPLEID", "SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
    inputBedDataFrame$MAINGROUP <- as.factor(inputBedDataFrame$MAINGROUP)
    inputBedDataFrame$SAMPLEID <- as.factor(inputBedDataFrame$SAMPLEID)
    inputBedDataFrame$SUBGROUP <- as.factor(inputBedDataFrame$SUBGROUP)
    inputBedDataFrame$FIGURE <- as.factor(inputBedDataFrame$FIGURE)
    inputBedDataFrame$ANOMALY <- as.factor(inputBedDataFrame$ANOMALY)
    inputBedDataFrame$POPULATION <- as.factor(inputBedDataFrame$POPULATION)

    # inputBedDataFrame <- data.frame(inputBedDataFrame, "KEY" = paste(inputBedDataFrame[, groupColumnID],"_",inputBedDataFrame$FIGURE,sep=""))
    if(length(groupColumnIDs)==2)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame[, groupColumnIDs[2]],"_",inputBedDataFrame$FIGURE,sep=""))
    }
    if(length(groupColumnIDs)==1)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame$FIGURE,sep=""))
    }
    inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)
    inputBedDataFrame <- subset(inputBedDataFrame, POPULATION != "Reference")
    pops <- unique(inputBedDataFrame$POPULATION)
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Case"] <- "Red"
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Control"] <- "Blue"

    for (anomaly in anomalies)
    {

      # tempDataFrame <- subset(inputBedDataFrame, ANOMALY == anomaly)
      # if(dim(tempDataFrame)[1]==0)
      #   next
      # tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID ~ KEY, value.var = "FREQ", sum)
      # mine.long <- tidyr::pivot_longer(data = tempDataFrame,
      #                           cols = -c(1),
      #                           names_to = "KEY",
      #                           values_to = "FREQ")
      # # head(mine.long)
      # fillLable <- paste0("Number of ", anomaly, sep="" )
      # mine.heatmap <- ggplot2::ggplot(data = mine.long,mapping = ggplot2::aes(x = KEY,y = SAMPLEID, fill = FREQ))
      # mine.heatmap <- mine.heatmap +  ggplot2::geom_tile()
      # mine.heatmap <- mine.heatmap + ggplot2::xlab(label = groupLabel)
      # mine.heatmap <- mine.heatmap + ggplot2::ylab(label = "SAMPLE NAME")
      # mine.heatmap <- mine.heatmap +  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
      #
      # # filename = paste( chartFolder,"/",paste( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, ".png",sep="")
      # # grDevices::png(file= filename, width=1200, height=1200)
      # # mine.heatmap
      # # grDevices::dev.off()
      #
      # filename = paste( chartFolder,"/",paste( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, ".png",sep="")
      #
      # try(
      #   ggplot2::ggsave(
      #     filename = filename,
      #     plot = mine.heatmap,
      #     device = NULL,
      #     path = NULL,
      #     scale = 1,
      #     width = NA,
      #     height = NA,
      #     units = c("in", "cm", "mm"),
      #     dpi = 1200,
      #     limitsize = TRUE
      #   )
      # )

      tempDataFrame <- subset(inputBedDataFrame, ANOMALY == anomaly)
      if(dim(tempDataFrame)[1]==0)
        next
      tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + POPULATION ~ KEY, value.var = "FREQ", sum)
      row.names(tempDataFrame) <- tempDataFrame$SAMPLEID

      # col<- colorRampPalette(c("violet","white","blue"))(1024)

      filename = paste( chartFolder,"/",paste( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, "_1.png",sep="")
      grDevices::png(file= filename, width=2000, height=2000)
      stats::heatmap(as.matrix(tempDataFrame[,3:dim(tempDataFrame)[2]]),
                     scale = "column",
                     RowSideColors =as.vector(tempDataFrame$POPULATION),
                     margins = c(25, 25))
      grDevices::dev.off()
    }
}
