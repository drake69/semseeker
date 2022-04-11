#' createHeatmap load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param envir semseekere working infos
#' @param inputBedDataFrame data frame to chart
#' @param anomalies vector of anomalies to manage
#' @param groupLabel main genomic area to char eg: gene
#' @param groupColumnIDs positions of the group coplumn id
#'
#' @return list of pivot by column identified with groupLabel and by Sample
#' @importFrom foreach %dopar%
createHeatmap <-
  function(envir, inputBedDataFrame, anomalies, groupLabel, groupColumnIDs ) {

    # parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c())

    chartFolder <- dir_check_and_create(envir$resultFolderChart,groupLabel)

    if (is.null(inputBedDataFrame))
      return()

    colnames(inputBedDataFrame) <- c("MAINGROUP","SAMPLEID", "SUBGROUP","FREQ","FIGURE","ANOMALY","POPULATION")
    inputBedDataFrame$MAINGROUP <- as.factor(inputBedDataFrame$MAINGROUP)
    inputBedDataFrame$SAMPLEID <- as.factor(inputBedDataFrame$SAMPLEID)
    inputBedDataFrame$SUBGROUP <- as.factor(inputBedDataFrame$SUBGROUP)
    inputBedDataFrame$FIGURE <- as.factor(inputBedDataFrame$FIGURE)
    inputBedDataFrame$ANOMALY <- as.factor(inputBedDataFrame$ANOMALY)
    inputBedDataFrame$POPULATION <- as.factor(inputBedDataFrame$POPULATION)

    # inputBedDataFrame <- data.frame(inputBedDataFrame, "KEY" = paste0(inputBedDataFrame[, groupColumnID],"_",inputBedDataFrame$FIGURE,sep=""))
    if(length(groupColumnIDs)==2)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste0(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame[, groupColumnIDs[2]],"_",inputBedDataFrame$FIGURE,sep=""))
    }
    if(length(groupColumnIDs)==1)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste0(inputBedDataFrame[, groupColumnIDs[1]],"_",inputBedDataFrame$FIGURE,sep=""))
    }
    inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)
    inputBedDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$POPULATION != "Reference")
    pops <- unique(inputBedDataFrame$POPULATION)
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Case"] <- "Cyan"
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Control"] <- "Blue"

    toExport <- c("anomalies", "inputBedDataFrame", "pops", "groupLabel", "chartFolder")
    g <- 0
    foreach::foreach(g = 1:length(anomalies), .export = toExport) %dopar%
    # for(g in 1:length(anomalies))
    # for (anomaly in anomalies)
    {

      anomaly <- anomalies[g]
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
      # # filename = paste0( chartFolder,"/",paste0( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, ".png",sep="")
      # # grDevices::png(file= filename, width=1200, height=1200)
      # # mine.heatmap
      # # grDevices::dev.off()
      #
      # filename = paste0( chartFolder,"/",paste0( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, ".png",sep="")
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

      tempDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$ANOMALY == anomaly)
      if(nrow(tempDataFrame)>0)
      {
        tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + POPULATION ~ KEY, value.var = "FREQ", sum)
        row.names(tempDataFrame) <- tempDataFrame$SAMPLEID

        mainTitle <- paste0( paste0( pops, collapse ="_Vs_")," ", groupLabel," ",anomaly, sep="")
        if(nrow(tempDataFrame)>1000 || ncol(tempDataFrame)>1000)
        {
          #reduce
          temp2 <- as.matrix(tempDataFrame[,3:dim(tempDataFrame)[2]])
          temp <- apply(temp2,2, sum)
          temp1 <- sort(temp, decreasing = T)
          limit <- temp1[1000]
          if(sum(temp1==limit)>1)
            limit <- limit + 1
          tempDataFrame <- data.frame(tempDataFrame[,1:2], temp2[,temp1 > limit])
          rm(temp)
          rm(temp1)
          rm(temp2)
          mainTitle <- paste0( mainTitle," (first 1000)",  sep="")
        }

        # col<- colorRampPalette(c("violet","white","blue"))(1024)

        # skip heatmap if no enough data are available
        tt <- tempDataFrame[,3:ncol(tempDataFrame)]
        if (!is.null(tt))
          if(!nrow(tt) < 2 & !ncol(tt) < 2)
          {
            filename = paste0( chartFolder,"/",paste0( pops, collapse ="_Vs_"),"_", groupLabel,"_",anomaly, ".png",sep="")
            grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = 144)
            stats::heatmap(as.matrix(tempDataFrame[,3:ncol(tempDataFrame)]),
                           # col= colorRampPalette(RColorBrewer::brewer.pal(8, "Blues")),
                           scale = "column",
                           RowSideColors =as.vector(tempDataFrame$POPULATION),
                           margins = c(25, 25),
                           main = mainTitle
            )
            grDevices::dev.off()
          }
      }
    }
    gc()
}
