#' create_heatmap load the multiple bed resulting from
#' analysis organized into files and folders per anomaly and produce a pivot
#'
#' @param envir semseekere working infos
#' @param inputBedDataFrame data frame to chart
#' @param anomalies vector of anomalies to manage
#' @param file_prefix main genomic area to char eg: gene
#' @param groupColumnLabels positions of the group coplumn id
#'
#' @return list of pivot by column identified with file_prefix and by Sample
#' @importFrom doRNG %dorng%
create_heatmap <-
  function(envir, inputBedDataFrame, anomalies, file_prefix, groupColumnLabels ) {

    chartFolder <- dir_check_and_create(envir$result_folderChart,file_prefix)

    if (is.null(inputBedDataFrame))
      return()

    # MUTATIONS "CHR" "SAMPLEID" "GENE" "ACCESSION" "GROUP" "POSITION" "VALUE" "FIGURE" "ANOMALY" "POPULATION"
    # DELTAS "CHR","VALUE","SAMPLEID","GENE","ACCESSION","GROUP","POSITION","FIGURE","ANOMALY","POPULATION"
    # colnames(inputBedDataFrame) <- c("MAINGROUP","SAMPLEID", "SUBGROUP","VALUE","FIGURE","ANOMALY","POPULATION")

    inputBedDataFrame$SAMPLEID <- as.factor(inputBedDataFrame$SAMPLEID)
    inputBedDataFrame$FIGURE <- as.factor(inputBedDataFrame$FIGURE)
    inputBedDataFrame$ANOMALY <- as.factor(inputBedDataFrame$ANOMALY)
    inputBedDataFrame$POPULATION <- as.factor(inputBedDataFrame$POPULATION)

    if(length(groupColumnLabels)==2)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste0(inputBedDataFrame[, groupColumnLabels[1]],"_",inputBedDataFrame[, groupColumnLabels[2]],sep=""))
    }
    if(length(groupColumnLabels)==1)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = inputBedDataFrame[, groupColumnLabels[1]])
    }
    inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)
    inputBedDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$POPULATION != "Reference")
    if(nrow(inputBedDataFrame) < 10)
      return()

    pops <- unique(inputBedDataFrame$POPULATION)
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Case"] <- "Cyan"
    levels(inputBedDataFrame$POPULATION)[levels(inputBedDataFrame$POPULATION)=="Control"] <- "Blue"

    figures <- unique(inputBedDataFrame$FIGURE)

    variables_to_export <- c("anomalies", "inputBedDataFrame", "pops", "file_prefix", "chartFolder","figures")
    g <- 0
    foreach::foreach(g = 1:length(anomalies), .export = variables_to_export) %dorng%
      {
        variables_to_export_nested <- c("anomalies", "inputBedDataFrame", "pops", "file_prefix", "chartFolder","g","figures")
        foreach::foreach(j = 1:length(figures), .export = variables_to_export_nested) %dorng%
          {
            figure <- figures[j]
            anomaly <- anomalies[g]
            tempDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$ANOMALY == anomaly & inputBedDataFrame$FIGURE == figure)
            if(!is.null(tempDataFrame))
              if(nrow(tempDataFrame)>2)
              {
                tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + POPULATION ~ KEY, value.var = "VALUE", sum)
                row.names(tempDataFrame) <- tempDataFrame$SAMPLEID

                mainTitle <- paste0( paste0( pops, collapse ="_Vs_")," ", file_prefix," ",anomaly, sep="")
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
                    filename = paste0( chartFolder,"/",paste0( pops, collapse ="_Vs_"),"_", file_prefix,"_",anomaly,"_",figure, ".png",sep="")
                    grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = 144)
                    stats::heatmap(as.matrix(tempDataFrame[,3:ncol(tempDataFrame)]),
                                   col = cm.colors(256),
                                   scale = "column",
                                   RowSideColors =as.vector(tempDataFrame$POPULATION),
                                   margins = c(25, 25),
                                   main = mainTitle
                    )
                    grDevices::dev.off()
                  }
              }
          }
      }

    gc()
  }
