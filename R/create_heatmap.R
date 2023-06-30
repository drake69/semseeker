#' create_heatmap load the multiple bed resulting from
#' analysis organized into files and folders per marker and produce a pivot
#'
#' @param inputBedDataFrame data frame to chart
#' @param markers vector of markers to manage
#' @param file_prefix main genomic area to char eg: gene
#' @param groupColumnLabels positions of the group coplumn id
#'
#' @return list of pivot by column identified with file_prefix and by Sample
#' @importFrom doRNG %dorng%
create_heatmap <-
  function( inputBedDataFrame, markers, file_prefix, groupColumnLabels ) {

    ssEnv <- get_session_info()

    chartFolder <- dir_check_and_create(ssEnv$result_folderChart,file_prefix)

    if (is.null(inputBedDataFrame) | plyr::empty(inputBedDataFrame) | nrow(inputBedDataFrame)==0)
      return()

    # MUTATIONS "CHR" "SAMPLEID" "GENE" "ACCESSION" "AREA" "POSITION" "VALUE" "FIGURE" "MARKER" "SAMPLE_GROUP"
    # DELTAS "CHR","VALUE","SAMPLEID","GENE","ACCESSION","AREA","POSITION","FIGURE","MARKER","SAMPLE_GROUP"
    # colnames(inputBedDataFrame) <- c("MAINGROUP","SAMPLEID", "SUBAREA","VALUE","FIGURE","MARKER","SAMPLE_GROUP")

    inputBedDataFrame$SAMPLEID <- as.factor(inputBedDataFrame$SAMPLEID)
    inputBedDataFrame$FIGURE <- as.factor(inputBedDataFrame$FIGURE)
    inputBedDataFrame$MARKER <- as.factor(inputBedDataFrame$MARKER)
    inputBedDataFrame$SAMPLE_GROUP <- as.factor(inputBedDataFrame$SAMPLE_GROUP)
    inputBedDataFrame$VALUE <- as.numeric(inputBedDataFrame$VALUE)

    if(length(groupColumnLabels)==2)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste0(inputBedDataFrame[, groupColumnLabels[1]],"_",inputBedDataFrame[, groupColumnLabels[2]],sep=""))
    }
    if(length(groupColumnLabels)==1)
    {
      inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = inputBedDataFrame[, groupColumnLabels[1]])
    }
    inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)

    inputBedDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$SAMPLE_GROUP != "Reference")
    if(nrow(inputBedDataFrame) < 10)
      return()

    pops <- unique(inputBedDataFrame$SAMPLE_GROUP)
    levels(inputBedDataFrame$SAMPLE_GROUP)[levels(inputBedDataFrame$SAMPLE_GROUP)=="Case"] <- "Cyan"
    levels(inputBedDataFrame$SAMPLE_GROUP)[levels(inputBedDataFrame$SAMPLE_GROUP)=="Control"] <- "Blue"

    figures <- unique(inputBedDataFrame$FIGURE)

    variables_to_export <- c("markers", "inputBedDataFrame", "pops", "file_prefix", "chartFolder", "figures",
                              "%dorng%", "j", "iter", "RNGseed", "checkRNGversion", "getRNG", "%||%", ".getDoParName",
                             "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
                             ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger",
                             "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider", ".RNGkind_length",
                             "tail", "RNGstr", "keys","doRNGseq", "%dopar%", "getDoPar","variables_to_export_nested","get_probes")

    # variables_to_export <- c("markers", "inputBedDataFrame", "pops", "file_prefix", "chartFolder","figures",
    #                          "%dorng%", "j", "iter", "RNGseed", "checkRNGversion", "getRNG", "%||%", ".getDoParName",
    #                          "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion", ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger", "nextRNG",
    #                          ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider", ".RNGkind_length", "tail", "RNGstr","keys")
    variables_to_export_nested <- c(variables_to_export, "markers", "inputBedDataFrame", "pops", "file_prefix", "chartFolder","g","figures","variables_to_export")
    # g <- 0
    i <- 0
    j <- 0

    keys <- expand.grid("markers"= markers,"figures"= figures)
    for (g in 1:length(markers))
    # foreach::foreach(g = 1:length(markers), .export = variables_to_export) %dorng%
      {
        foreach::foreach(j = 1:nrow(keys), .export = variables_to_export_nested) %dorng%
          # for(j in 1:nrow(keys))
          {

            figure <- as.character(keys[j,"figures"])
            marker <- as.character(keys[j, "markers"])
            tempDataFrame <- subset(inputBedDataFrame, inputBedDataFrame$MARKER == marker & inputBedDataFrame$FIGURE == figure)
            tempDataFrame$VALUE <- as.numeric(tempDataFrame$VALUE)
            if(!plyr::empty(tempDataFrame))
              if(nrow(tempDataFrame)>2)
              {
                tempDataFrame <- reshape2::dcast(data = tempDataFrame, SAMPLEID + SAMPLE_GROUP ~ KEY, value.var = "VALUE", sum)
                row.names(tempDataFrame) <- tempDataFrame$SAMPLEID

                mainTitle <- paste0( paste0( pops, collapse ="_Vs_")," ", file_prefix," ",marker, sep="")
                if(ncol(tempDataFrame)>1000)
                {
                  tempDataFrame <- tempDataFrame[,1:min(ncol(tempDataFrame),1005)]
                  mainTitle <- paste0( mainTitle," (first 1000)",  sep="")
                }

                # col<- colorRampPalette(c("violet","white","blue"))(1024)
                # skip heatmap if no enough data are available
                tt <- as.data.frame(tempDataFrame[,3:ncol(tempDataFrame)])
                if (!plyr::empty(tt))
                  if(nrow(tt) > 2 & ncol(tt) > 2)
                  {
                    filename = paste0( chartFolder,"/",paste0( pops, collapse ="_Vs_"),"_", file_prefix,"_",marker,"_",figure, ".png",sep="")
                    grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = 144)
                    stats::heatmap(as.matrix(tempDataFrame[,3:ncol(tempDataFrame)]),
                                   col = grDevices::cm.colors(256),
                                   scale = "column",
                                   RowSideColors =as.vector(tempDataFrame$SAMPLE_GROUP),
                                   margins = c(25, 25),
                                   main = mainTitle
                    )
                    grDevices::dev.off()
                  }
              }
          }
      }
    message("INFO: ", Sys.time(), " Heatmap created." )

  }
