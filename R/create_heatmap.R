#' create_heatmap load the multiple bed resulting from
#' analysis organized into files and folders per marker and produce a pivot
#'
#' @return nothing
#' @importFrom doRNG %dorng%
create_heatmap <- function() {

  ssEnv <- get_session_info()
  sample_group_comb <- utils::combn(ssEnv$keys_sample_groups, 2)
  for (g in 1:sample_group_comb)
  {
    localKeys <- reshape::expand.grid.df(ssEnv$keys_areas_subareas_markers_figures, sample_group_comb )
    variables_to_export <- c("markers", "annotatedData", "sample_group_comb", "chartFolder", "figures",
      "%dorng%", "j", "iter", "RNGseed", "checkRNGversion", "getRNG", "%||%", ".getDoParName",
      "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", "doRNGversion",
      ".getRNG", ".getRNGattribute", "hasRNG", "isNumber", "isReal", "isInteger",
      "nextRNG", ".foreachGlobals", "RNGkind", "setRNG", "RNGprovider", ".RNGkind_length",
      "tail", "RNGstr", "localKeys","doRNGseq", "%dopar%", "getDoPar","variables_to_export_nested","get_probes")

    variables_to_export_nested <- c(variables_to_export, "markers", "annotatedData", "sample_group_comb", "chartFolder","g","figures","variables_to_export")
    # g <- 0
    i <- 0
    j <- 0

    foreach::foreach(j = 1:nrow(localKeys), .export = variables_to_export_nested) %dorng%
    # for(j in 1:nrow(localKeys))
    {
      area <- localKeys[j,"AREA"]
      subarea <- localKeys[j,"SUBAREA"]
      marker <- localKeys[j,"MARKER"]
      figure <- localKeys[j,"FIGURE"]
      chartFolder <- dir_check_and_create(ssEnv$result_folderChart, area)

      annotatedData <-  read_annotated_bed(figure,marker,area,subarea)
      annotatedData <- subset(annotatedData, annotatedData$SAMPLE_GROUP %in% sample_group_comb)
      annotatedData <- subset(annotatedData, annotatedData$VALUE != 0 )
      annotatedData$SAMPLEID <- as.factor(annotatedData$SAMPLEID)
      annotatedData$FIGURE <- as.factor(annotatedData$FIGURE)
      annotatedData$MARKER <- as.factor(annotatedData$MARKER)
      annotatedData$SAMPLE_GROUP <- as.factor(annotatedData$SAMPLE_GROUP)
      annotatedData$VALUE <- as.numeric(annotatedData$VALUE)
      annotatedData$SUBAREA <- subarea

      if (is.null(annotatedData) | plyr::empty(annotatedData) | nrow(annotatedData)==0 || nrow(annotatedData) < 10)
        next

      if(!plyr::empty(annotatedData))
        if(nrow(annotatedData)>2)
        {
          annotatedData <- reshape2::dcast(data = annotatedData, SAMPLEID + SAMPLE_GROUP ~ SUBAREA, value.var = "VALUE", sum)
          row.names(annotatedData) <- annotatedData$SAMPLEID

          mainTitle <- paste0( paste0( sample_group_comb, collapse ="_Vs_")," ",marker, sep="")
          if(ncol(annotatedData)>1000)
          {
            annotatedData <- annotatedData[,1:min(ncol(annotatedData),1005)]
            mainTitle <- paste0( mainTitle," (first 1000)",  sep="")
          }

          # col<- colorRampPalette(c("violet","white","blue"))(1024)
          # skip heatmap if no enough data are available
          tt <- as.data.frame(annotatedData[,3:ncol(annotatedData)])
          if (!plyr::empty(tt))
            if(nrow(tt) > 2 & ncol(tt) > 2)
            {
              filename = paste0( chartFolder,"/",paste0( sample_group_comb, collapse ="_Vs_"),"_",marker,"_",figure, ".png",sep="")
              grDevices::png(file= filename, width=2480, height = 2480, pointsize = 15, res = 144)
              stats::heatmap(as.matrix(annotatedData[,3:ncol(annotatedData)]),
                col = grDevices::cm.colors(256),
                scale = "column",
                RowSideColors =as.vector(annotatedData$SAMPLE_GROUP),
                margins = c(25, 25),
                main = mainTitle
              )
              grDevices::dev.off()
            }

      }
  }
    message("INFO: ", Sys.time(), " Heatmap created." )

  }
}
