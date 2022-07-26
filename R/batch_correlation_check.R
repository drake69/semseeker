batch_correlation_check <- function(envir) {

  i <- g <- 1
  sample_sheet <- utils::read.csv2(file.path(envir$result_folderData,"sample_sheet_result.csv"))
  localKeys <- expand.grid("FIGURE"=envir$keys_figures_default[,1],"ANOMALY"=c("DELTAS","DELTAQ") ,"EXT"="bedgraph")
  batch_analysis_folder <-dir_check_and_create(envir$result_folderData,"Batch_Analysis")


  summary_cor <- foreach::foreach(i = 1:nrow(localKeys), .combine = rbind ) %dorng%
    {
      key <- localKeys[i,]
      populations <- unique(sample_sheet$Sample_Group)
      sub_export <- c("populations", "envir","key")
      total_data <- foreach::foreach(g = 1:length(populations), .combine = rbind, .export = sub_export ) %dorng%
        {
          pop <- populations[g]
          tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(pop) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
          file_to_read <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), "fst")
          fst::read_fst(file_to_read, as.data.table = T)
          # total_data <- fst::read_fst(file_to_read, as.data.table = T)
        }
      total_data$KEY <- paste(total_data$CHR, total_data$START, total_data$END, sep="_")
      total_data$VALUE <- as.numeric(total_data$VALUE)
      # create pivot
      tempDataFrame <- reshape2::dcast(data = total_data, formula = SAMPLEID ~ KEY, value.var = "VALUE",
                                       fun.aggregate = sum, drop = TRUE)


      rownames(tempDataFrame) <- tempDataFrame$SAMPLEID
      tempDataFrame <- t(tempDataFrame[,-1])
      # calculate PCA
      pca_result <- stats::prcomp(tempDataFrame, scale=FALSE)

      # get result of variables
      res.var <- factoextra::get_pca_var(pca_result)
      tt <- as.data.frame(res.var$contrib)
      tt$Sample_ID <- rownames(tt)
      res.var <- merge(tt, sample_sheet[,c("Sample_ID","Batch_ID")], by="Sample_ID")

      result_file <- file_path_build(batch_analysis_folder, c("pca_contrib", as.character(key$ANOMALY), as.character(key$FIGURE)), "csv")
      utils::write.csv2(res.var,result_file,row.names = F)

      # calculate correlation
      res.var <- res.var[,-1]

      cor.result <- Hmisc::rcorr(x = as.matrix(res.var), type = c("pearson","spearman"))
      result_cor <- data.frame("p.value"=cor.result$P[,"Batch_ID"], "rho"=cor.result$r[,"Batch_ID"])
      result_cor <- result_cor[order(result_cor$p.value),]

      # save data checked
      result_cor <- result_cor[-(which(rownames(result_cor)=="Batch_ID")),]

      result_file <- file_path_build(batch_analysis_folder, c("batch_cor", as.character(key$ANOMALY), as.character(key$FIGURE)), "csv")
      utils::write.csv2(result_cor,result_file,row.names = F)

      result_cor <- subset(result_cor, result_cor$p.value < 0.05)
      result_cor$anomaly <- key$ANOMALY
      result_cor$figure <- key$FIGURE
      result_cor
    }
  result_file <- file_path_build(batch_analysis_folder, c("result","cor"), "csv")
  utils::write.csv2(summary_cor,result_file,row.names = F)

}
