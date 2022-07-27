batch_correlation_check <- function(envir) {

  i <- g <- 1
  sample_sheet <- utils::read.csv2(file.path(envir$result_folderData,"sample_sheet_result.csv"))
  localKeys <- expand.grid("FIGURE"=envir$keys_figures_default[,1],"ANOMALY"=c("DELTAS","DELTAQ","MUTATIONS","LESIONS"))

  batch_analysis_folder <-dir_check_and_create(envir$result_folderData,"Batch_Analysis")


  summary_cor <- foreach::foreach(i = 1:nrow(localKeys), .combine = rbind ) %dorng%
    # for(i in 1:nrow(localKeys))
    {
      key <- localKeys[i,]
      populations <- unique(sample_sheet$Sample_Group)
      sub_export <- c("populations", "envir","key")
      total_data_for <- foreach::foreach(g = 1:length(populations), .combine = rbind, .export = sub_export ) %dorng%
        # for(g in 1:length(populations))
        {
          pop <- populations[g]
          tempresult_folderData <-dir_check_and_create(envir$result_folderData,c(as.character(pop) ,paste(as.character(key$ANOMALY),"_",as.character(key$FIGURE),sep="")))
          file_to_read <- file_path_build(tempresult_folderData, c("MULTIPLE", as.character(key$ANOMALY), as.character(key$FIGURE)), "fst")
          temp <- fst::read_fst(file_to_read, as.data.table = T)
          if(key$ANOMALY=="MUTATIONS" | key$ANOMALY=="LESIONS")
            temp$VALUE <- 1
          # if(exists("total_data_for"))
          #   total_data_for <- rbind(total_data_for,temp)
          # else
          #   total_data_for <- temp
          temp
        }

      total_data <- total_data_for
      rm(total_data_for)

      if(plyr::empty(total_data))
        NULL

      total_data$KEY <- paste(total_data$CHR, total_data$START, total_data$END, sep="_")
      total_data$VALUE <- as.numeric(total_data$VALUE)
      # create pivot
      tempDataFrame <- reshape2::dcast(data = total_data, formula = SAMPLEID ~ KEY, value.var = "VALUE",
                                       fun.aggregate = sum, drop = TRUE)


      rownames(tempDataFrame) <- tempDataFrame$SAMPLEID
      tempDataFrame <- t(tempDataFrame[,-1])
      # calculate PCA
      pca_result <- stats::prcomp(tempDataFrame, scale=FALSE)
      pca_summary <- summary(pca_result)
      qq <- data.frame("proportion"=(pca_summary$importance["Cumulative Proportion",]))
      qq$dim <- gsub("PC","Dim.",rownames(qq))

      # get result of variables
      res.var <- factoextra::get_pca_var(pca_result)
      pca_contrib <- as.data.frame(res.var$contrib)
      pca_contrib$Sample_ID <- rownames(pca_contrib)
      pca_contrib <- merge(pca_contrib, sample_sheet[,c("Sample_ID","Batch_ID")], by="Sample_ID")



      result_file <- file_path_build(batch_analysis_folder, c("pca_contrib", as.character(key$ANOMALY), as.character(key$FIGURE)), "csv")
      utils::write.csv2(pca_contrib,result_file,row.names = F)

      # calculate correlation
      pca_contrib <- pca_contrib[,-1]

      cor.result <- Hmisc::rcorr(x = as.matrix(pca_contrib), type = c("pearson","spearman"))
      result_cor <- data.frame("p.value"=cor.result$P[,"Batch_ID"], "rho"=cor.result$r[,"Batch_ID"])
      result_cor <- result_cor[order(result_cor$p.value),]

      # save data checked
      result_cor <- result_cor[-(which(rownames(result_cor)=="Batch_ID")),]
      result_cor$dim <- rownames(result_cor)

      result_cor <- merge(result_cor, qq, by="dim")

      result_file <- file_path_build(batch_analysis_folder, c("batch_cor", as.character(key$ANOMALY), as.character(key$FIGURE)), "csv")
      utils::write.csv2(result_cor,result_file,row.names = F)

      result_cor <- subset(result_cor, result_cor$p.value < 0.05)
      if(!plyr::empty(result_cor))
      {
        result_cor$anomaly <- key$ANOMALY
        result_cor$figure <- key$FIGURE
      }

      # if(exists("summary_cor"))
      #   summary_cor <- rbind(result_cor,summary_cor)
      # else
      #   summary_cor <- result_cor
      result_cor
    }
  result_file <- file_path_build(batch_analysis_folder, c("result","cor"), "csv")
  utils::write.csv2(summary_cor,result_file,row.names = F)

}
