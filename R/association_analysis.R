#' Association analysis of SEMseeker's results
#'
#' @param inference_details
#' independent variable: deve essere nalla sample sheet passata a semseeker
#' quando lo abbiamo eseguito la prima volta
#' tipo di regressioni: gaussian, poisson, binomial,quantreg_tau_runs(both as number) eg quantreg_0.25_2000
#' tipi di test: wilcoxon, stats::t.test,
#' tipi di correlazioni: pearson, kendall, spearman
#' MUTATIONS_* ~ tcdd_mother + exam_age
#' transformation to be applied to dependent variable
#' (mutations and lesions): scale, log, log2, log10, exp, none, quantile_quantiles(as number) eg quantile_3
#' depth analysis:
#' 1: sample level
#' 2: type level (gene, DMR, cpgisland) (includes 1)
#' 3: genomic area: gene, body, gene tss1550, gene whole, gene tss200,  (includes 1 and 2)
#' filter_p_value report after adjusting saves only significant nominal p-value
#' @param result_folder where semseeker's results are stored, the root folder
#' @param maxResources percentage of max system's resource to use
#' @param parallel_strategy which strategy to use for parallel execution see future vignete: possible values, none, multisession,sequential, multicore, cluster
#' @param ... other options to filter elaborations
#'
#' @importFrom doRNG %dorng%
#' @export
association_analysis <- function(inference_details,result_folder, maxResources = 90, parallel_strategy  = "multisession", ...)
{

  j <- 0
  k <- 0
  z <- 0
  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)

  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
  }

  figures <- ssEnv$keys_figures[,1]
  markers <- ssEnv$keys_markers[,1]
  sample_groups <- c("Reference","Control","Case")

  create_excel_pivot()

  # if(exists("figures")) figures,  if(exists("markers"))  markers,if(exists("areas"))  areas

  inference_details <- unique(inference_details)
  # variables_to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable", "inference_details", "ssEnv", "%dorng%", "k", "iter", "RNGseed", "checkRNGversion",
  #                          "getRNG", "%||%", ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", ".getRNGattribute", "isNumber", "isReal", "isInteger", ".foreachGlobals", "RNGprovider", ".RNGkind_length", "tail",
  #                          "file_path_build", "dir_check_and_create", "apply_stat_model", "doRNGversion", ".getRNG", "hasRNG", "nextRNG", "RNGkind", "setRNG", "RNGstr")

  # foreach::foreach(z  =  1:nrow(inference_details), .export  =  variables_to_export) %dorng%
  for(z in 1:nrow(inference_details))
  {
    inference_detail <- inference_details[z,]

    filter_p_value <- if(!is.null(inference_detail$filter_p_value)) inference_detail$filter_p_value else TRUE

    covariates <- inference_detail$covariates
    covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
    family_test <- inference_detail$family_test
    if( is.null(family_test) || length(family_test)  ==  0)
    {
      message("WARNING: ", Sys.time(), " One test family_test is missed! Skipped.")
    }
    else
    {
      transformation <- inference_detail$transformation
      if(is.null(transformation) || length(transformation)  ==  0)
      {
        transformation <- NULL
      }

      independent_variable <- gsub(" ","", inference_detail$independent_variable)

      if(independent_variable %in% covariates)
      {
        stop("ERROR: ", Sys.time(), " The independent variable is also present as covariate!")
      }

      if( is.null(independent_variable) || length(independent_variable)  ==  0)
      {
        message("WARNING: ", Sys.time(), " One indipendent variable is missed! Skipped.")
      }
      else
      {
        depth_analysis <- inference_detail$depth_analysis
        if( is.null(depth_analysis) || length(depth_analysis)  ==  0)
        {
          depth_analysis <- 1
          message("WARNING: ", Sys.time(), " Missed depth analysis inference forced to 1.")
        }

        study_summary <-   utils::read.csv2(file_path_build( ssEnv$result_folderData, "sample_sheet_result","csv"))

        # transform independent variable as factor
        if(family_test  ==  "binomial" || family_test  ==  "wilcoxon" || family_test  ==  "t.test")
          study_summary[,independent_variable] <- as.factor(study_summary[,independent_variable])

        result_folderPivot <- dir_check_and_create(ssEnv$result_folderData,"Pivots")

        file_result_prefix <- paste(depth_analysis, as.character(independent_variable),sep = "_")

        study_summary <- subset(study_summary, study_summary$Sample_Group !=  "Reference")
        #########################################################################################################
        #########################################################################################################

        if (!(independent_variable %in% colnames(study_summary)))
        {
          message("WARNING: ", Sys.time(), " This indipendent variabile:", independent_variable, " is missed! Skipping")
        }
        else
        {
          if(is.null(covariates) || length(covariates)  ==  0)
          {
            sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable)])
          }
          else
          {
            #tramsform in factor not numeric covariate
            for(cc in 1:length(covariates))
            {
              cname <- covariates[cc]
              if(!is.numeric(stats::na.omit(study_summary[,cname])))
              {
                study_summary[,cname] <- as.factor(study_summary[,cname])
              }
            }
            sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable, covariates)])
          }

          result  =  data.frame (
            "INDIPENDENT.VARIABLE" = "",
            "MARKER"  =  "",
            "FIGURE"  =  "",
            "AREA"  =  "",
            "SUBAREA"  =  "",
            "AREA_OF_TEST"  =  "",
            "PVALUE"  =  "",
            "PVALUEADJ"  =  "",
            "TEST"  =  "",
            "BETA"  =  "",
            "STD.ERROR" = "",
            "AIC"  =  "",
            "RESIDUALS.SUM"  =  "",
            "FAMILY" = "",
            "R.PACK" = "",
            "TRANSFORMATION"  =  "",
            "COVARIATES"  =  "",
            "SHAPIRO.PVALUE"  =  "",
            "BREUSCH-PAGAN.PVALUE" = "",
            "BARTLETT.PVALUE"  =  "",
            "CASE.LABEL" = "",
            "COUNT.CASE","",
            "MEAN.CASE"  = "",
            "SD.CASE" = "",
            "CONTROL.LABEL" = "",
            "COUNT.CONTROL" = "",
            "MEAN.CONTROL" = "",
            "SD.CONTROL" = "",
            "RHO" = "",
            "CI.LOWER" =  "",
            "CI.UPPER" =  "",
            "CI.LOWER.ADJUSTED" =   NA,
            "CI.UPPER.ADJUSTED" =   NA,
            "N.PERMUTATIONS"  =  ""
          )
          result <- result[-1,]

          if(!is.null(covariates) && !length(covariates)  ==  0)
          {
            sample_names <-   sample_names[, c("Sample_ID", independent_variable, covariates)]
            colnames(sample_names) <- c("Sample_ID", independent_variable, covariates)
          }
          else
          {
            sample_names <-   sample_names[, c("Sample_ID", independent_variable)]
            colnames(sample_names) <- c("Sample_ID", independent_variable)
          }
          ######################################################################################################
          # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
          # binomiale che si vuole usare per la regressione logistica


          for (a in length(markers))
          {
            keys <- expand.grid("MARKER" =  markers[a], "FIGURE" =  figures)
            # # aggiungiamo temporaneamente BETA MEAN
            # levels(keys$FIGURE) <- c(levels(keys$FIGURE), "MEAN")
            # keys$FIGURE[keys$MARKER=="BETA"] <- "MEAN"
            keys <- unique(keys)
            cols <- paste0(keys$MARKER,"_",keys$FIGURE, sep = "")
            # temporaneamente filtriamo per le colonne esistenti
            cols <- cols[cols %in% colnames(study_summary)]
            keys$AREA  =  "SAMPLE_GROUP"
            keys$SUBAREA  =  "SAMPLE"
            iters <- length(cols)

            if(!is.null(covariates) && !length(covariates)  ==  0)
              study_summary <- study_summary[, c(independent_variable, covariates, cols) ]

            fileNameResults <- inference_file_name(inference_detail, markers[a])

            # clean keys from already done association
            if(file.exists(fileNameResults))
            {
              old_results <- utils::read.csv2(fileNameResults, header  =  T)
              keys_markers_figures_areas_done <- unlist(apply(unique(old_results [, c("MARKER","FIGURE","AREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
              keys_to_be_done <- unlist(apply(keys[, c("MARKER","FIGURE","AREA")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
              keys <- keys[!( keys_to_be_done %in% keys_markers_figures_areas_done), ]
            }

            if(nrow(keys)>0)
              for(j in seq_along(nrow(keys)))
              {
                g_start <- 2 + length(covariates)
                result_temp <- apply_stat_model(tempDataFrame  =  study_summary[, c(independent_variable, covariates, cols[j])], g_start  =  g_start , family_test  =  family_test, covariates  =  covariates,
                  key  =  keys[j,], transformation  =  transformation, dototal  =  FALSE, session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,  ...)
                result <- rbind(result, result_temp)
              }

            result <- result[order(result$PVALUEADJ),]

            study_summaryToPlot <- study_summary

            if(independent_variable  ==  "Sample_Group")
              study_summaryToPlot$Sample_Group  <- stats::relevel(as.factor(study_summaryToPlot$Sample_Group), "Control")

            if(family_test  ==  "binomial" || family_test  ==  "wilcoxon")
            {

              chartFolder <- dir_check_and_create(ssEnv$result_folderChart,"sample_group_COMPARISON")

              if(sum(grepl(pattern = "MUTATIONS",x  =  ssEnv$keys_markers))>0)
              {
                filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "MUTATIONS"),"png")
                grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
                graphics::par(mfrow = c(1,3))
                graphics::boxplot(MUTATIONS_HYPO~ study_summaryToPlot[,independent_variable],main = "Hypo Mutations", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(MUTATIONS_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of Mutations", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(MUTATIONS_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper Mutations", data  =  study_summaryToPlot, cex = 2)
                grDevices::dev.off()
              }

              if(sum(grepl(pattern = "LESIONS",x  =  ssEnv$keys_markers))>0)
              {
                filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "LESIONS"),"png")
                grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
                graphics::par(mfrow = c(1,3))
                graphics::boxplot(LESIONS_HYPO~study_summaryToPlot[,independent_variable],main = "Hypo Lesions", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(LESIONS_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of Lesions", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(LESIONS_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper Lesions", data  =  study_summaryToPlot, cex = 2)
                grDevices::dev.off()
              }

              if(sum(grepl(pattern = "DELTAS",x  =  ssEnv$keys_markers))>0)
              {
                study_summaryToPlot$DELTAS_HYPO <- as.numeric(study_summaryToPlot$DELTAS_HYPO)
                study_summaryToPlot$DELTAS_BOTH <- as.numeric(study_summaryToPlot$DELTAS_BOTH)
                study_summaryToPlot$DELTAS_HYPER <- as.numeric(study_summaryToPlot$DELTAS_HYPER)
                filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "DELTAS"),"png")
                grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
                graphics::par(mfrow = c(1,3))
                graphics::boxplot(DELTAS_HYPO~study_summaryToPlot[,independent_variable],main = "Hypo DELTAS Mean", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(DELTAS_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of DELTAS", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(DELTAS_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper DELTAS Mean", data  =  study_summaryToPlot, cex = 2)
                grDevices::dev.off()
              }

              if(sum(grepl(pattern = "DELTAQ",x  =  ssEnv$keys_markers))>0)
              {
                filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "DELTAQ"),"png")
                grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
                graphics::par(mfrow = c(1,3))
                graphics::boxplot(DELTAQ_HYPO~study_summaryToPlot[,independent_variable],main = "Hypo DELTAQ Mean", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(DELTAQ_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of DELTAQ", data  =  study_summaryToPlot, cex = 2)
                graphics::boxplot(DELTAQ_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper DELTAQ Mean", data  =  study_summaryToPlot, cex = 2)
                grDevices::dev.off()
              }
            }


            if(depth_analysis >1)
            {

              keys <- ssEnv$keys_markers_figures_areas
              # clean keys from already done association
              if(file.exists(fileNameResults))
              {
                old_results <- utils::read.csv2(fileNameResults, header  =  T)
                keys_markers_figures_areas_done <- unlist(apply(unique(old_results [, c("MARKER","FIGURE","AREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys_to_be_done <- unlist(apply(keys[, c("MARKER","FIGURE","AREA")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys <- keys[!( keys_to_be_done %in% keys_markers_figures_areas_done), ]
              }

              nkeys <- nrow(keys)
              variables_to_export_nested <- c("variables_to_export", "keys", "result_folderPivot", "sample_names", "independent_variable", "covariates",
                "family_test", "transformation", "ssEnv", "depth_analysis","file_path_build", "apply_stat_model")
              if(nrow(keys)>0)
                # result_temp_foreach <- foreach::foreach(k  =  1:nkeys, .combine  =  rbind, .export  =  variables_to_export_nested) %dorng%
                for (k in 1:nkeys)
                {
                  if(exists("tempDataFrame"))
                    rm(list  =  c("tempDataFrame"))
                  key <- keys [k,]
                  pivot_subfolder <- dir_check_and_create(result_folderPivot, key$MARKER)
                  fname <-file_path_build( pivot_subfolder ,c(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA),"csv")
                  if (file.exists(fname))
                  {
                    message("INFO: ", Sys.time(), " Starting to read pivot:", fname,".")
                    tempDataFrame <- utils::read.csv2(fname, sep  =  ";")
                    #assign the area name (eg gene...) to the rows
                    # has SAMPLEID as name but is the genomic area
                    row.names(tempDataFrame) <- tempDataFrame$SAMPLEID
                    if(length(areas_selection)>0)
                      tempDataFrame <- tempDataFrame[ tempDataFrame[,1] %in% areas_selection, ]
                    #removes area name (eg. gene...)
                    tempDataFrame <- tempDataFrame[,-1]
                    if(plyr::empty(tempDataFrame) | nrow(tempDataFrame)==0)
                      next
                    # max_row_count <- ceiling(10^6/ncol(tempDataFrame))
                    # batch_count <- ceiling(nrow(tempDataFrame)/max_row_count)
                    # for(h in 0:batch_count)
                    # {
                    # tempDataFrameBatch <- tempDataFrame[ (1+h*max_row_count):min((h+1)*max_row_count,nrow(tempDataFrame)), ]
                    tempDataFrameBatch <- tempDataFrame
                    tempDataFrameBatch <- t(tempDataFrameBatch)
                    tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                    tempDataFrameBatch$Sample_ID <- rownames(tempDataFrameBatch)
                    message("INFO: ", Sys.time(), " Read pivot:", fname, " with ", ncol(tempDataFrameBatch), " rows.")
                    if(nrow(tempDataFrameBatch)>1)
                    {
                      tempDataFrameBatch <- subset(tempDataFrameBatch, "SAMPLE_GROUP"  !=   "Reference")
                      tempDataFrameBatch <- subset(tempDataFrameBatch, "SAMPLE_GROUP"  !=   0)
                      tempDataFrameBatch <-  merge( x  =   sample_names, y  =   tempDataFrameBatch,  by.x  =  "Sample_ID",  by.y  =  "Sample_ID" , all.x  =  TRUE)
                      tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                      tempDataFrameBatch$SAMPLE_GROUP <- sample_names[, independent_variable]
                      tempDataFrameBatch[is.na(tempDataFrameBatch)] <- 0
                      tempDataFrameBatch[, independent_variable] <- sample_names[, independent_variable]
                      tempDataFrameBatch <- tempDataFrameBatch[, !(names(tempDataFrameBatch) %in% c("SAMPLE_GROUP","Sample_ID"))]
                      cols <- (gsub(" ", "_", colnames(tempDataFrameBatch)))
                      cols <- (gsub("-", "_", cols))
                      cols <- (gsub(":", "_", cols))
                      cols <- (gsub("/", "_", cols))
                      cols <- (gsub("'", "_", cols))
                      tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                      if(length(colnames(tempDataFrameBatch)) !=  length(cols))
                        stop("ERROR: I'm stopping here data t associate are not correct, file a bug!")
                      colnames(tempDataFrameBatch) <- cols
                      g_start <- 2 + length(covariates)
                      result_temp_local_batch <- apply_stat_model(tempDataFrame  =  tempDataFrameBatch, g_start  =  g_start, family_test  =  family_test, covariates  =  covariates,
                        key  =  key, transformation =  transformation, dototal  =  TRUE,
                        session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,  ...)

                      if(!exists("result_temp_foreach"))
                        result_temp_foreach <- result_temp_local_batch
                      else
                        result_temp_foreach <- rbind(result_temp_local_batch, result_temp_foreach)
                      # result_temp_local_batch
                    }
                  }
                }
            }
          }
        }
      }
    }

    if(exists("result_temp_foreach"))
      result <- rbind(result, result_temp_foreach)

    result <- unique(result)
    result[,"PVALUEADJ_ALL_BH"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BH")
    result[,"PVALUEADJ_ALL_BY"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BY")
    result[,"PVALUEADJ_ALL_FDR"] <- stats::p.adjust(result[,"PVALUE"],method  =  "fdr")
    result <- result[order(result$PVALUEADJ),]

    if(filter_p_value)
      result <- subset(result, result$PVALUE < 0.05 | result$PVALUEADJ < 0.05)

    if(nrow(result)>0)
    {
      if(file.exists(fileNameResults))
      {
        old_results <- utils::read.csv2(fileNameResults, header  =  T)
        result <- plyr::rbind.fill(result, old_results)
      }
      result[,"PVALUEADJ_ALL_BH"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BH")
      result[,"PVALUEADJ_ALL_BY"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BY")
      result <- result[order(result$PVALUEADJ),]
      result <- result[,colSums(is.na(result))<nrow(result)]
      utils::write.csv2(result,fileNameResults , row.names  =  FALSE)
    }

  }

  close_env()
}
