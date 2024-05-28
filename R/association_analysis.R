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
#' DEPTH analysis:
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
  arguments <- list(...)
  areas_selection <- c()
  if(!is.null(arguments[["areas_selection"]]))
  {
    areas_selection <-arguments$areas_selection
    arguments[["areas_selection"]] <- NULL
  }

  ssEnv <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, start_fresh = FALSE, ...)

  localKeys <- ssEnv$keys_markers_figures
  sample_groups <- c("Reference","Control","Case")

  create_excel_pivot()

  inference_details <- unique(inference_details)
  # variables_to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable", "inference_details", "ssEnv", "%dorng%", "k", "iter", "RNGseed", "checkRNGversion",
  #                          "getRNG", "%||%", ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", ".getRNGattribute", "isNumber", "isReal", "isInteger", ".foreachGlobals", "RNGprovider", ".RNGkind_length", "tail",
  #                          "file_path_build", "dir_check_and_create", "apply_stat_model", "doRNGversion", ".getRNG", "hasRNG", "nextRNG", "RNGkind", "setRNG", "RNGstr")

  # foreach::foreach(z  =  1:nrow(inference_details), .export  =  variables_to_export) %dorng%
  for(z in 1:nrow(inference_details))
  {
    start_time <- Sys.time()
    processed_items <- 0
    inference_detail <- inference_details[z,]
    collinearity_check <- ifelse(is.null(inference_detail$collinearity_check),FALSE,inference_detail$collinearity_check)
    filter_p_value <- if(!is.null(inference_detail$filter_p_value)) inference_detail$filter_p_value else TRUE

    covariates <- inference_detail$covariates
    covariates <- if(length(covariates) !=  0 && !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split  =  "+", fixed  =  T)))
    # remove covariates with length zero
    covariates <- covariates[lengths(covariates)  !=  0]
    covariates <- covariates[covariates  !=  ""]
    family_test <- inference_detail$family_test
    if (validate_family_test(family_test))
    {
      transformation <- inference_detail$transformation
      if(is.null(transformation) || length(transformation)  ==  0)
        transformation <- NULL

      independent_variable <- gsub(" ","", inference_detail$independent_variable)

      if(independent_variable %in% covariates)
      {
        log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"), " The independent variable is also present as covariate!")
        next
      }

      if( is.null(independent_variable) || length(independent_variable)  ==  0)
      {
        log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " One indipendent variable is missed! Skipped.")
      }
      else
      {
        depth_analysis <- inference_detail$depth_analysis
        if( is.null(depth_analysis) || length(depth_analysis)  ==  0)
        {
          depth_analysis <- 1
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missed DEPTH analysis inference forced to 1.")
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
          log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " This indipendent variabile:", independent_variable, " is missed! Skipping")
        }
        else
        {
          if(is.null(covariates) || length(covariates)  ==  0)
          {
            sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable)])
          }
          else
          {
            #transform in factor not numeric covariate
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

          if(!is.null(covariates) && !length(covariates)  ==  0)
          {
            sample_names <-   sample_names[, c("Sample_ID", independent_variable, covariates)]
            colnames(sample_names) <- c("Sample_ID", independent_variable, covariates)
          }
          else
          {
            sample_names <-   unique(sample_names[, c("Sample_ID", independent_variable)])
            colnames(sample_names) <- c("Sample_ID", independent_variable)
          }
          ######################################################################################################
          # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
          # binomiale che si vuole usare per la regressione logistica
          # # browser()
          # min_covariates_length <- ifelse(grepl("mediation-quantreg", family_test), 2, 1 )
          # min_covariates_length <- ifelse(grepl("mediation-linear", family_test), 2, 1 )
          # min_covariates_length <- ifelse(grepl("mediation-ridge", family_test), length(covariates),1 )
          if(collinearity_check)
            # check collinearity of covariates
          {
            collinearity_score <- calculate_collinearity_score(study_summary[,covariates])
            if(collinearity_score > 0.7)
            {
              log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Collinearity score of covariates is too high!")
              next
            }
          }
          markers <- unique(localKeys$MARKER)
          if (exists("results"))
            rm(results)
          for (a in 1:length(markers) )
          {
            if (exists("results"))
              rm(results)

            keys <- localKeys[localKeys$MARKER==markers[a],]
            keys <- unique(keys)
            cols <- keys$COMBINED
            if (sum(cols %in% colnames(study_summary))!=0)
            {
              # temporaneamente filtriamo per le colonne esistenti
              cols <- cols[cols %in% colnames(study_summary)]
              keys$AREA  =  "SAMPLE_GROUP"
              keys$SUBAREA  =  "SAMPLE"
              g_end <- length(cols)

              if(!is.null(covariates) && !length(covariates)  ==  0)
                study_summary_local <- study_summary[, c(independent_variable, covariates, cols,"Sample_Group") ]
              else
                study_summary_local <- study_summary

              fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference)

              #
              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  > 3
              # clean keys from already done association
              if(file_good)
              {
                old_results <- unique(utils::read.csv2(fileNameResults, header  =  T))
                old_results_filtered <- old_results[old_results$DEPTH == 1,]
                keys_markers_figures_areas_done <- unlist(apply(unique(old_results_filtered [, c("MARKER","FIGURE","AREA","SUBAREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys_to_be_done <- unlist(apply(keys[, c("MARKER","FIGURE","AREA","SUBAREA")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys <- keys[!( keys_to_be_done %in% keys_markers_figures_areas_done), ]
              }

              if(nrow(keys)>0)
                for(j in 1:nrow(keys))
                {
                  key  =  keys[j,]
                  key$FIGURE <- as.character(key$FIGURE)
                  key$MARKER <- as.character(key$MARKER)
                  g_start <- 2 + length(covariates)
                  column_selectors <- c(independent_variable, covariates, key$COMBINED)
                  # remove empty items
                  column_selectors <- column_selectors[column_selectors != ""]
                  processed_items <- processed_items + ncol(study_summary_local) - g_start
                  result_temp <- apply_stat_model(tempDataFrame  =  study_summary_local[, column_selectors], g_start  =  g_start , family_test  =  family_test, covariates  =  covariates,
                    key  =  key, transformation  =  transformation, dototal  =  FALSE, session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,  ...)
                  if (exists("results"))
                    results <- plyr::rbind.fill(results, result_temp)
                  else
                    results <- result_temp

                  # study_summaryToPlot <- study_summary_local
                  # if(independent_variable  ==  "Sample_Group")
                  #   study_summaryToPlot$Sample_Group  <- stats::relevel(as.factor(study_summaryToPlot$Sample_Group), "Control")
                  # chartFolder <- dir_check_and_create(ssEnv$result_folderChart,"SAMPLE_GROUP_COMPARISON")
                  # filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), key$MARKER, key$FIGURE),ssEnv$plot_format)
                  # grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res = 300)
                  # # graphics::par(mfrow = c(1,3))
                  # study_summaryToPlot[,paste(key$MARKER,"_", key$FIGURE, sep="")] <- as.numeric(study_summaryToPlot[,paste(key$MARKER,"_", key$FIGURE, sep="")])
                  # formula <- as.formula(paste(key$MARKER,"_", key$FIGURE,"~",independent_variable, sep=""))
                  # title <- paste(key$MARKER," ", key$FIGURE, sep="")
                  # graphics::boxplot( formula , main = title, data  =  study_summaryToPlot, cex = 2)
                  # grDevices::dev.off()
                }

              if (exists("results"))
                if (!is.null(dim(results)))
                  if (nrow(results)>0)
                    if("PVALUE_ADJ" %in% colnames(results))
                      results <- results[order(results$PVALUE_ADJ),]

              if (exists("old_results"))
              {
                if (exists("results"))
                  results <- plyr::rbind.fill(results, old_results)
                else
                  results <- old_results
                rm(old_results)
              }

              results[results == ""] <- NA
              # remove columns where all rows are NA
              results <- results[, colSums(is.na(results)) < nrow(results)]
              # utils::write.csv2(results, fileNameResults , row.names  =  FALSE)
            }


            #
            # execute for all the areas
            if(depth_analysis>1)
            {
              # browser()

              localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
              keys <- localKeys_1[localKeys_1$MARKER==markers[a],]

              # clean keys from already done association
              fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference)
              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  >10
              dototal <- TRUE
              # clean keys from already done association
              if(file_good)
              {
                old_results <- unique(utils::read.csv2(fileNameResults, header  =  T))
                old_results_filtered <- old_results[(old_results$DEPTH  !=  1),]
                dototal <- !any(old_results_filtered$DEPTH == 2)
                old_results_filtered <- old_results_filtered[(old_results_filtered$DEPTH !=2),]
                keys_markers_figures_areas_done <- unlist(apply(unique(old_results_filtered [, c("MARKER","FIGURE","AREA","SUBAREA")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
                keys_to_be_done <- unlist(apply(keys[, c("MARKER","FIGURE","AREA","SUBAREA")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
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
                  fname <- file_path_build( pivot_subfolder ,c(key$MARKER, key$FIGURE, key$AREA,key$SUBAREA),"csv", add_gz=TRUE)
                  if (file.exists(fname))
                  {
                    areas_selection_temp <- areas_selection
                    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting to read pivot:", fname,".")
                    tempDataFrame <- utils::read.table(gzfile(fname), sep  =  ";", header  =  T)
                    tempDataFrame <- plyr::rename(tempDataFrame, c("SAMPLEID"  =  "Sample_ID"))
                    #removes area name (eg. sample_group name)
                    tempDataFrame <- tempDataFrame[-1,]
                    # # browser()
                    if(length(areas_selection_temp)>0)
                    {
                      #
                      # check if areas_selection_temp is a range
                      if (any(grepl(":",areas_selection_temp)))
                      {
                        areas_selection_temp <- areas_selection_temp[grepl(":",areas_selection_temp)][1]
                        areas_selection_temp <- unlist(strsplit(areas_selection_temp,":"))
                        min_col <- min(as.numeric(areas_selection_temp[1]), nrow(tempDataFrame))
                        max_col <- min(as.numeric(areas_selection_temp[2]), nrow(tempDataFrame))
                        areas_selection_temp <- seq(from  =  min_col, to  = max_col )
                        tempDataFrame <- tempDataFrame[areas_selection_temp, ]
                      }
                      else
                      {
                        # areas_selection_temp <- gsub(":","_",gsub("'","_",gsub("-","_",areas_selection_temp)))
                        # names_area <- gsub(":","_",gsub("'","_",gsub("-","_",tempDataFrame[,1])))
                        tempDataFrame <- tempDataFrame[ tempDataFrame[,1] %in% areas_selection_temp[,1], ]
                      }
                    }
                    if(is.null(dim(tempDataFrame)))
                      next
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
                    colnames(tempDataFrameBatch) <- tempDataFrameBatch[1,]
                    tempDataFrameBatch <- tempDataFrameBatch[-1,]
                    tempDataFrameBatch$Sample_ID <- rownames(tempDataFrameBatch)
                    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Read pivot:", fname, " with ", ncol(tempDataFrameBatch), " rows.")
                    if(nrow(tempDataFrameBatch)>1)
                    {
                      tempDataFrameBatch <-  merge( x  =   sample_names, y  =   tempDataFrameBatch,  by.x  =  "Sample_ID",  by.y  =  "Sample_ID" , all.x  =  TRUE)
                      tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                      tempDataFrameBatch[is.na(tempDataFrameBatch)] <- 0
                      tempDataFrameBatch <- tempDataFrameBatch[,-1]
                      cols <- (gsub(" ", "_", colnames(tempDataFrameBatch)))
                      cols <- (gsub("-", "_", cols))
                      cols <- (gsub(":", "_", cols))
                      cols <- (gsub("/", "_", cols))
                      cols <- (gsub("'", "_", cols))
                      tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                      if(length(colnames(tempDataFrameBatch)) !=  length(cols))
                        stop("ERROR: I'm stopping here data to associate are not correct, file a bug!")
                      colnames(tempDataFrameBatch) <- cols
                      g_start <- 2 + length(covariates)
                      processed_items <- processed_items + ncol(tempDataFrameBatch) - g_start
                      result_temp_local_batch <- apply_stat_model(tempDataFrame  =  tempDataFrameBatch, g_start  =  g_start, family_test  =  family_test, covariates  =  covariates,
                        key  =  key, transformation =  transformation, dototal  =  dototal,
                        session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,  ...)

                      if(!exists("results"))
                        results <- result_temp_local_batch
                      else
                        results <- plyr::rbind.fill(results, result_temp_local_batch)
                    }
                    save_result(results,fileNameResults, family_test, filter_p_value )
                  }
                  if (exists("results"))
                    save_result(results,fileNameResults, family_test, filter_p_value )
                  association_analysis_log(cbind(inference_detail,keys[k,]), start_time, Sys.time(), processed_items)
                }

              if (exists("old_results"))
              {
                if (exists("results"))
                  results <- plyr::rbind.fill(results, old_results)
                else
                  results <- old_results
                rm(old_results)
                save_result(results,fileNameResults, family_test, filter_p_value )
              }
            }
            if(exists("results"))
              save_result(results,fileNameResults, family_test, filter_p_value )

            if(exists("result_temp_foreach"))
              rm(result_temp_foreach)
            if(exists("result_temp_local_batch"))
              rm(result_temp_local_batch)
            if(exists("tempDataFrame"))
              rm(tempDataFrame)

          }
        }
      }
    }
    total_time <- difftime(Sys.time(), start_time, units  =  "mins")
    end_time <- Sys.time()
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Finished processing association analysis for ", processed_items, " items in ", total_time, " minutes." )
  }
  close_env()
}

save_result <- function(results=NULL,fileNameResults, family_test, filter_p_value ){

  #
  if(!exists("results") | is.null(results))
    return

  results <- unique(results)

  pvalue_columns <- colnames(results)[grepl("PVALUE", colnames(results)) & !grepl("_ADJ", colnames(results))]

  # remove all existing column adjusted all pvalues
  results <- results[,!grepl("_ADJ_ALL_", colnames(results))]

  if (exists("results") & length(pvalue_columns)>0)
  {
    for (p in 1:length(pvalue_columns))
    {
      methods <- c("BY", "fdr","BH")
      for(method in 1:3)
      {
        col_p <- toupper(paste0(pvalue_columns[p], "_ADJ_ALL_", methods[method]))
        results[,col_p] <- stats::p.adjust(results[,pvalue_columns[p]],method  =  methods[method])
        colnames(results) <- toupper(colnames(results))
      }
    }

    if (nrow(results)>0)
      results <- results[order(results$PVALUE_ADJ_ALL_FDR),]

    if(filter_p_value)
      results <- subset(results, results$PVALUE < as.numeric(ssEnv$alpha) | results$PVALUE_ADJ_ALL_FDR < as.numeric(ssEnv$alpha))
  }


  results$DEPTH <- 3
  results[results$SUBAREA=="SAMPLE","DEPTH"] <- 1
  results[results$AREA_OF_TEST=="TOTAL","DEPTH"] <- 2
  # replace empty with NA
  results[results == ""] <- NA
  results[results == " "] <- NA
  # remove columns where all rows are NA
  results <- results[, colSums(is.na(results)) < nrow(results)]
  utils::write.csv2(results,fileNameResults , row.names  =  FALSE)


}

association_analysis_log <- function(inference_detail, start_time, end_time, processed_items)
{
  ssEnv <- get_session_info()
  log_folder <- ssEnv$session_folder
  association_file <- paste0(log_folder, "/association_analysis.csv")
  inference_detail$node_name <- as.character(Sys.info()["nodename"])
  inference_detail$session_id <- ssEnv$session_id
  # inference_detail$start_time <- as.character(format(start_time, "%a %b %d %X %Y"))
  # inference_detail$end_time <- as.character(format(end_time, "%a %b %d %X %Y"))
  inference_detail$processed_time <- as.numeric(difftime(end_time, start_time, units  =  "mins"))
  inference_detail$processed_items <- processed_items
  if(!file.exists(association_file))
  {
    write.table(inference_detail, file  =  association_file, sep  =  ",", row.names  =  FALSE, col.names  =  TRUE)
  } else
  {
    # convert all columns of inference detail as character
    tryCatch({
      association_file_data <- read.table(association_file, header  =  TRUE, sep  =  ",", stringsAsFactors  =  FALSE)
      association_file_data <- plyr::rbind.fill(inference_detail, association_file_data)
      write.table(association_file_data, file  =  association_file, sep  =  ",", row.names  =  FALSE, col.names  =  TRUE)
    }, error = function(e) {
      # print("ERROR: I'm stopping here, data to associate are not correct, file a bug!")
    })
  }

}
