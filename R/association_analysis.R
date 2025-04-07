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
association_analysis <- function(inference_details,result_folder, maxResources = 90,
  parallel_strategy  = "multicore",start_fresh = FALSE, ...)
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

  log_event("BANNER: ", format(Sys.time(), "%a %b %d %X %Y"), " SemSeeker will perform the association analysys for project \n in ", ssEnv$result_folderData)

  if(start_fresh)
    unlink(ssEnv$result_folderInference, recursive = TRUE)
  dir_check_and_create(ssEnv$result_folderInference,c())

  localKeys <- ssEnv$keys_markers_figures
  sample_groups <- c("Reference","Control","Case")

  deltaX_get()
  annotate_position_pivots()

  inference_details <- unique(inference_details)
  # variables_to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable", "inference_details", "ssEnv", "%dorng%", "k", "iter", "RNGseed", "checkRNGversion",
  #                          "getRNG", "%||%", ".getDoParName", "getDoParName", "getDoBackend", "setDoBackend", "RNGtype", "showRNG", ".getRNGattribute", "isNumber", "isReal", "isInteger", ".foreachGlobals", "RNGprovider", ".RNGkind_length", "tail",
  #                          "file_path_build", "dir_check_and_create", "apply_stat_model", "doRNGversion", ".getRNG", "hasRNG", "nextRNG", "RNGkind", "setRNG", "RNGstr")

  # foreach::foreach(z  =  1:nrow(inference_details), .export  =  variables_to_export) %dorng%
  for(z in 1:nrow(inference_details))
  {
    results <- data.frame()
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
    family_test <- as.character(inference_detail$family_test)
    study_summary <-   study_summary_get(inference_detail$samples_sql_condition)

    # dummify covariates dummy
    # covariates <- dummify_covariates(covariates, study_summary)

    # browser()
    results <- data.frame()
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


        # transform independent variable as factor
        if(family_test  ==  "binomial" || family_test  ==  "wilcoxon" || family_test  ==  "t.test")
          study_summary[,independent_variable] <- as.factor(study_summary[,independent_variable])

        file_result_prefix <- paste(depth_analysis, as.character(independent_variable),sep = "_")

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

          # remove samples with missed values
          sample_names <- sample_names[complete.cases(sample_names),]
          if(nrow(sample_names)  ==  0)
          {
            log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " No samples with complete data for the analysis! Skipped.")
            next
          }

          ######################################################################################################
          # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
          # binomiale che si vuole usare per la regressione logistica
          #
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
          for (a in 1:length(markers))
          {
            results <- data.frame()
            keys <- localKeys[localKeys$MARKER==markers[a],]
            keys <- unique(keys)
            cols <- keys$COMBINED
            fileNameResults <- inference_file_name(inference_detail, markers[a], ssEnv$result_folderInference,prefix= ifelse(areas_selection==c(),"",paste(areas_selection, "_", sep = "")))
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

              #
              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  > 3
              old_results <- data.frame()
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
                  if(any(is.na(study_summary_local[, column_selectors])))
                  {
                    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missing values in the data frame!")
                    # remove rows with missing values
                    study_summary_local <- study_summary_local[complete.cases(study_summary_local[, column_selectors]),]
                    # browser()
                  }
                  result_temp <- apply_stat_model(tempDataFrame  =  study_summary_local[, column_selectors], g_start  =  g_start , family_test  =  family_test, covariates  =  covariates,
                    key  =  key, transformation  =  transformation, dototal  =  FALSE, session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,inference_detail$samples_sql_condition,  ...)
                  results <- plyr::rbind.fill(results, result_temp)
                }

              if (!is.null(dim(results)))
                if (nrow(results)>0)
                  if("PVALUE_ADJ" %in% colnames(results))
                    results <- results[order(results$PVALUE_ADJ),]

              if (exists("old_results"))
              {
                results <- plyr::rbind.fill(results, old_results)
                rm(old_results)
              }

              results[results == ""] <- NA
              # remove columns where all rows are NA
              results <- results[, colSums(is.na(results)) < nrow(results)]
              # utils::write.csv2(results, fileNameResults , row.names  =  FALSE)
            }


            # execute for all the areas
            if(depth_analysis>1)
            {
              localKeys_1 <- ssEnv$keys_areas_subareas_markers_figures
              keys <- localKeys_1[localKeys_1$MARKER==markers[a],]

              file_good <- file.exists(fileNameResults) && file.info(fileNameResults)$size  >10
              dototal <- TRUE
              old_results <- data.frame()
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
              variables_to_export_nested <- c("variables_to_export", "keys", "sample_names", "independent_variable", "covariates",
                "family_test", "transformation", "ssEnv", "depth_analysis","file_path_build", "apply_stat_model")
              if(nrow(keys)>0)
                # result_temp_foreach <- foreach::foreach(k  =  1:nkeys, .combine  =  rbind, .export  =  variables_to_export_nested) %dorng%
                for (k in 1:nkeys)
                {
                  if(exists("tempDataFrame"))
                    rm(list  =  c("tempDataFrame"))
                  key <- keys [k,]
                  if (key$AREA=="POSITION")
                    next
                  pivot_filename <- pivot_file_name_parquet(key$MARKER, key$FIGURE, key$AREA, key$SUBAREA)
                  if (file.exists(pivot_filename))
                  {
                    areas_selection_temp <- areas_selection
                    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting to read pivot:", pivot_filename,".")
                    tempDataFrame <- arrow::read_parquet(pivot_filename)
                    log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Read pivot:", pivot_filename, " with ", nrow(tempDataFrame), " rows.")
                    tempDataFrame[is.na(tempDataFrame)] <- 0
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
                        tempDataFrame <- tempDataFrame[ tempDataFrame[,1] %in% areas_selection_temp, ]
                      }
                    }
                    if(is.null(dim(tempDataFrame)))
                      next
                    if(plyr::empty(tempDataFrame) | nrow(tempDataFrame)==0)
                      next

                    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"), " Starting to execute required test for:", key$MARKER, key$FIGURE, key$AREA, key$SUBAREA,".")

                    chunk_size <- 50000  # Define a chunk size
                    for (i in seq(1, nrow(tempDataFrame), by = chunk_size)) {
                      {
                        chunk_indices <- i:min(i + chunk_size - 1, nrow(tempDataFrame))
                        tempDataFrameBatch <- as.data.frame(tempDataFrame)[chunk_indices,]
                        rownames(tempDataFrameBatch) <- tempDataFrameBatch[,1]
                        tempDataFrameBatch <- tempDataFrameBatch[,-1]
                        tempDataFrameBatch <- t(tempDataFrameBatch)
                        tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                        tempDataFrameBatch$Sample_ID <- rownames(tempDataFrameBatch)
                        log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Transposed pivot:", pivot_filename, " with ", ncol(tempDataFrameBatch) -1 , " columns.")
                        if(nrow(tempDataFrameBatch)>1)
                        {
                          tempDataFrameBatch <-  merge( x  =   sample_names, y  =   tempDataFrameBatch,  by.x  =  "Sample_ID",  by.y  =  "Sample_ID" , all.x  =  TRUE)
                          log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Merged pivot:", pivot_filename, " with ", ncol(tempDataFrameBatch), " columns.")
                          tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                          tempDataFrameBatch[is.na(tempDataFrameBatch)] <- 0
                          tempDataFrameBatch <- tempDataFrameBatch[,-1]
                          cols <- (colnames(tempDataFrameBatch))
                          tempDataFrameBatch <- as.data.frame(tempDataFrameBatch)
                          if(length(colnames(tempDataFrameBatch)) !=  length(cols))
                            stop("ERROR: I'm stopping here data to associate are not correct, file a bug!")
                          colnames(tempDataFrameBatch) <- cols
                          # two are the sample_id and the sample_group
                          g_start <- 2 + length(covariates)
                          processed_items <- processed_items + ncol(tempDataFrameBatch) - g_start
                          if(any(is.na(tempDataFrameBatch)))
                          {
                            log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " Missing values in the data frame!")
                            # next
                          }
                          result_temp_local_batch <- apply_stat_model(tempDataFrame  =  tempDataFrameBatch, g_start  =  g_start, family_test  =  family_test, covariates  =  covariates,
                            key  =  key, transformation =  transformation, dototal  =  dototal,
                            session_folder =  ssEnv$session_folder, independent_variable, depth_analysis,inference_detail$samples_sql_condition,  ...)

                          results <- plyr::rbind.fill(results, result_temp_local_batch)
                        }
                      }
                    }
                  }
                  else
                    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"), " File not found:", pivot_filename,".")
                  association_analysis_log(cbind(inference_detail,keys[k,]), start_time, Sys.time(), processed_items)
                  if(nrow(results)!=0)
                    results <- subset(results, MARKER == key$MARKER)
                }

              if (exists("old_results"))
              {
                results <- plyr::rbind.fill(results, old_results)
                rm(old_results)
              }
            }

            if(exists("result_temp_foreach"))
              rm(result_temp_foreach)
            if(exists("result_temp_local_batch"))
              rm(result_temp_local_batch)
            if(exists("tempDataFrame"))
              rm(tempDataFrame)

            # remove any column containing in the name samples_sql_condition
            results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
            # results$samples_sql_condition <- inference_detail$samples_sql_condition
            association_analysis_save_results(results,fileNameResults, family_test, filter_p_value )

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

association_analysis_save_results <- function(results=NULL,fileNameResults, family_test, filter_p_value ){


  if(nrow(results)==0)
    return()

  ssEnv <- get_session_info()

  # there is a bug which mantain more family test in the same results file
  # so we need to filter the results
  # browser()
  colnames(results) <- toupper(colnames(results))
  results <- subset(results, FAMILY.TEST==as.character(family_test))

  # check if results is empty
  if(is.null(results))
    return()

  if(nrow(results)==0)
    return()

  results <- results[,!grepl("SAMPLES_SQL_CONDITION", colnames(results))]
  results <- unique(results)

  pvalue_columns <- colnames(results)[grepl("PVALUE", colnames(results)) & !grepl("_ADJ", colnames(results))]

  # remove all existing column adjusted all pvalues
  results <- results[,!grepl("_ADJ_ALL_", colnames(results))]

  if (exists("results") & length(pvalue_columns)>0)
  {
    for (p in 1:length(pvalue_columns))
    {
      col_p <- toupper(paste0(pvalue_columns[p], "_ADJ_ALL_", ssEnv$multiple_test_adj))
      results[,col_p] <- stats::p.adjust(results[,pvalue_columns[p]],method  =  ssEnv$multiple_test_adj)
      colnames(results) <- toupper(colnames(results))
    }

    pvalue_adj_colname <- name_cleaning(paste0("PVALUE_ADJ_ALL_", ssEnv$multiple_test_adj))

    if (nrow(results)>0)
      results <- results[order(results[,pvalue_adj_colname]),]

    if(filter_p_value)
      results <- subset(results, results$PVALUE < as.numeric(ssEnv$alpha) | results[,pvalue_adj_colname] < as.numeric(ssEnv$alpha))
  }

  if(nrow(results)==0)
    return()


  results$DEPTH <- 3
  # replace NA of SUBAREA with TOTAL
  results[is.na(results$SUBAREA),"SUBAREA"] <- "TOTAL"
  results[results$SUBAREA=="SAMPLE","DEPTH"] <- 1
  selector <- grepl("TOTAL",results$AREA_OF_TEST)
  results[selector,"DEPTH"] <- 2
  # replace empty with NA
  results[results == ""] <- NA
  results[results == " "] <- NA
  # remove columns where all rows are NA
  results <- results[, colSums(is.na(results)) < nrow(results)]

  # check if exists at least a column with PVALUE
  if(!any(grepl("PVALUE", colnames(results))))
    return()

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
    utils::write.csv2(inference_detail, file  =  association_file, row.names  =  FALSE, col.names  =  TRUE)
  } else
  {
    # convert all columns of inference detail as character
    tryCatch({
      association_file_data <- utils::read.csv2(association_file, header  =  TRUE, stringsAsFactors  =  FALSE)
      association_file_data <- plyr::rbind.fill(inference_detail, association_file_data)
      utils::write.csv2(association_file_data, file  =  association_file, row.names  =  FALSE, col.names  =  TRUE)
    }, error = function(e) {
      # print("ERROR: I'm stopping here, data to associate are not correct, file a bug!")
    })
  }

}
