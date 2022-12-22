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
  envir <- init_env( result_folder =  result_folder, maxResources =  maxResources, parallel_strategy  =  parallel_strategy, ...)


  figures <- envir$keys_figures[,1]
  anomalies <- envir$keys_anomalies[,1]
  populations <- c("Reference","Control","Case")

  if(sum(envir$keys_metaareas[,1]  ==  "PROBE")  ==  1)
  {
    subGroups <- c("")
    probes_prefix <- "PROBES"
    mainGroupLabel <-  "PROBE"
    subGroupLabel <- "GROUP"
    create_excel_pivot (envir = envir, populations  =   populations, figures  =   figures,anomalies  =   anomalies, subGroups  =   subGroups, probes_prefix  =    probes_prefix, mainGroupLabel  =   mainGroupLabel, subGroupLabel  =   subGroupLabel)
  }

  if(sum(envir$keys_metaareas[,1]  ==  "CHR")  ==  1)
  {
    subGroups <- c("CHR")
    probes_prefix <- "PROBES_CHR_"
    mainGroupLabel <-  "CHR"
    subGroupLabel <- "GROUP"
    create_excel_pivot (envir = envir, populations  =   populations, figures  =   figures,anomalies  =   anomalies, subGroups  =   subGroups, probes_prefix  =    probes_prefix, mainGroupLabel  =   mainGroupLabel, subGroupLabel  =   subGroupLabel)
  }

  if(sum(envir$keys_metaareas[,1]  ==  "GENE")  ==  1)
  {
    subGroups <- envir$gene_subareas[,1]
    probes_prefix  =  "PROBES_Gene_"
    mainGroupLabel  =   "GENE"
    subGroupLabel = "GROUP"

    create_excel_pivot (envir = envir, populations  =   populations, figures  =   figures,anomalies  =   anomalies, subGroups  =   subGroups, probes_prefix  =    probes_prefix, mainGroupLabel  =   mainGroupLabel, subGroupLabel  =   subGroupLabel)
  }


  if(sum(envir$keys_metaareas[,1]  ==  "ISLAND")  ==  1)
  {
    probes_prefix <- "PROBES_Island_"
    subGroups <- envir$island_subareas[,1]
    mainGroupLabel <- "ISLAND"
    subGroupLabel <- "RELATION_TO_CPGISLAND"
    create_excel_pivot (envir = envir, populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
  }

  if(sum(envir$keys_metaareas[,1]  ==  "DMR")  ==  1)
  {
    subGroups <- c("DMR")
    probes_prefix  =  "PROBES_DMR_"
    mainGroupLabel  =   "DMR"
    subGroupLabel = "GROUP"
    create_excel_pivot (envir = envir,populations, figures, anomalies, subGroups, probes_prefix, mainGroupLabel, subGroupLabel)
  }

  # if(exists("figures")) figures,  if(exists("anomalies"))  anomalies,if(exists("metaareas"))  metaareas

  inference_details <- unique(inference_details)
  # variables_to_export <- c("n", "working_data", "sig.formula", "tau", "lqm_control", "estimate", "independent_variable", "inference_details", "envir", "%dorng%", "k", "iter", "RNGseed", "checkRNGversion",
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
      message("Warning: one test family_test is missed! Skipped.")
    }
    else
    {
      transformation <- inference_detail$transformation
      if(is.null(transformation) || length(transformation)  ==  0)
      {
        transformation <- NULL
      }

      independent_variable <- gsub(" ","", inference_detail$independent_variable)
      if( is.null(independent_variable) || length(independent_variable)  ==  0)
      {
        message("Warning: one indipendent variable is missed! Skipped.")
      }
      else
      {
        depth_analysis <- inference_detail$depth_analysis
        if( is.null(depth_analysis) || length(depth_analysis)  ==  0)
        {
          depth_analysis <- 1
          message("Warning: missed depth analysis inference forced to 1.")
        }

        study_summary <-   utils::read.csv2(file_path_build( envir$result_folderData, "sample_sheet_result","csv"))

        # transform independent variable as factor
        if(family_test  ==  "binomial" || family_test  ==  "wilcoxon" || family_test  ==  "t.test")
          study_summary[,independent_variable] <- as.factor(study_summary[,independent_variable])

        result_folderPivot <- dir_check_and_create(envir$result_folderData,"Pivots")

        file_result_prefix <- paste0( depth_analysis,"_", as.character(independent_variable),"_",sep = "")

        study_summary <- subset(study_summary, study_summary$Sample_Group !=  "Reference")
        #########################################################################################################
        #########################################################################################################

        if (!(independent_variable %in% colnames(study_summary)))
        {
          message(" This indipendent variabile:", independent_variable, " is missed! Skipping")
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
            "ANOMALY"  =  "",
            "FIGURE"  =  "",
            "GROUP"  =  "",
            "SUBGROUP"  =  "",
            "AREA_OF_TEST"  =  "",
            "PVALUE"  =  "",
            "PVALUEADJ"  =  "",
            "TEST"  =  "",
            "BETA"  =  "",
            "AIC"  =  "",
            "RESIDUALS.SUM"  =  "",
            "FAMILY"  =  "",
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


          keys <- expand.grid("ANOMALY" =  envir$keys_anomalies[,1], "FIGURE" =  envir$keys_figures[,1])
          cols <- paste0(keys$ANOMALY,"_",keys$FIGURE, sep = "")
          keys$GROUP  =  "POPULATION"
          keys$SUBGROUP  =  "SAMPLE"
          iters <- length(cols)

          if(!is.null(covariates) && !length(covariates)  ==  0)
            study_summary <- study_summary[, c(independent_variable, covariates, cols) ]


          # define file results filename
          if(is.null(covariates) || length(covariates)  ==  0)
          {
            file_suffix <- "_test_result"
          }
          else
          {
            file_suffix <- "_test_corrected_result"
          }
          fileNameResults <- file_path_build(envir$result_folderInference,c(file_result_prefix , as.character(transformation), as.character(family_test), file_suffix),"csv")

          # clean keys from already done association
          if(file.exists(fileNameResults))
          {
            old_results <- utils::read.csv2(fileNameResults, header  =  T)
            keys_anomalies_figures_areas_done <- unlist(apply(unique(old_results [, c("ANOMALY","FIGURE","GROUP")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
            keys_to_be_done <- unlist(apply(keys[, c("ANOMALY","FIGURE","GROUP")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
            keys <- keys[!( keys_to_be_done %in% keys_anomalies_figures_areas_done), ]
          }

          if(nrow(keys)>0)
            for(j in seq_along(nrow(keys)))
            {
              g_start <- 2 + length(covariates)
              result_temp <- apply_stat_model(tempDataFrame  =  study_summary[, c(independent_variable, covariates, cols[j])], g_start  =  g_start , family_test  =  family_test, covariates  =  covariates,
                                              key  =  keys[j,], transformation  =  transformation, dototal  =  FALSE, logFolder =  envir$logFolder, independent_variable, depth_analysis, envir, ...)
              result <- rbind(result, result_temp)
            }

          result <- result[order(result$PVALUEADJ),]

          study_summaryToPlot <- study_summary

          if(independent_variable  ==  "Sample_Group")
            study_summaryToPlot$Sample_Group  <- stats::relevel(as.factor(study_summaryToPlot$Sample_Group), "Control")

          if(family_test  ==  "binomial" || family_test  ==  "wilcoxon")
          {

            chartFolder <- dir_check_and_create(envir$result_folderChart,"POPULATION_COMPARISON")

            if(sum(grepl(pattern = "MUTATIONS",x  =  envir$keys_anomalies))>0)
            {
              filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "MUTATIONS"),"png")
              grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
              graphics::par(mfrow = c(1,3))
              graphics::boxplot(MUTATIONS_HYPO~ study_summaryToPlot[,independent_variable],main = "Hypo Mutations", data  =  study_summaryToPlot, cex = 2)
              graphics::boxplot(MUTATIONS_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of Mutations", data  =  study_summaryToPlot, cex = 2)
              graphics::boxplot(MUTATIONS_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper Mutations", data  =  study_summaryToPlot, cex = 2)
              grDevices::dev.off()
            }

            if(sum(grepl(pattern = "LESIONS",x  =  envir$keys_anomalies))>0)
            {
              filename  =  file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "LESIONS"),"png")
              grDevices::png(file =  filename, width = 2480,height = 2480, pointsize  =  15, res  =  144)
              graphics::par(mfrow = c(1,3))
              graphics::boxplot(LESIONS_HYPO~study_summaryToPlot[,independent_variable],main = "Hypo Lesions", data  =  study_summaryToPlot, cex = 2)
              graphics::boxplot(LESIONS_BOTH~study_summaryToPlot[,independent_variable], main = "Both Type of Lesions", data  =  study_summaryToPlot, cex = 2)
              graphics::boxplot(LESIONS_HYPER~study_summaryToPlot[,independent_variable],main = "Hyper Lesions", data  =  study_summaryToPlot, cex = 2)
              grDevices::dev.off()
            }

            if(sum(grepl(pattern = "DELTAS",x  =  envir$keys_anomalies))>0)
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

            if(sum(grepl(pattern = "DELTAQ",x  =  envir$keys_anomalies))>0)
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

            keys <- envir$keys_anomalies_figures_areas
            # clean keys from already done association
            if(file.exists(fileNameResults))
            {
              old_results <- utils::read.csv2(fileNameResults, header  =  T)
              keys_anomalies_figures_areas_done <- unlist(apply(unique(old_results [, c("ANOMALY","FIGURE","GROUP")]), 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
              keys_to_be_done <- unlist(apply(keys[, c("ANOMALY","FIGURE","GROUP")], 1, function(x) paste(x, collapse  =  "_", sep  =  "")))
              keys <- keys[!( keys_to_be_done %in% keys_anomalies_figures_areas_done), ]
            }

            nkeys <- nrow(keys)
            variables_to_export_nested <- c("variables_to_export", "keys", "result_folderPivot", "sample_names", "independent_variable", "covariates", "family_test", "transformation", "envir", "depth_analysis","file_path_build", "apply_stat_model")
            if(nrow(keys)>0)
              # result_temp_foreach <- foreach::foreach(k  =  1:nkeys, .combine  =  rbind, .export  =  variables_to_export_nested) %dorng%
              for (k in 1:nkeys)
              {
                if(exists("tempDataFrame"))
                  rm(list  =  c("tempDataFrame"))
                key <- keys [k,]
                pivot_subfolder <- dir_check_and_create(result_folderPivot, key$ANOMALY)
                fname <-file_path_build( pivot_subfolder ,c(key$ANOMALY, key$FIGURE, key$GROUP,key$SUBGROUP),"csv")
                if (file.exists(fname))
                {
                  message(Sys.time()," Starting to read pivot:", fname,".")
                  tempDataFrame <- utils::read.csv(fname, sep  =  ";")
                  #assign the area name (eg gene...) to the rows
                  row.names(tempDataFrame) <- tempDataFrame$SAMPLEID
                  #removes area name (eg. gene...)
                  tempDataFrame <- tempDataFrame[,-1]
                  # max_row_count <- ceiling(10^6/ncol(tempDataFrame))
                  # batch_count <- ceiling(nrow(tempDataFrame)/max_row_count)
                  # for(h in 0:batch_count)
                  # {
                  # tt <- tempDataFrame[ (1+h*max_row_count):min((h+1)*max_row_count,nrow(tempDataFrame)), ]
                  tt <- tempDataFrame
                  tt <- t(tt)
                  tt <- as.data.frame(tt)
                  tt$Sample_ID <- rownames(tt)
                  message(Sys.time()," Read pivot:", fname, " with ", ncol(tt), " rows.")
                  if(nrow(tt)>1)
                  {
                    tt <- subset(tt, "POPULATION"  !=   "Reference")
                    tt <- subset(tt, "POPULATION"  !=   0)
                    tt <-  merge( x  =   sample_names, y  =   tt,  by.x  =  "Sample_ID",  by.y  =  "Sample_ID" , all.x  =  TRUE)
                    tt <- as.data.frame(tt)
                    tt$POPULATION <- sample_names[, independent_variable]
                    tt[is.na(tt)] <- 0
                    tt[, independent_variable] <- sample_names[, independent_variable]
                    tt <- tt[, !(names(tt) %in% c("POPULATION","Sample_ID"))]
                    cols <- (gsub(" ", "_", colnames(tt)))
                    cols <- (gsub("-", "_", cols))
                    cols <- (gsub(":", "_", cols))
                    cols <- (gsub("/", "_", cols))
                    cols <- (gsub("'", "_", cols))
                    tt <- as.data.frame(tt)
                    if(length(colnames(tt)) !=  length(cols))
                      browser()
                    colnames(tt) <- cols
                    g_start <- 2 + length(covariates)
                    result_temp_local_batch <- apply_stat_model(tempDataFrame  =  tt, g_start  =  g_start, family_test  =  family_test, covariates  =  covariates,
                                                                key  =  key, transformation =  transformation, dototal  =  TRUE,
                                                                logFolder =  envir$logFolder, independent_variable, depth_analysis, envir, ...)

                    if(!exists("result_temp_foreach"))
                      result_temp_foreach <- result_temp_local_batch
                    else
                      result_temp_foreach <- rbind(result_temp_local_batch, result_temp_foreach)
                    # result_temp_local_batch
                  }
                  # }
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
    result <- result[order(result$PVALUEADJ),]

    if(filter_p_value)
      result <- subset(result, result$PVALUE < 0.05 | result$PVALUEADJ < 0.05)

    if(nrow(result)>0)
    {
      if(file.exists(fileNameResults))
      {
        old_results <- utils::read.csv2(fileNameResults, header  =  T)
        result <- rbind(result, old_results)
      }
      result[,"PVALUEADJ_ALL_BH"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BH")
      result[,"PVALUEADJ_ALL_BY"] <- stats::p.adjust(result[,"PVALUE"],method  =  "BY")
      result <- result[order(result$PVALUEADJ),]
      utils::write.csv2(result,fileNameResults , row.names  =  FALSE)
    }

  }

  future::plan( future::sequential)
}
