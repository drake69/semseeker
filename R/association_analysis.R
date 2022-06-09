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
#' (mutations and lesions): scale, log, log2, log10, exp, johnson, none, quantile_quantiles(as number) eg quantile_3
#' depth analysis:
#' 1: sample level
#' 2: type level (gene, DMR, cpgisland) (includes 1)
#' 3: genomic area: gene, body, gene tss1550, gene whole, gene tss200,  (includes 1 and 2)
#' filter_p_value report after adjusting saves only significant nominal p-value
#' @param result_folder where semseeker's results are stored, the root folder
#' @param maxResources percentage of max system's resource to use
#'
#' @return
#' @export
#'
#' @examples
association_analysis <- function(inference_details,result_folder, maxResources)
{

  envir <- init_env(result_folder, maxResources)

  # covariates, family_test, transformation = NULL, independent_variable, depth_analysis=3
  # covariates= covariates, transformation = "log",independent_variable = "Sample_Group", depth_analysis=2
  inference_details <- unique(inference_details)
  for(i in 1:nrow(inference_details))
  {
    #browser()
    inference_detail <- inference_details[i,]

    filter_p_value <- if(!is.null(inference_detail$filter_p_value)) inference_detail$filter_p_value else TRUE

    covariates <- inference_detail$covariates
    covariates <- if(length(covariates)!=0 & !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split = "+", fixed = T)))
    family_test <- inference_detail$family_test
    if( is.null(family_test) || length(family_test)==0)
    {
      message("Warning: one test family_test is missed! Skipped.")
      #browser()
      next
    }
    # browser()
    transformation <- inference_detail$transformation
    if(is.null(transformation) | length(transformation)==0)
    {
      transformation <- NULL
    }

    independent_variable <- gsub(" ","", inference_detail$independent_variable)
    if( is.null(independent_variable) || length(independent_variable)==0)
    {
      message("Warning: one indipendent variable is missed! Skipped.")
      #browser()
      next
    }

    depth_analysis <- inference_detail$depth_analysis
    if( is.null(depth_analysis) || length(depth_analysis)==0)
    {
      depth_analysis <- 1
      message("Warning: missed depth analysis inference forced to 1.")
    }

    study_summary <-   utils::read.csv2(file_path_build( envir$result_folderData, "sample_sheet_result","csv"))

    # transform independent variable as factor
    if(family_test=="binomial" | family_test=="wilcoxon" | family_test=="t.test")
      study_summary[,independent_variable] <- as.factor(study_summary[,independent_variable])

    result_folderPivot <- dir_check_and_create(envir$result_folderData,"Pivots")

    file_result_prefix <- paste0( depth_analysis,"_", as.character(independent_variable),"_",sep="")

    study_summary <- subset(study_summary, study_summary$Sample_Group!="Reference")
    #########################################################################################################
    #########################################################################################################

    if (!(independent_variable %in% colnames(study_summary)))
    {
      message(" This indipendent variabile:", independent_variable, " is missed! Skipping")
      #browser()
      next
    }

    if(is.null(covariates) || length(covariates)==0)
    {
      sample_names <- data.frame(study_summary[, c("Sample_ID", independent_variable)])
    } else
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
    # sample_names <- data.frame(study_summary)
    result = data.frame (
      "INDIPENDENT.VARIABLE"="",
      "ANOMALY" = "",
      "FIGURE" = "",
      "GROUP" = "",
      "SUBGROUP" = "",
      "AREA_OF_TEST" = "",
      "PVALUE" = "",
      "PVALUEADJ" = "",
      "TEST" = "",
      "BETA" = "",
      "AIC" = "",
      "RESIDUALS.SUM" = "",
      "FAMILY" = "",
      "TRANSFORMATION" = "",
      "COVARIATES" = "",
      "SHAPIRO.PVALUE" = "",
      "BARTLETT.PVALUE" = "",
      "CASE.LABEL"="",
      "COUNT.CASE","",
      "MEAN.CASE" ="",
      "SD.CASE"="",
      "CONTROL.LABEL"="",
      "COUNT.CONTROL"="",
      "MEAN.CONTROL"="",
      "SD.CONTROL"="",
      "RHO"="",
      "CI.LOWER"="",
      "CI.UPPER"="",
      "N.PERMUTATIONS" =""
    )
    result <- result[-1,]

    # resultRun <- result

    if(!is.null(covariates) && !length(covariates)==0)
    {
      sample_names <-   sample_names[, c("Sample_ID", independent_variable, covariates)]
      colnames(sample_names) <- c("Sample_ID", independent_variable, covariates)
    } else
    {
      sample_names <-   sample_names[, c("Sample_ID", independent_variable)]
      colnames(sample_names) <- c("Sample_ID", independent_variable)
    }
    ######################################################################################################
    # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
    # binomiale che si vuole usare per la regressione logistica


    keys <- expand.grid("ANOMALY"=c("LESIONS","MUTATIONS"), "FIGURE"=c("HYPER","HYPO","BOTH"))
    cols <- paste0(keys$ANOMALY,"_",keys$FIGURE, sep="")
    keys$GROUP = "POPULATION"
    keys$SUBGROUP = "SAMPLE"
    iters <- length(cols)

    if(!is.null(covariates) && !length(covariates)==0)
      study_summary <- study_summary[, c(independent_variable, covariates, cols) ]

    # #browser()

    for(i in 1:nrow(keys))
    {
      # i <- 1
      # family_test <- "poisson"
      # transformation <- NULL
      g_start <- 2 + length(covariates)
      result_temp <- apply_stat_model(tempDataFrame = study_summary[, c(independent_variable, covariates, cols[i])], g_start = g_start , family_test = family_test, covariates = covariates,
                                 key = keys[i,], transformation = transformation, dototal = FALSE, logFolder= envir$logFolder, independent_variable, depth_analysis)
      # browser()
      result <- rbind(result, result_temp)
    }

    result <- result[order(result$PVALUEADJ),]

    # #browser()
    study_summaryToPlot <- study_summary

    if(independent_variable=="Sample_Group")
       study_summaryToPlot$Sample_Group  <- stats::relevel(as.factor(study_summaryToPlot$Sample_Group), "Control")

    if(family_test=="binomial")
    {

      chartFolder <- dir_check_and_create(envir$result_folderChart,"POPULATION_COMPARISON")

      filename = file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "MUTATIONS"),"png")
      grDevices::png(file= filename, width=2480,height=2480, pointsize = 15, res = 144)
      graphics::par(mfrow=c(1,3))
      graphics::boxplot(MUTATIONS_HYPO~ study_summaryToPlot[,independent_variable],main="Hypo Mutations", data = study_summaryToPlot, cex=2)
      graphics::boxplot(MUTATIONS_BOTH~study_summaryToPlot[,independent_variable], main="Both Type of Mutations", data = study_summaryToPlot, cex=2)
      graphics::boxplot(MUTATIONS_HYPER~study_summaryToPlot[,independent_variable],main="Hyper Mutations", data = study_summaryToPlot, cex=2)
      grDevices::dev.off()

      filename = file_path_build(chartFolder,c(file_result_prefix,as.character(transformation), "LESIONS"),"png")
      grDevices::png(file= filename, width=2480,height=2480, pointsize = 15, res = 144)
      graphics::par(mfrow=c(1,3))
      graphics::boxplot(LESIONS_HYPO~study_summaryToPlot[,independent_variable],main="Hypo Lesions", data = study_summaryToPlot, cex=2)
      graphics::boxplot(LESIONS_BOTH~study_summaryToPlot[,independent_variable], main="Both Type of Lesions", data = study_summaryToPlot, cex=2)
      graphics::boxplot(LESIONS_HYPER~study_summaryToPlot[,independent_variable],main="Hyper Lesions", data = study_summaryToPlot, cex=2)
      grDevices::dev.off()
    }

    # browser()
    if(depth_analysis >1)
    {

      keys <- envir$keys_anomalies_figures_areas
      # #browser()
      # areas <- rbind(genes, islands, dmrs)
      nkeys <- nrow(keys)
      # #browser()

      # parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c("apply_stat_model","file_path_build"))
      to_export <- c("keys", "file_path_build", "result_folderPivot", "sample_names", "independent_variable", "covariates", "apply_stat_model", "family_test",
                     "transformation", "envir", "depth_analysis","BCApval")

      result_temp_foreach <- foreach::foreach(i = 1:nkeys, .combine = rbind, .export = to_export) %dorng%
      # for (i in 1:nkeys)
      {
        # i <- 25
        if(exists("tempDataFrame"))
          rm(list = c("tempDataFrame"))
        key <- keys [i,]
        # print(key)
        fname <-file_path_build( result_folderPivot ,c(key$ANOMALY, key$FIGURE, key$GROUP,key$SUBGROUP),"csv")
        # print(fname)
        if (file.exists(fname))
        {
          tempDataFrame <- utils::read.csv(fname, sep = ";")
          row.names(tempDataFrame) <- tempDataFrame$X
          tempDataFrame <- tempDataFrame[,-1]
          tempDataFrame <- t(tempDataFrame)
          if(nrow(tempDataFrame)>1)
          {
            message(nrow(tempDataFrame))
            # tempDataFrame <- tempDataFrame[,-2]
            tempDataFrame <- as.data.frame(tempDataFrame)
            tempDataFrame <- subset(tempDataFrame, "POPULATION" != "Reference")
            tempDataFrame <- subset(tempDataFrame, "POPULATION" != 0)

            # #browser()
            tempDataFrame <-  merge( x =  sample_names, y =  tempDataFrame,  by.x = "Sample_ID",  by.y = "SAMPLEID" , all.x = TRUE)
            tempDataFrame$POPULATION <- sample_names[, independent_variable]
            tempDataFrame[is.na(tempDataFrame)] <- 0
            #  we want to preserve the NA in the indipendent variables to be removed by the models
            tempDataFrame[, independent_variable] <- sample_names[, independent_variable]
            tempDataFrame <- tempDataFrame[, !(names(tempDataFrame) %in% c("POPULATION","Sample_ID"))]

            cols <- (gsub(" ", "_", colnames(tempDataFrame[])))
            cols <- (gsub("-", "_", cols))
            cols <- (gsub(":", "_", cols))
            cols <- (gsub("/", "_", cols))
            cols <- (gsub("'", "_", cols))

            colnames(tempDataFrame) <- cols

            tempDataFrame <- as.data.frame(tempDataFrame)
            # tempDataFrame$Sample_Group <- as.factor(tempDataFrame$Sample_Group)
            # tempDataFrame$Sample_Group  <- stats::relevel(tempDataFrame$Sample_Group, "Control")

            g_start <- 2 + length(covariates)

            result_temp_local <- apply_stat_model(tempDataFrame = tempDataFrame, g_start = g_start, family_test = family_test, covariates = covariates, key = key, transformation= transformation, dototal = TRUE,
                                            logFolder= envir$logFolder, independent_variable =  independent_variable, depth_analysis = depth_analysis)

            # #browser()
            # n_adj <- iters - g_start
            # result <- rbind(result, result_temp)
            result_temp_local
          }
        }
      }
    }

    if(exists("result_temp_foreach"))
      result <- rbind(result, result_temp_foreach)

    # browser()
    result <- unique(result)
    result[,"PVALUEADJ_ALL"] <- round(stats::p.adjust(result[,"PVALUE"],method = "BH"),5)

    result <- result[order(result$PVALUEADJ),]

    if(filter_p_value)
        result <- subset(result, result$PVALUE < 0.05 | result$PVALUEADJ < 0.05)

    if(nrow(result)>0)
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        file_suffix <- "_test_result"
      } else
      {
        file_suffix <- "_test_corrected_result"
      }

      # resultRun <- rbind(result, resultRun)
      fileName <- file_path_build(envir$result_folderInference,c(file_result_prefix , as.character(transformation), as.character(family_test), file_suffix),"csv")
      utils::write.csv2(result,fileName , row.names = FALSE)
    }

    # fileName <- file_path_build(result_folder_inference,"20220326","csv")
    # write.csv2(resultRun,fileName , row.names = FALSE)

    # res.pvalue <- subset(result, PVALUE < 0.05)
    # res.pvalue$beta_gt1 <- res.pvalue$BETA>1
    # res.pvalue$beta_gt1 <- as.numeric(res.pvalue$beta_gt1)

    # source("/home/lcorsaro/Desktop/Progetti/r-studio/smarties/R/microarray/epigenetics/epimutation_analysis/qqplot_inferential.R")
    # result <- utils::read_csv(file.path(result_folder_inference,paste0(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")))
    # qqunif.plot(diffMethTable_site_cmp1$diffmeth.p.val, result_folder_inference =  report.dir, filePrefix ="diff_meth_sites")

    # case_vs_control_binomial_regression_corrected_result <-
    #   read.csv2(
    #     "/home/lcorsaro/Desktop/experiments_data/DIOSSINA_DESIO/3_semseeker_result/Pivots/case_vs_control_binomial_regression_corrected_result_1.csv"
    #   )
    # keys <- unique(result[, c("ANOMALY", "FIGURE", "GROUP", "SUBGROUP")])
    #
    # for (i in 1:dim(keys)[1])
    # {
    #   # i <- 2
    #   key <- keys[i, ]
    #   diffmeth.p.val <-
    #     subset(result,
    #            ANOMALY == key$ANOMALY)
    #   diffmeth.p.val <-
    #     subset(diffmeth.p.val, FIGURE == key$FIGURE)
    #   diffmeth.p.val <- subset(diffmeth.p.val, GROUP == key$GROUP)
    #
    #   diffmeth.p.val <-
    #     subset(diffmeth.p.val, SUBGROUP == key$SUBGROUP)
    #
    #   diffmeth.p.val <- subset(diffmeth.p.val,PVALUE !=0 )
    #   if (dim(diffmeth.p.val)[1] <= 1)
    #     return
    #   #######inserisco nella funzione i pvalues ottenuti dalla differential (non aggiustati)
    #
    #   file_prefix <- paste0("case_vs_control_binomial_regression_corrected_result","_", key$ANOMALY,"_", key$FIGURE,"_", key$GROUP,"_", key$SUBGROUP,"_", sep="")
    #   qqunifPlot(diffmeth.p.val$PVALUE,
    #               result_folder_inference = result_folder_inference,
    #               filePrefix = file_prefix)
    # }
    #
    #
    #
    # qqunif.plot(pvalues)
  }
}
