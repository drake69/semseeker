inferenceAnalysis <- function(inferenceDetails)
{

  # covariates, family_test, transformation = NULL, independentVariable, depthAnalysis=3
  # covariates= covariates, transformation = "log",independentVariable = "Sample_Group", depthAnalysis=2
  inferenceDetails <- unique(inferenceDetails)
  for(i in 1:nrow(inferenceDetails))
  {
    #browser()
    inferenceDetail <- inferenceDetails[i,]

    filterpvalue <- if(!is.null(inferenceDetail$filterpvalue)) inferenceDetail$filterpvalue else TRUE

    covariates <- inferenceDetail$covariates
    covariates <- if(length(covariates)!=0 & !is.null(covariates)) unlist(t(strsplit( gsub(" ","",covariates),split = "+", fixed = T)))
    family_test <- inferenceDetail$family_test
    if( is.null(family_test) || length(family_test)==0)
    {
      message("Warning: one test family_test is missed! Skipped.")
      #browser()
      next
    }
    # browser()
    transformation <- inferenceDetail$transformation
    if(is.null(transformation) | length(transformation)==0)
    {
      transformation <- NULL
    }

    independentVariable <- gsub(" ","", inferenceDetail$independentVariable)
    if( is.null(independentVariable) || length(independentVariable)==0)
    {
      message("Warning: one indipendent variable is missed! Skipped.")
      #browser()
      next
    }

    depthAnalysis <- inferenceDetail$depthAnalysis
    if( is.null(depthAnalysis) || length(depthAnalysis)==0)
    {
      depthAnalysis <- 1
      message("Warning: missed depth analysis inference forced to 1.")
    }

    studySummary <-   utils::read.csv2(file_path_build(resultFolderData, "sample_sheet_result","csv"))
    resultFolderDataPivot <- dir_check_and_create(resultFolderData,"Pivots")

    file_result_prefix <- paste0( depthAnalysis,"_", independentVariable,"_",sep="")

    studySummary <- subset(studySummary, Sample_Group!="Reference")
    #########################################################################################################
    #########################################################################################################

    if(!dir.exists(resultFolderDataInference))
      dir.create(resultFolderDataInference)

    if (!(independentVariable %in% colnames(studySummary)))
    {
      message(" This indipendent variabile:", independentVariable, " is missed! Skipping")
      #browser()
      next
    }

    if(is.null(covariates) || length(covariates)==0)
    {
      sample_names <- data.frame(studySummary[, c("Sample_ID", independentVariable)])
    } else
    {
      sample_names <- data.frame(studySummary[, c("Sample_ID", independentVariable, covariates)])
    }
    # sample_names <- data.frame(studySummary)
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
      "COUNT.CASE","",
      "MEAN.CASE" ="",
      "SD.CASE"="",
      "COUNT.CONTROL"="",
      "MEAN.CONTROL"="",
      "SD.CONTROL"="",
      "RHO"=""
    )
    result <- result[-1,]



    if(!is.null(covariates) && !length(covariates)==0)
    {
      sample_names <-   sample_names[, c("Sample_ID", independentVariable, covariates)]
      colnames(sample_names) <- c("Sample_ID", independentVariable, covariates)
    } else
    {
      sample_names <-   sample_names[, c("Sample_ID", independentVariable)]
      colnames(sample_names) <- c("Sample_ID", independentVariable)
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
      studySummary <- studySummary[, c(independentVariable, covariates, cols) ]

    # #browser()

    for(i in 1:nrow(keys))
    {
      # i <- 1
      # family_test <- "poisson"
      # transformation <- NULL
      g_start <- 2 + length(covariates)
      result_temp <- apply_model(tempDataFrame = studySummary[, c(independentVariable, covariates, cols[i])], g_start = g_start , family = family_test, covariates = covariates,
                                 key = keys[i,], transformation = transformation, dototal = FALSE, logFolder= logFolder, independentVariable, depthAnalysis)
      result <- rbind(result, result_temp)
    }

    result <- result[order(result$PVALUEADJ),]

    # #browser()
    studySummaryToPlot <- studySummary
    studySummaryToPlot$Sample_Group  <- stats::relevel(as.factor(studySummaryToPlot$Sample_Group), "Control")

    if(family_test=="binomial")
    {

      chartFolder <- dir_check_and_create(resultFolderDataChart,"POPULATION_COMPARISON")

      filename <- file_path_build(chartFolder,file_result_prefix, "MUTATIONS.png")
      grDevices::png(file= filename, width=2480,height=2480)
      graphics::par(mfrow=c(1,3))
      graphics::boxplot(MUTATIONS_HYPO~ studySummaryToPlot[,independentVariable],main="Hypo Mutations", data = studySummaryToPlot, cex=2)
      graphics::boxplot(MUTATIONS_BOTH~studySummaryToPlot[,independentVariable], main="Both Type of Mutations", data = studySummaryToPlot, cex=2)
      graphics::boxplot(MUTATIONS_HYPER~studySummaryToPlot[,independentVariable],main="Hyper Mutations", data = studySummaryToPlot, cex=2)
      grDevices::dev.off()

      filename = paste0( chartFolder,"/", file_result_prefix, "LESIONS.png",sep="")
      grDevices::png(file= filename, width=2480,height=2480)
      graphics::par(mfrow=c(1,3))
      graphics::boxplot(LESIONS_HYPO~studySummaryToPlot[,independentVariable],main="Hypo Lesions", data = studySummaryToPlot, cex=2)
      graphics::boxplot(LESIONS_BOTH~studySummaryToPlot[,independentVariable], main="Both Type of Lesions", data = studySummaryToPlot, cex=2)
      graphics::boxplot(LESIONS_HYPER~studySummaryToPlot[,independentVariable],main="Hyper Lesions", data = studySummaryToPlot, cex=2)
      grDevices::dev.off()
    }

    #browser()
    if(depthAnalysis >1)
    {

      keys <- annotation_keys()
      # #browser()
      # areas <- rbind(genes, islands, dmrs)
      nkeys <- dim(keys)[1]
      # #browser()

      parallel::clusterExport(envir=environment(), cl = computationCluster, varlist =c("apply_model"))
      result_temp_foreach <- foreach::foreach(i = 1:nkeys, .combine = rbind) %dopar%
      # for (i in 1:nkeys)
      {
        # i <- 25
        if(exists("tempDataFrame"))
          rm(list = c("tempDataFrame"))
        key <- keys [i,]
        # print(key)
        fname <-file_path_build(resultFolderData,c(key$ANOMALY, key$FIGURE, key$GROUP,key$SUBGROUP),"csv")
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
            tempDataFrame$POPULATION <- sample_names[, independentVariable]
            tempDataFrame[is.na(tempDataFrame)] <- 0
            #  we want to preserve the NA in the indipendent variables to be removed by the models
            tempDataFrame[, independentVariable] <- sample_names[, independentVariable]
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

            result_temp <- apply_model(tempDataFrame = tempDataFrame, g_start = g_start, family = family_test, covariates = covariates, key = key, transformation= transformation, dototal = TRUE, logFolder= logFolder, independentVariable, depthAnalysis)

            # #browser()
            # n_adj <- iters - g_start
            # result <- rbind(result, result_temp)
            result_temp
          }
        }
      }
    }

    if(exists("result_temp_foreach"))
      result <- rbind(result, result_temp_foreach)

    # browser()
    result <- unique(result)
    # result[ result$AREA_OF_TEST=="TOTAL" & result$ANOMALY=="LESIONS" ,"PVALUEADJ" ] <- round(stats::p.adjust(result[ result$AREA_OF_TEST=="TOTAL" & result$ANOMALY=="LESIONS","PVALUE" ],method = "BH"),3)
    # result[ result$AREA_OF_TEST=="TOTAL" & result$ANOMALY=="MUTATIONS" ,"PVALUEADJ" ] <- round(stats::p.adjust(result[ result$AREA_OF_TEST=="TOTAL" & result$ANOMALY=="MUTATIONS","PVALUE" ],method = "BH"),3)

    result <- result[order(result$PVALUEADJ),]

    if(filterpvalue)
        result <- subset(result, PVALUE < 0.05 | PVALUEADJ < 0.05)

    if(nrow(result)>0)
    {
      if(is.null(covariates) || length(covariates)==0)
      {
        file_suffix <- "_test_result"
      } else
      {
        file_suffix <- "_test_corrected_result"
      }

      fileName <- file_path_build(resultFolderDataInference,c(file_result_prefix , transformation, family_test, file_suffix),"csv")
      utils::write.csv2(result,fileName , row.names = FALSE)
    }


    # res.pvalue <- subset(result, PVALUE < 0.05)
    # res.pvalue$beta_gt1 <- res.pvalue$BETA>1
    # res.pvalue$beta_gt1 <- as.numeric(res.pvalue$beta_gt1)

    # source("/home/lcorsaro/Desktop/Progetti/r-studio/smarties/R/microarray/epigenetics/epimutation_analysis/qqplot_inferential.R")
    # result <- utils::read.csv(file.path(resultFolderDataInference,paste0(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")))
    # qqunif.plot(diffMethTable_site_cmp1$diffmeth.p.val, resultFolderDataInference =  report.dir, filePrefix ="diff_meth_sites")

    # case_vs_control_binomial_regression_corrected_result <-
    #   utils::read.csv2(
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
    #               resultFolderDataInference = resultFolderDataInference,
    #               filePrefix = file_prefix)
    # }
    #
    #
    #
    # qqunif.plot(pvalues)
  }
}
