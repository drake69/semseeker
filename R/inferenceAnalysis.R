inferenceAnalysis <- function(studySummary, resultFolder, logFolder, covariates, family, transformation = NULL)
{

  # browser()

  workingFolder <- file.path(resultFolder, "Pivots/")
  inferenceFolder <-  file.path(resultFolder, "Inference/")
  file_result_prefix <- "case_vs_control_"

  studySummary <- subset(studySummary, Sample_Group!="Reference")
  #########################################################################################################
  #########################################################################################################

  if(!dir.exists(inferenceFolder))
    dir.create(inferenceFolder)

  if(is.null(covariates) || length(covariates)==0)
  {
    sample_names <- data.frame(studySummary[, c("Sample_ID", "Sample_Group")])
  } else
  {
    sample_names <- data.frame(studySummary[, c("Sample_ID", "Sample_Group", covariates)])
  }
  # sample_names <- data.frame(studySummary)
  result = data.frame (
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
    "MEAN.CASE" ="",
    "MEAN.CONTROL"="",
    "SD.CASE"="",
    "SD.CONTROL"=""
  )

  result <- result[-1,]
  # sample_names[, c("Sample_Group")]  <-  toupper(sample_names[, c("Sample_Group")])
  # sample_names <- sample_names[complete.cases(sample_names),]
  sample_names <- subset(sample_names, sample_names[, "Sample_Group"]  != "Reference")

  if(!is.null(covariates) && !length(covariates)==0)
  {
    sample_names <-   sample_names[, c("Sample_ID", "Sample_Group", covariates)]
    colnames(sample_names) <- c("Sample_ID", "Sample_Group", covariates)
  } else
  {
    sample_names <-   sample_names[, c("Sample_ID", "Sample_Group")]
    colnames(sample_names) <- c("Sample_ID", "Sample_Group")
  }
  ######################################################################################################
  # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
  # binomiale che si vuole usare per la regressione logistica



  keys <- expand.grid("ANOMALY"=c("LESIONS","MUTATIONS"), "FIGURE"=c("HYPER","HYPO","BOTH"))
  cols <- paste(keys$ANOMALY,"_",keys$FIGURE, sep="")
  keys$GROUP = "POPULATION"
  keys$SUBGROUP = "SAMPLE"
  iters <- length(cols)

  if(!is.null(covariates) && !length(covariates)==0)
    studySummary <- studySummary[, c("Sample_Group", covariates, cols) ]

  # browser()

  for(i in 1:nrow(keys))
  {
    # i <- 1
    # family <- "poisson"
    # transformation <- NULL
    g_start <- 2 + length(covariates)
    result_temp <- apply_model(df = studySummary[, c("Sample_Group", covariates, cols[i])], g_start = g_start , family = family, covariates = covariates, key = keys[i,], transformation = transformation, dototal = FALSE, logFolder= logFolder)
    result <- rbind(result, result_temp)
  }

  result <- result[order(result$PVALUEADJ),]

  # browser()

  studySummaryToPlot <- studySummary
  studySummaryToPlot$Sample_Group  <- relevel(as.factor(studySummaryToPlot$Sample_Group), "Control")

  chartFolder <- paste(inferenceFolder, "/Charts/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  chartFolder <- paste(inferenceFolder, "/Charts/POPULATION_COMPARISON/", sep="")
  if (chartFolder != "" && !dir.exists(chartFolder)) {
    dir.create(chartFolder)
  }

  filename = paste( chartFolder,"/MUTATIONS.png",sep="")
  grDevices::png(file= filename, width=2000, height=2000)
  par(mfrow=c(1,3))
  graphics::boxplot(MUTATIONS_HYPO~Sample_Group,main="Hypo Mutations", data = studySummaryToPlot, cex=2)
  graphics::boxplot(MUTATIONS_BOTH~Sample_Group, main="Both Type of Mutations", data = studySummaryToPlot, cex=2)
  graphics::boxplot(MUTATIONS_HYPER~Sample_Group,main="Hyper Mutations", data = studySummaryToPlot, cex=2)
  grDevices::dev.off()

  filename = paste( chartFolder,"/LESIONS.png",sep="")
  grDevices::png(file= filename, width=2000, height=2000)
  par(mfrow=c(1,3))
  graphics::boxplot(LESIONS_HYPO~Sample_Group,main="Hypo Lesions", data = studySummaryToPlot, cex=2)
  graphics::boxplot(LESIONS_BOTH~Sample_Group, main="Both Type of Lesions", data = studySummaryToPlot, cex=2)
  graphics::boxplot(LESIONS_HYPER~Sample_Group,main="Hyper Lesions", data = studySummaryToPlot, cex=2)
  grDevices::dev.off()


  anomalies <- c("MUTATIONS", "LESIONS")
  figures <- c("HYPO", "HYPER", "BOTH")

  group =  "GENE"
  subGroups <-
    c("Body",
      "TSS1500",
      "5'UTR",
      "TSS200",
      "1stExon",
      "3'UTR",
      "ExonBnd",
      "Whole")

  # genes <- data.frame("group" = group, "area" = subGroups)

  keys <-
    expand.grid(
      "ANOMALY" = anomalies,
      "FIGURE" = figures,
      "GROUP" = group,
      "SUBGROUP" = subGroups
    )

  group <- "ISLAND"
  subGroups <- c("N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "Island", "Whole")

  # islands <- data.frame("group" = group, "area" = subGroups)

  keys <-
    rbind(
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      ),
      keys
    )


  group =  "DMR"
  subGroups <- c("DMR")

  # dmrs <- data.frame("group" = group, "area" = subGroups)

  keys <-
    rbind(
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      ),
      keys
    )

  # browser()
  # areas <- rbind(genes, islands, dmrs)
  nkeys <- dim(keys)[1]
  for (i in 1:nkeys)
  {
    # i <- 25
    rm(list = c("df"))
    key <- keys [i,]
    # print(key)
    fname <-file.path(workingFolder,paste("/", key$ANOMALY, "_", key$FIGURE, "_", key$GROUP, "_", key$SUBGROUP, ".csv" , sep =""))
    print(fname)
    if (!file.exists(fname))
    {
      print("File Not found!")
      next
    }
    df <- read.csv(fname, sep = ";")
    row.names(df) <- df$X
    df <- df[,-1]
    df <- t(df)
    if(dim(df)[1]<=1)
      next
    # df <- df[,-2]
    df <- as.data.frame(df)
    df <- subset(df, POPULATION != "Reference")
    df <- subset(df, POPULATION != 0)

    df <-  merge( x =  sample_names, y =  df,  by.x = "Sample_ID",  by.y = "SAMPLEID" , all.x = TRUE)
    df$POPULATION <- sample_names$Sample_Group
    df[is.na(df)] <- 0
    df <- df[, !(names(df) %in% c("POPULATION","Sample_ID"))]

    cols <- (gsub(" ", "_", colnames(df[])))
    cols <- (gsub("-", "_", cols))
    cols <- (gsub(":", "_", cols))
    cols <- (gsub("/", "_", cols))
    cols <- (gsub("'", "_", cols))

    colnames(df) <- cols

    df <- as.data.frame(df)
    # df$Sample_Group <- as.factor(df$Sample_Group)
    # df$Sample_Group  <- relevel(df$Sample_Group, "Control")

    g_start <- 2 + length(covariates)


    result_temp <- apply_model(df = df, g_start = g_start, family = family, covariates = covariates, key = key, transformation= transformation, dototal = TRUE, logFolder= logFolder)

    # n_adj <- iters - g_start
    result <- rbind(result, result_temp)
  }


  result[ result$AREA_OF_TEST=="TOTAL","PVALUEADJ" ] <- round(p.adjust(result[ result$AREA_OF_TEST=="TOTAL","PVALUE" ],method = "BH"),3)
  result <- result[order(result$PVALUEADJ),]

  if(is.null(covariates) || length(covariates)==0)
  {
    file_suffix <- "_regression_result.csv"
  } else
  {
    file_suffix <- "_regression_corrected_result.csv"
  }

  fileName <- paste(file_result_prefix , transformation, "_", family, file_suffix, sep = "")
  write.csv2(result, file.path(
    inferenceFolder,
    fileName
  ), row.names = FALSE)


  return(fileName)
  # res.pvalue <- subset(result, PVALUE < 0.05)
  # res.pvalue$beta_gt1 <- res.pvalue$BETA>1
  # res.pvalue$beta_gt1 <- as.numeric(res.pvalue$beta_gt1)

  # source("/home/lcorsaro/Desktop/Progetti/r-studio/smarties/R/microarray/epigenetics/epimutation_analysis/qqplot_inferential.R")
  # result <- read_csv(file.path(inferenceFolder,paste(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")))
  # qqunif.plot(diffMethTable_site_cmp1$diffmeth.p.val, inferenceFolder =  report.dir, filePrefix ="diff_meth_sites")

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
  #     next
  #   #######inserisco nella funzione i pvalues ottenuti dalla differential (non aggiustati)
  #
  #   file_prefix <- paste0("case_vs_control_binomial_regression_corrected_result","_", key$ANOMALY,"_", key$FIGURE,"_", key$GROUP,"_", key$SUBGROUP,"_", sep="")
  #   qqunifPlot(diffmeth.p.val$PVALUE,
  #               inferenceFolder = inferenceFolder,
  #               filePrefix = file_prefix)
  # }
  #
  #
  #
  # qqunif.plot(pvalues)

}
