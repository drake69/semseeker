inferenceAnalysisWithCorrection <- function(sampleSheet, resultFolder, logFolder, covariates)
{

  # browser()

  nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
  outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
  print(outFile)
  computation_cluster <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
  doParallel::registerDoParallel(computation_cluster)

  # options(digits = 22)
  parallel::clusterExport(envir=environment(), cl = computation_cluster, varlist =c())

  workingFolder <- file.path(resultFolder, "Pivots/")
  resultFolder <-  file.path(resultFolder, "Inference/")
  file_result_prefix <- "case_vs_control_"

  sampleSheet <- subset(sampleSheet, Sample_Group!="Reference")
  #########################################################################################################
  #########################################################################################################

  if(!dir.exists(resultFolder))
    dir.create(resultFolder)

  # sampleSheet$Sample_Group <- as.numeric(sampleSheet[,"Sample_Group"] == "Case")
  sample_names <- data.frame(sampleSheet[, c("Sample_ID", "Sample_Group", covariates)])
  # sample_names <- data.frame(sampleSheet)
  result = data.frame (
    "ANOMALY" = "",
    "FIGURE" = "",
    "GROUP" = "",
    "SUBGROUP" = "",
    "AREA_OF_TEST" = "",
    "PVALUE" = "",
    "PVALUEADJ" = "",
    "TEST" = "",
    "BETA" = ""
  )

  result <- result[-1,]
  # sample_names[, c("Sample_Group")]  <-  toupper(sample_names[, c("Sample_Group")])
  # sample_names <- sample_names[complete.cases(sample_names),]
  sample_names <- subset(sample_names, sample_names[, "Sample_Group"]  != "Reference")

  sample_names <-   sample_names[, c("Sample_ID", "Sample_Group", covariates)]
  colnames(sample_names) <- c("Sample_ID", "Sample_Group", covariates)
  ######################################################################################################
  # sample_names deve avere due colonne la prima con il nome del campione e la seconda con la variabile categorica
  # binomiale che si vuole usare per la regressione logistica

  anomalies <- c("MUTATIONS", "LESIONS")
  figures <- c("HYPO", "HYPER")

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
    rm(df)
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
    df$Sample_Group[df$Sample_Group=="Control"] <- 0
    df$Sample_Group[df$Sample_Group=="Case"] <- 1
    df$Sample_Group <- as.numeric(df$Sample_Group)

    g_start <- dim(sample_names)[2] + 1
    df[,g_start:ncol(df)] = as.data.frame(sapply(df[, g_start:ncol(df)] , as.numeric))

    iters <- length(cols)

    result_temp <- foreach::foreach(g = g_start:iters, .combine = rbind) %dopar%
      # for (g in g_start:iters)
      {
        #g <- 2
        covariates_model <- paste(paste0(c(cols[g], covariates),collapse="+", sep=""))
        sig.formula <- as.formula(paste0("Sample_Group","~", covariates_model, sep=""))
        result.glm  <- glm( sig.formula, family = "binomial", data = df)
        pvalue <- summary( result.glm )$coeff[-1, 4][1]
        pvalueadjusted <- pvalue
        beta <- exp(summary( result.glm )$coeff[-1, 1][1])
        data.frame (
          "ANOMALY" = key$ANOMALY,
          "FIGURE" = key$FIGURE,
          "GROUP" = key$GROUP,
          "SUBGROUP" = key$SUBGROUP,
          "AREA_OF_TEST" = cols[g],
          "PVALUE" = pvalue,
          "PVALUEADJ" = pvalueadjusted,
          "TEST" = "SINGLE_AREA",
          "BETA" = beta
        )
      }
    result_temp$PVALUEADJ <- p.adjust(result_temp$PVALUE,method = "BH", n = dim(result_temp)[1])
    result <- rbind(result, result_temp)
  }

  write.csv2(result, file.path(
    resultFolder,
    paste(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")
  ), row.names = FALSE)

  case_summary <- read.csv(file.path(workingFolder, "../Case/summary.csv"))
  ctrl_summary <- read.csv(file.path(workingFolder, "../Control/summary.csv"))
  # reference_summary <-  read.csv(file.path(workingFolder, "../Reference/summary.csv"))
  # summary <- rbind(case_summary, ctrl_summary, reference_summary)

  summary <- rbind(case_summary, ctrl_summary)

  summary$Mutations_Total <- summary$Hyper_Mutations + summary$Hypo_Mutations
  summary$Lesions_Total <- summary$Hyper_Lesions + summary$Hypo_Lesions

  summary$Sample_Group <- as.factor(summary$Sample_Group)
  summary$Sample_Group <- relevel(summary$Sample_Group, "Control")

  m1 <- glm(Hyper_Lesions ~ Sample_Group, family = "poisson", data = summary)
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "LESIONS",
    "FIGURE" = "HYPER",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)

  m1 <- glm(summary$Hypo_Lesions ~ summary$Sample_Group, family = "poisson")
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "LESIONS",
    "FIGURE" = "HYPO",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)

  m1 <- glm(summary$Hypo_Mutations ~ summary$Sample_Group, family = "poisson")
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "MUTATIONS",
    "FIGURE" = "HYPO",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)

  m1 <- glm(summary$Hyper_Mutations ~ summary$Sample_Group, family = "poisson")
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "MUTATIONS",
    "FIGURE" = "HYPER",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)

  m1 <- glm(summary$Mutations_Total ~ summary$Sample_Group, family = "poisson")
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "MUTATIONS",
    "FIGURE" = "TOTAL",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)

  m1 <- glm(summary$Lesions_Total ~ summary$Sample_Group, family = "poisson")
  pvalue <- summary(m1)$coeff[-1, 4]
  beta <- exp(summary(m1)$coeff[-1, 1])
  result_temp <- data.frame (
    "ANOMALY" = "LESIONS",
    "FIGURE" = "TOTAL",
    "GROUP" = "TOTAL",
    "SUBGROUP" = "TOTAL",
    "AREA_OF_TEST" = "SAMPLE",
    "PVALUE" = pvalue,
    "PVALUEADJ" = 0.05,
    "TEST" = "SINGLE_AREA",
    "BETA" = beta
  )
  result <- rbind(result, result_temp)


  result_temp$PVALUEADJ <- p.adjust(result_temp$PVALUE,method = "BH", n = dim(result_temp)[1])

  write.csv2(result, file.path(
    resultFolder,
    paste(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")
  ), row.names = FALSE)


  # res.pvalue <- subset(result, PVALUE < 0.05)
  # res.pvalue$beta_gt1 <- res.pvalue$BETA>1
  # res.pvalue$beta_gt1 <- as.numeric(res.pvalue$beta_gt1)

  # source("/home/lcorsaro/Desktop/Progetti/r-studio/smarties/R/microarray/epigenetics/epimutation_analysis/qqplot_inferential.R")
  # result <- read_csv(file.path(resultFolder,paste(file_result_prefix , "binomial_regression_corrected_result.csv", sep = "")))
  # qqunif.plot(diffMethTable_site_cmp1$diffmeth.p.val, resultFolder =  report.dir, filePrefix ="diff_meth_sites")

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
  #               resultFolder = resultFolder,
  #               filePrefix = file_prefix)
  # }
  #
  #
  #
  # qqunif.plot(pvalues)

}
