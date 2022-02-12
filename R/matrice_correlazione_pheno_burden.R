# library(readr)
# sample_sheet_result <- read_delim("~/Desktop/experiments_data/DIOSSINA_DESIO/20211018/sample_sheet_result.csv",delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = "."), trim_ws = TRUE)
# sample_sheet_result$Sample_ID <- paste0(sample_sheet_result$Sentrix_Position,"_", sample_sheet_result$Sentrix_ID, sep="")
#
# sample_sheet_result_pheno <- utils::read.csv("Desktop/experiments_data/DIOSSINA_DESIO/Samplesheet_Completed_Diossina_With_Phenotype.csv")
# sample_sheet_result_pheno$Sample_ID <- paste0(sample_sheet_result_pheno$Sentrix_Position,"_", sample_sheet_result_pheno$Sentrix_ID, sep="")
#
# sample_sheet_result <- merge(sample_sheet_result, sample_sheet_result_pheno, by.x="Sample_ID", by.y="Sample_ID", all.x = T)
# sample_sheet_result <- as.data.frame(sample_sheet_result)
# sample_sheet_result$mother_exposition_difference <- sample_sheet_result$mother_age_delivery - sample_sheet_result$mother_exposure_age
# sample_sheet_result$Sample_Group <- sample_sheet_result$Sample_Group.x
#
# sample_sheet_result[is.na(sample_sheet_result$exam_age),"exam_age"] <-18
#
# sperm_count_limit <- c(15,200)
# sperm_count_integration <- runif(7,15,200)
#
# write_csv(sample_sheet_result,"~/Desktop/experiments_data/DIOSSINA_DESIO/20211018/sample_sheet_result_complete.csv")
#
# sample_sheet_result[is.na(sample_sheet_result$sperm_count),"sperm_count"] <- sperm_count_integration
# colToCorrelate <- c("exam_age","sperm_count","tcdd_mother","mother_exposure_age","mother_age_delivery","mother_exposition_difference")
#
# library(readxl)
# # library(Hmisc)
#
#
# # bootstrap cluster
# logFolder <- paste0(tempdir(),"/log")
#
# if (logFolder != "" && !dir.exists(logFolder)) {
#   dir.create(logFolder)
# }
# nCore <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
# outFile <- paste0(logFolder, "/cluster_r.out", sep = "")
# computationCluster  <-  parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1, type = "PSOCK", outfile = outFile)
# doParallel::registerDoParallel(computationCluster)
#
#
# subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
# mainGroupLabel =  "GENE"
# figures <- c("HYPO", "HYPER", "BOTH")
# anomalies <- c("MUTATIONS","LESIONS")
# keys <- expand.grid("anomalies"= anomalies, "figures"=figures, "group"="GENE", "subGroup"=subGroups)
#
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# library(foreach)
# parallel::clusterExport(envir=environment(), cl = computationCluster, varlist = list(  "read_excel"))
#
# mmm <- data.frame("var1"="",
#                      "var2"="",
#                      "p.value.spearman"="",
#                      "p.value.kendall"="",
#                      "p.value.pearson"="",
#                      "cor.spearman"="",
#                      "cor.kendall"="",
#                      "cor.pearson"="",
#                      "stats::mean"="",
#                      "sd"="",
#                      "key"="",
#                      "p.adj.area.sperman"="",
#                      "p.adj.area.kendall"="",
#                      "p.adj.area.pearson"="",
#                      "shapiro.pvalue"="",
#                      "bartlett.pvalue"="")
# mmm <- mmm[-1,]
#
# for(j in c("GENE","ISLAND","DMR"))
# {
#   # j <- "DMR"
#   message(j)
#   keys <- excel_sheets(paste0("~/Desktop/experiments_data/DIOSSINA_DESIO/20211018/Pivots/",j,".xlsx", sep=""))
#   for(i in 1:length(keys))
#   # sheetList <- foreach::foreach(i=1:length(keys), .combine= rbind ) %dopar%
#   {
#     # i <-1
#     # key <- paste0(keys[i,"anomalies"],"_",keys[i,"figures"],"_",keys[i,"group"],"_",keys[i,"subGroup"] ,sep="")
#     key <- keys[i]
#     # key <- "LESIONS_HYPER_GENE_ExonBnd"
#     message(key)
#     tempDataFrame <- read_excel(paste0("~/Desktop/experiments_data/DIOSSINA_DESIO/20211018/Pivots/",j,".xlsx", sep=""), sheet = key)
#     tempDataFrame <- as.data.frame(t(tempDataFrame))
#     colnames(tempDataFrame) <- tempDataFrame[1,]
#     tempDataFrame$Sample_ID <- rownames(tempDataFrame)
#     tempDataFrame <- tempDataFrame[-1,]
#
#     tempDataFrame <- merge(sample_sheet_result[,c("Sample_Group","Sample_ID", colToCorrelate)] , tempDataFrame, by.x="Sample_ID",by.y="Sample_ID", all.x=T)
#     tempDataFrame <-tempDataFrame[!is.na(tempDataFrame$sperm_count),]
#
#     areaBurden <- colnames(tempDataFrame)[!(colnames(tempDataFrame) %in% c("Sample_Group","Sample_ID","POPULATION", colToCorrelate))]
#
#     tt <- expand.grid(areaBurden, colToCorrelate)
#
#     hh <- foreach::foreach(g=1:nrow(tt), .combine=rbind) %dopar%
#       {
#         # g <- 2
#         c1 <- as.character(tt[g,1])
#         var1 <- as.numeric(tempDataFrame[,c1])
#         var1[is.na(var1)] <-0
#         c2 <- as.character(tt[g,2])
#         var2 <- as.numeric(tempDataFrame[,c2])
#         res.spearman <- cor.test(var1,var2, method="spearman")
#         res.kendall <- cor.test(var1,var2, method="kendall")
#
#         bartlett.p.value <- NA
#         try (bartlett.p.value  <- bartlett.test(var1,tempDataFrame[,"Sample_Group"])$p.value)
#         shapiro.p.value <- NA
#         try (shapiro.p.value <- shapiro.test(var1)$p.value)
#
#         res.pearson.p.value <- NA
#         res.pearson.estimate <- NA
#         try
#         {
#           res.pearson <- cor.test(var1,var2, method="pearson")
#           res.pearson.p.value <- res.pearson$p.value
#           res.pearson.estimate <- res.pearson$estimate
#         }
#
#         data.frame("var1"=c1,
#                    "var2"=c2,
#                    "p.value.spearman"=res.spearman$p.value,
#                    "p.value.kendall"=res.kendall$p.value,
#                    "p.value.pearson"=res.pearson.p.value,
#                    "cor.spearman"=res.spearman$estimate,
#                    "cor.kendall"=res.kendall$estimate,
#                    "cor.pearson"=res.pearson.estimate,
#                    "stats::mean"=stats::mean(var1),
#                    "sd"=stats::sd(var1),
#                    "key"=key,
#                    "p.adj.area.spearman"="",
#                    "p.adj.area.kendall"="",
#                    "p.adj.area.pearson"="",
#                    "shapiro.pvalue"=shapiro.p.value,
#                    "bartlett.pvalue"=bartlett.p.value
#         )
#       }
#
#     hh$p.adj.area.kendall <- stats::p.adjust(hh$p.value.kendall,"BH")
#     hh$p.adj.area.spearman <- stats::p.adjust(hh$p.value.spearman,"BH")
#     hh$p.adj.area.pearson <- stats::p.adjust(hh$p.value.pearson,"BH")
#     mmm <- rbind(mmm,hh)
#     #
#     # mcorr <- rcorr(as.matrix(tempDataFrame[, !(colnames(tempDataFrame) %in% c("Sample_ID","POPULATION")) ] ),type="pearson")
#     # rm(tempDataFrame)
#     # mmm <- flattenCorrMatrix(mcorr$r, mcorr$P)
#     # rm(mcorr)
#     gc()
#   }
# }
#
#
# mmm$p.adj.kendall <- stats::p.adjust(mmm$p.value.kendall,"BH")
# mmm$p.adj.spearman <- stats::p.adjust(mmm$p.value.spearman,"BH")
# mmm$p.adj.pearson <- stats::p.adjust(hh$p.value.pearson,"BH")
#
# write_csv(mmm,"~/Desktop/experiments_data/DIOSSINA_DESIO/20211018/Inference/correlation.csv")
#
# genomica_area_correlation <- subset(mmm, p.adj.kendall < 0.05)
#
# genomica_area_correlation <- subset(mmm, var2=="sperm_count")
# library(ggplot2)
# genomica_area_correlation_1 <- subset(genomica_area_correlation, key=="MUTATIONS_BOTH_GENE_Whole")
# p1 <- ggplot(data=genomica_area_correlation_1, aes(x=cor.kendall , y=-log10(p.value.kendall))) +
#   geom_point() + theme_minimal() + geom_vline(xintercept=c(0), col="blue") +  geom_hline(yintercept=-log10(0.05), col="blue") +
#   ggtitle("Mutations Gene Whole Both") + theme(plot.title = element_text(hjust = 0.5 ))
#
# genomica_area_correlation_2 <- subset(genomica_area_correlation, key=="MUTATIONS_BOTH_ISLAND_Whole")
# p2 <- ggplot(data=genomica_area_correlation_2, aes(x=cor.kendall , y=-log10(p.value.kendall))) +
#   geom_point() + theme_minimal() + geom_vline(xintercept=c(0), col="blue") +  geom_hline(yintercept=-log10(0.05), col="blue") +
#   ggtitle("Mutations Island Whole Both") + theme(plot.title = element_text(hjust = 0.5 ))
#
# genomica_area_correlation_3 <- subset(genomica_area_correlation, key=="MUTATIONS_BOTH_DMR_Whole")
# p3 <- ggplot(data=genomica_area_correlation_3, aes(x=cor.kendall , y=-log10(p.value.kendall))) +
#   geom_point() + theme_minimal() + geom_vline(xintercept=c(0), col="blue") +  geom_hline(yintercept=-log10(0.05), col="blue") +
#   ggtitle("Mutations DMR Whole Both") + theme(plot.title = element_text(hjust = 0.5 ))
#
# lay <- rbind(
#   c(1,1,2,2),
#   c(NA,3,3,NA)
# )
#
# filename = paste0( chartFolder,"/","scatter_plot_rnbeads.png",sep="")
# grDevices::png(file= filename, width=2480, height = 2480)
# gridExtra::grid.arrange(p1,p2,p3,layout_matrix = lay)
# grDevices::dev.off()
#
# parallel::stopCluster(computationCluster)
# gc()
#
#
