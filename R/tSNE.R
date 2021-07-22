#
# finalBed$GROUP <- finalBed$RELATION_TO_CPGISLAND
# finalBed$GENE <- finalBed$ISLAND
#
# bodyGene <- subset(finalBed, finalBed$GROUP=="Island")
#
# test1<-reshape2::dcast(data = bodyGene, POPULATION + GENE + SAMPLEID ~ FIGURE , value.var = "freq", sum)
# View(test1)
# plot(test1)
#
# utils::write.table(test1, file = "~/Desktop/test_sni.csv" , quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
#
# test2 <- reshape2::dcast(data = test1, POPULATION + SAMPLEID + HYPO + HYPER ~ GENE  , value.var = c("GENE"), fun.aggregate = length)
# View(test2)
#
# test4 <- data.frame(test2[,1:4], "REGIONS_COUNT" =rowSums(test2[,-c(1:4)]))
# View(test4)
#
#
# test3 <- reshape2::dcast(data = test4, REGIONS_COUNT + POPULATION + HYPO + HYPER ~ SAMPLEID  , value.var = c("SAMPLEID"), fun.aggregate = length)
# View(test3)
#
# test5 <- data.frame(test3[,1:4], "SAMPLES_COUNT" =rowSums(test3[,-c(1:4)]))
# View(test5)
#
# ggplot(test5, aes(HYPO, HYPER,REGIONS_COUNT,  colour = POPULATION)) + geom_point()
#
# library(plotly)
# fig <- plot_ly(x=test5$HYPO, y=test5$HYPER, z=test5$REGIONS_COUNT, type="scatter3d", mode="markers", size=test5$SAMPLES_COUNT , color= test5$POPULATION )
# fig <- plot_ly(x=test5$HYPO, y=test5$HYPER, z=test5$REGIONS_COUNT, type="scatter3d", mode="markers", color= 1/test5$SAMPLES_COUNT )
# # fig <- plot_ly(x=test5$HYPO, y=test5$HYPER, z=test5$REGIONS_COUNT, type="scatter3d", mode="markers", color= 1/test5$SAMPLES_COUNT alpha=test5$REGIONS_COUNT)
# # fig <- layout(fig, xaxis = list(type = "log"), yaxis = list(type = "log"))
# fig
#
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# #
# # # The following initializes usage of Bioc devel
# # BiocManager::install(version='devel')
# #
# # BiocManager::install("M3C")
# #
# # library(M3C)
# # tsne(mydata= test1 ,colvec=c('gold'), perplex = 1,seed = 100)
#
#
# library(Rtsne)
# # normalize_input(test1)
#
# iris_unique <- unique(test1) # Remove duplicates
# # iris_unique
# # iris_matrix <- as.matrix(iris_unique[,4:5])
# # X <- normalize_input(iris_matrix)
# # colMeans(X)
# # range(X)
# # Rtsne(iris_unique, dims = 2, initial_dims = 50,
# #       perplexity = 30, theta = 0.5, check_duplicates = TRUE,
# #       pca = TRUE, partial_pca = FALSE, max_iter = 1000,
# #       verbose = getOption("verbose", FALSE), is_distance = FALSE, Y_init = NULL, pca_center = TRUE, pca_scale = FALSE,
# #       normalize = TRUE, stop_lying_iter = ifelse(is.null(Y_init), 250L, 0L), mom_switch_iter = ifelse(is.null(Y_init), 250L, 0L), momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 7)
#
#
# library(Rtsne)
# inputBedDataFrame <- subset(geneBed, POPULATION != "Reference")
# inputBedDataFrame <- data.frame(inputBedDataFrame,"KEY" = paste(inputBedDataFrame[, 1],"_",inputBedDataFrame[, 3],"_",inputBedDataFrame$FIGURE,sep=""))
# inputBedDataFrame$KEY <- as.factor(inputBedDataFrame$KEY)
#
# tempDataFrame <- subset(inputBedDataFrame, ANOMALY == "LESIONS")
# tempDataFrame <- reshape2::dcast(data = tempDataFrame, KEY ~ SAMPLEID, value.var = "freq", sum)
# rownames(tempDataFrame) <- tempDataFrame$KEY
# stats::heatmap(as.matrix(tempDataFrame[,2:dim(tempDataFrame)[2]]),scale = "column", margins = c(5, 5))

# mynumbers <- Rtsne(tempDataFrame)
#
# plot(mynumbers$Y,col=c("Yellow", "Blue", "Red"), asp=1)
#
# plot(mynumbers$Y, asp = 1, pch = 20, col = "red",
#      cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5,
#      xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2",
#      main = "2D t-SNE projection")
#
# mynumbers <- Rtsne(tempDataFrame, dims = 3, initial_dims = 50,
#       perplexity = 30, theta = 0.5, check_duplicates = TRUE,
#       pca = TRUE, partial_pca = FALSE, max_iter = 10,
#       verbose = getOption("verbose", FALSE), is_distance = FALSE,  pca_center = TRUE, pca_scale = FALSE,
#       normalize = TRUE, momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 7)

# #levels(as.factor(iris_unique$POPULATION))
# plot(mynumbers$Y,col=c("Yellow", "Blue", "Red"), asp=1)
# #plot(mynumbers$Y)
#
#
#
#
#

#
# library(M3C)
#
# pca(tempData,legendtextsize = 10,axistextsize = 10,dotsize=2)
#
# tsne(tempDataFrame)
# tsne(tempData)
# browseVignettes("M3C")


