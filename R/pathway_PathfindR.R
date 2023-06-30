# pathway_pathfindR <- function(resultFolder, inferenceFile, pvalue, path_db) {
#   prefix <- stringi::stri_rand_strings(1, 7, pattern = "[A-Za-z0-9]")
#   # "~/Desktop/What_s_going_On/papers_writing/breast_study/TCGA-BRCA/analysis/white_tumor/epimutation/Inference/3_er_status_by_ihc_none_gaussian_test_corrected_result.csv"
#   results_inference <- read.csv(paste(resultFolder, "/Inference/", inferenceFile, sep = ""), sep = ";", dec = ",")
#   results_inference <- subset(results_inference, results_inference$PVALUEADJ_ALL_BY < pvalue)
#   # features <- unique(results_inference[, c("AREA_OF_TEST","AREA")])
#   # table(features$AREA)

#   # results_inference <- read.csv("~/Desktop/What_s_going_On/papers_writing/breast_study/TCGA-BRCA/analysis/white_tumor/epimutation/Inference/gene_3_er_status_by_ihc_none_gaussian_test_corrected_result.csv", sep=";", dec=",")
#   # results_inference <- results_inference[results_inference$PVALUEADJ_ALL_BH < pvalue,]
#   results_inference <- subset(results_inference, AREA == "GENE" & AREA_OF_TEST != "TOTAL")
#   # results_inference$PVALUEADJ_ALL_BH <- p.adjust(results_inference$PVALUEADJ)
#   keys <- unique(results_inference[, c("SUBAREA", "AREA", "MARKER", "FIGURE")])
#   seq <- 0
#   browser()
#   for (k in 1:length(path_db))
#   {
#     for (i in 1:nrow(keys))
#     {
#       # k <- 1
#       # i <- 10
#       gene_set <- results_inference[
#         results_inference$SUBAREA == keys[i, ]$SUBAREA &
#           results_inference$AREA == keys[i, ]$AREA &
#           results_inference$MARKER == keys[i, ]$MARKER &
#           results_inference$FIGURE == keys[i, ]$FIGURE,
#         c("AREA_OF_TEST", "BETA", "PVALUEADJ_ALL_BH"),
#       ]
#       # table(is.na(results_inference))
#       gene_set <- na.omit(gene_set)

#       applied_model <- paste(c(unique(results_inference$INDIPENDENT.VARIABLE), unlist(strsplit(unique(results_inference$COVARIATES), " "))), collapse = "_")
#       try({
#         output_temp <- pathfindR::run_pathfindR(gene_set, path_db[k],
#           max_gset_size = nrow(gene_set), iterations = 10,
#           output_dir = paste(paste(resultFolder, "/Pathway/", prefix, "_finder_", sep = ""), seq, sep = "")
#         )
#         seq <- seq + 1
#         if (nrow(output_temp) == 0) {
#           next()
#         }
#         output_temp$key <- paste(keys[i, ]$FIGURE, keys[i, ]$MARKER, keys[i, ]$AREA, keys[i, ]$SUBAREA, sep = "_")
#         output_temp$seq <- seq
#         output_temp$gene_count <- nrow(gene_set)
#         output_temp$source <- path_db[k]
#         output_temp$order <- 1:nrow(output_temp)
#         if (exists("pathway_result")) {
#           pathway_result <- rbind(pathway_result, output_temp)
#         } else {
#           pathway_result <- output_temp
#         }
#         write.csv(pathway_result, paste(resultFolder, "/Pathway/", prefix, "_", gsub("[.]", "", pvalue), "_", applied_model, "_pathway_result.csv", sep = ""))
#       })
#     }
#   }
# }
