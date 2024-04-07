# markers <- c("MUTATIONS","LESIONS","DELTAS","DELTAQ","DELTAR","DELTARQ","BETA")
# disease <- "TCDD"
# study ="TCDD"
# adjust_per_area <- F
# adjust_globally <- F
# pvalue_columns <- c("PVALUE","PVALUEADJ","PVALUEADJ_ALL_BH")
# pvalue_limit <- c(0.05,0.05,0.05)
# # pvalue_columns <- c("PVALUEADJ_ALL_BH")
# # pvalue_limit <- c(0.05)
# adjustment_method <- "BH"
# path_db <- c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF", "cell_markers", "mmu_KEGG")
# phenolyzer_folder_bin <- "/usr/local/lib/phenolyzer/"
#
# for (id in 1:nrow(inference_details))
# {
#   inference_detail <- inference_details[id,]
#   for (i in 1:length(pvalue_columns)){
#     # i <- 1
#     pvalue_column <- pvalue_columns[i]
#     pvalue <- pvalue_limit[i]
#
#     # score_limit <- c(0.2,1)
#     # disgenetplus2r_user="lcorsaro69@gmail.com"
#     # disgenetplus2r_password="sitnuq-qussa5-vytpoR"
#     # phenotype_disgenetplus2r(study= study, disease=disease , score_limit = score_limit,disgenetplus2r_user, disgenetplus2r_password,
#     #   pvalue = pvalue, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column = pvalue_column,
#     #   inference_details= inference_detail, result_folder=result_folder, maxResources = 90, parallel_strategy  = "multicore", showprogress=TRUE)
#
#     pathway_pathfindR(study,
#       path_db,  iterations = 10, exclude_beta= FALSE,
#       pvalue = pvalue, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column= pvalue_column,
#       inference_details=inference_detail,result_folder=result_folder, maxResources = 90, parallel_strategy  = "multicore", showprogress=TRUE, markers = markers)
#
#     pathway_WebGestalt( study = study,
#       types=c("BP","MF"),  enrich_methods = c("ORA"),
#       pvalue = pvalue, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column=pvalue_column,
#       inference_details = inference_detail,result_folder, maxResources = 90, parallel_strategy  = "multicore", showprogress=TRUE)
#
#     phenotype_phenolyzer <- phenotype_phenolyzer(study,
#       disease,phenolyzer_folder_bin,minimum_score = 0.5,
#       pvalue =pvalue, adjust_per_area = F, adjust_globally = F,adjustment_method = "BH", pvalue_column=pvalue_column,
#       inference_details=inference_detail,result_folder, maxResources = 90, parallel_strategy  = "multicore", showprogress=TRUE)
#
#   }
# }
