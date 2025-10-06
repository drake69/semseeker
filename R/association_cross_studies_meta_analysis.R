association_cross_studies_meta_analysis <- function(inference_details,statistic_parameter="BETA", pvalue_column="PVALUE_ADJ_ALL_BH",studies,
  studies_base_folder, result_folder)
{

  library(meta)
  library(tidyverse)
  keys <- na.omit(unique(results_inference[,c("FIGURE","SUBAREA")]))
  for (k in 1:nrow(keys))
  {
    # k <- 1
    results_inference_for <- subset(results_inference, FIGURE==keys[k,"FIGURE"] & SUBAREA==keys[k,"SUBAREA"])
    results_inference_for <- na.omit(results_inference_for[,c("BETA","STD.ERROR","STUDY",pvalue_column,"AREA_OF_TEST")])
    areas <- na.omit(unique(results_inference[, "AREA_OF_TEST"]))
    for (g in seq_along(areas))
    {
      # g <- 1
      first_area <- areas[g]
      if(is.na(first_area))

      results_inference_subset <- subset(results_inference_for, AREA_OF_TEST==first_area)
      # if( is.null(results_inference_subset[,"STUDY"]))
      #   next
      studies_count <- length(unique(results_inference_subset[,"STUDY"]))
      if(studies_count<2)
        next
      meta_model <- metagen(
        TE = results_inference_subset[,statistic_parameter],   # l'effetto di ogni studio
        seTE = results_inference_subset$STD.ERROR,  # lo standard error dell'effetto di ogni studio
        studlab = results_inference_subset$STUDY,  # il nome di ogni studio
        data = results_inference_subset,  # il dataframe che contiene tutti i dati
        comb.fixed = FALSE, # utilizza il modello a effetti casuali
        hakn = FALSE, # non applica la correzione di Hartung-Knapp
        TE.targ = 1, # specifica che la stima dell'effetto è il coefficiente di regressione
        pval = results_inference_subset[,pvalue_column],  # specifica la colonna del p-value nel dataframe
        ncpus = 9
      )
      meta_analysis_results <- summary(meta_model)
      # pval <- meta_analysis_results$pval
      beta <- meta_analysis_results$random$TE

      common.ci.lower <- meta_analysis_results$common$lower
      common.ci.upper <- meta_analysis_results$common$upper
      common.pval <- meta_analysis_results$pval.common

      random.ci.lower <- meta_analysis_results$random$lower
      random.ci.upper <- meta_analysis_results$random$upper
      random.pval <- meta_analysis_results$pval.random

      #test of heteroareaity
      pval.Q <- meta_analysis_results$pval.Q
      pval.fixed <- meta_analysis_results$pval.fixed
      k.study <- meta_analysis_results$k.study
      tau2 <- meta_analysis_results$tau2

      studies <- paste(sort(meta_analysis_results$studlab), collapse = " ")
      result <- data.frame(tau2,beta,common.ci.lower,common.ci.upper, common.pval, random.ci.lower,random.ci.upper, random.pval, pval.Q, pval.fixed,k.study, studies, first_area, studies_count)
      if (exists("final_result"))
        final_result <- rbind(final_result, result)
      else
        final_result <- result
    }
  }
}
