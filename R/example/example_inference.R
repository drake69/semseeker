result_folder = "/dati/experiments_data/DIOXIN/20220507/semseeker/"
logFolder <- result_folder

# independent variable: deve essere nalla sample sheet passata a semseeker quando lo abbiamo eseguito la prima volta
# tipo di regressioni: gaussian, poisson, binomial,
# tipi di test: wilcoxon, stats::t.test,
# tipi di correlazioni: pearson, kendall, spearman

# MUTATIONS_* ~ tcdd_mother + exam_age

# transformation to dependet variable: mutations and lesions: scale, log, log2, log10, exp, johnson, none
# depth analysis
# 1: sample level
# 2: type level (gene, DMR, cpgislan) (includes 1)
# 3: genomica area: gene, body, gene tss1550, gene whole, gene tss200,  (includes 1 and 2)
# filter_p_value report after adjusting saves only significative nominal p-value

inference_details <- expand.grid("independent_variable"= c("tcdd_mother"),
                                "covariates"=c("exam_age+tcdd_father","breast_feeding"),
                                "family_test"=c("gaussian","wilcoxon"),
                                "transformation"="log",
                                "depth_analysis"=1,
                                "filter_p_value" = FALSE)


association_analysis(inference_details = inference_details)
