# test_that("createExcelPivot", {
#
#   init_env()
#   populations <- c("Reference","Control","Case")
#   figures <- c("HYPO", "HYPER", "BOTH")
#   anomalies <- c("MUTATIONS","LESIONS")
#
#   subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
#   probesPrefix = "PROBES_Gene_"
#   mainGroupLabel =  "GENE"
#   subGroupLabel="GROUP"
#   createExcelPivot ( populations =  populations, figures =  figures,anomalies =  anomalies, subGroups =  subGroups, probesPrefix =   probesPrefix, mainGroupLabel =  mainGroupLabel, subGroupLabel =  subGroupLabel)
#
#
#   probesPrefix <- "PROBES_Island_"
#   subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
#   mainGroupLabel <- "ISLAND"
#   subGroupLabel <- "RELATION_TO_CPGISLAND"
#   createExcelPivot ( populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
#
#   subGroups <- c("DMR")
#   probesPrefix = "PROBES_DMR_"
#   mainGroupLabel =  "DMR"
#   subGroupLabel="GROUP"
#   createExcelPivot (populations, figures, anomalies, subGroups, probesPrefix, mainGroupLabel, subGroupLabel)
#
#   doParallel::stopImplicitCluster()
#   parallel::stopCluster(computationCluster)
#
# })
