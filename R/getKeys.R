# getKeys <- function ()
# {
#
#   ssEnv$keys_populations <- c("Reference","Control","Case")
#   ssEnv$keys_figures <- c("HYPO", "HYPER", "BOTH")
#   ssEnv$keys_anomalies <- c("MUTATIONS","LESIONS")
#
#   ssEnv$keys <- expand.grid("FIGURE"=ssEnv$keys_figures,"ANOMALY"=ssEnv$keys_anomalies, "POPULATION"=ssEnv$keys_populations)
#
# }

#
# getProbesKey <- function()
# {
#   probes_subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
#
#   probes_Prefix = "PROBES_Gene_"
#   probes_MainGroupLabel =  "GENE"
#   probes_SubGroupLabel="GROUP"
#   probes <- expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups)
#
#   # probes.27k
#   # probes.450k
#   # probes.850k
#
#   probes_subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
#
#   probes_Prefix <- "PROBES_Island_"
#   probes_MainGroupLabel <- "ISLAND"
#   probes_SubGroupLabel <- "RELATION_TO_CPGISLAND"
#   probes <- rbind(expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups), probes)
#
#   probes_subGroups <- c("DMR")
#
#   probes_Prefix = "PROBES_DMR_"
#   probes_MainGroupLabel =  "DMR"
#   probes_SubGroupLabel="GROUP"
#
#   probes <- rbind(expand.grid("prefix"=probes_Prefix,"maingrouplable"= probes_MainGroupLabel,"subgrouplable"= probes_SubGroupLabel,"subgroups"= probes_subGroups), probes)
#
# }
#
# getProbes <- function(genomicArea)
# {
#
#   probesKey <- getProbesKey()
#   probesKey <- subset(probesKey, MainGroupLabel == genomicArea)
#
#   get(paste0(probes$probesPrefix, genomicPart ,sep=""))
#   return()
# }
#
