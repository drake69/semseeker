# getKeys <- function ()
# {
#
#   keys_populations <- c("Reference","Control","Case")
#   keys_figures <- c("HYPO", "HYPER", "BOTH")
#   keys_anomalies <- c("MUTATIONS","LESIONS")
#
#   keys <- expand.grid("FIGURE"=keys_figures,"ANOMALY"=keys_anomalies, "POPULATION"=keys_populations)
#
# }

#
# getProbesKey <- function()
# {
#   probes.subGroups <- c("Body","TSS1500","5UTR","TSS200","1stExon","3UTR","ExonBnd","Whole")
#
#   probes.Prefix = "PROBES_Gene_"
#   probes.MainGroupLabel =  "GENE"
#   probes.SubGroupLabel="GROUP"
#   probes <- expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups)
#
#   # probes.27k
#   # probes.450k
#   # probes.850k
#
#   probes.subGroups <- c("N_Shore","S_Shore","N_Shelf","S_Shelf","Island", "Whole")
#
#   probes.Prefix <- "PROBES_Island_"
#   probes.MainGroupLabel <- "ISLAND"
#   probes.SubGroupLabel <- "RELATION_TO_CPGISLAND"
#   probes <- rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)
#
#   probes.subGroups <- c("DMR")
#
#   probes.Prefix = "PROBES_DMR_"
#   probes.MainGroupLabel =  "DMR"
#   probes.SubGroupLabel="GROUP"
#
#   probes <- rbind(expand.grid("prefix"=probes.Prefix,"maingrouplable"= probes.MainGroupLabel,"subgrouplable"= probes.SubGroupLabel,"subgroups"= probes.subGroups), probes)
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
