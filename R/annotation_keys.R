annotation_keys <- function( selection_area = c("GENE","ISLAND","DMR"))
{

  anomalies <- c("MUTATIONS", "LESIONS")
  figures <- c("HYPO", "HYPER", "BOTH")

  ssEnv$keys <- data.frame(
    "ANOMALY" = "",
    "FIGURE" = "",
    "GROUP" = "",
    "SUBGROUP" = ""
  )
  ssEnv$keys <- ssEnv$keys[-1,]

  group =  "GENE"
  subGroups <-c("Body","TSS1500","5'UTR","TSS200","1stExon","3'UTR","ExonBnd","Whole")

  ssEnv$keys_gene <-
    expand.grid(
      "ANOMALY" = anomalies,
      "FIGURE" = figures,
      "GROUP" = group,
      "SUBGROUP" = subGroups
    )

  if("GENE" %in% selection_area)
    ssEnv$keys <- rbind(ssEnv$keys, ssEnv$keys_gene)

  group <- "ISLAND"
  subGroups <- c("N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "Island", "Whole")

  ssEnv$keys_island <-
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      )

  if("ISLAND" %in% selection_area)
    ssEnv$keys <- rbind(ssEnv$keys, ssEnv$keys_island)

  group =  "DMR"
  subGroups <- c("DMR")

  ssEnv$keys_dmr <-
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      )

  if("DMR" %in% selection_area)
    ssEnv$keys <- rbind(ssEnv$keys, ssEnv$keys_island)

  return (ssEnv$keys)
}
