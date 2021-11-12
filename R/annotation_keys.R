annotation_keys <- function( selection_area = c("GENE","ISLAND","DMR"))
{

  anomalies <- c("MUTATIONS", "LESIONS")
  figures <- c("HYPO", "HYPER", "BOTH")

  keys <- data.frame(
    "ANOMALY" = "",
    "FIGURE" = "",
    "GROUP" = "",
    "SUBGROUP" = ""
  )
  keys <- keys[-1,]

  group =  "GENE"
  subGroups <-c("Body","TSS1500","5'UTR","TSS200","1stExon","3'UTR","ExonBnd","Whole")

  keys_gene <-
    expand.grid(
      "ANOMALY" = anomalies,
      "FIGURE" = figures,
      "GROUP" = group,
      "SUBGROUP" = subGroups
    )

  if("GENE" %in% selection_area)
    keys <- rbind(keys, keys_gene)

  group <- "ISLAND"
  subGroups <- c("N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "Island", "Whole")

  keys_island <-
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      )

  if("ISLAND" %in% selection_area)
    keys <- rbind(keys, keys_island)

  group =  "DMR"
  subGroups <- c("DMR")

  keys_dmr <-
      expand.grid(
        "ANOMALY" = anomalies,
        "FIGURE" = figures,
        "GROUP" = group,
        "SUBGROUP" = subGroups
      )

  if("DMR" %in% selection_area)
    keys <- rbind(keys, keys_island)

  return (keys)
}
