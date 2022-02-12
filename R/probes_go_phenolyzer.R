probes_go_association_phenolizer <- function(terms, onlySeed = TRUE, resultFolder)
{
   init_env(resultFolder)
  for (term in terms)
  {
    # term <- "HP:0000798"
    prio_genes <- phenolyzer_call(term)
    if(onlySeed)
      prio_genes <- subset(prio_genes, prio_genes$Status=="SeedGene")
    if(is.null(prio_genes))
      next
    # geneMutated <- geneMutated[order(geneMutated$Case, decreasing = TRUE),]
    probes <- get("PROBES_Gene_Whole")
    selectegGenes <- probes [ string_normalize(probes$GENE) %in% string_normalize(prio_genes$Gene), ]
    selectedProbes <- unique(selectegGenes$PROBE)
    # cg07057042 seed FALSE
    # cg03310027 seed TRUE
  }
  return(selectedProbes)
}
