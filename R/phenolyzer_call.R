phenolyzer_call <- function(term)
{

  # cpan App::cpanminus
  # biotools:phenolyzer
  # cpanm Module::Name
  # cd /usr/local/lib/
  # git clone https://github.com/WGLab/phenolyzer.git
  # sudo cpan Bio::OntologyIO
  # sudo cpan Graph::Directed

  pkgFolder <- system.file(package="semseeker")

  cleanedTerm <- string_normalize(term)
  cacheFolder <-dir_check_and_create(pkgFolder,c("cache","Phenolyzer", cleanedTerm))

  filePrioritisedGenes <- file_path_build(cacheFolder,cleanedTerm,"final_gene_list")
  if(!file.exists(filePrioritisedGenes))
  {
    #cd hypertension
    #perl /usr/local/lib/phenolyzer/disease_annotation.pl hypertension -prediction -phenotype -logistic -out hypertension -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15
    # pcommand <- paste0("cd ", resFolder ," && perl /usr/local/lib/phenolyzer/disease_annotation.pl '" ,term, "' -p -ph -logistic -out '" ,term, "' -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25", sep="")
    pcommand <-
      paste0(
        "cd ",
        cacheFolder ,
        " && perl /usr/local/lib/phenolyzer/disease_annotation.pl '" ,
        term,
        "' -prediction -phenotype -logistic -out '" ,
        cleanedTerm ,
        "' -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 -nproc 15",
        sep = ""
      )

    system(pcommand)
  }

  prio_genes <-  read.csv(filePrioritisedGenes, sep = "\t")
  if(nrow(prio_genes)==0)
    return (NULL)
  return (prio_genes)

}
