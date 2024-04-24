disease_get <- function(disease, vocabulary="HPO")
{

  # check if exisists pehnotype.hpoa file if not exit
  if (!file.exists("~/Downloads/phenotype.hpoa"))
  {
    log_event("Phenolizer: phenotype.hpoa file not found")
    return(data.frame("DISEASE"="")[-1,])
  }
  # read file and skip first 4 rows
  diseases <- read.csv2("~/Downloads/phenotype.hpoa", header = TRUE, stringsAsFactors = FALSE, sep ="\t", skip = 4)
  message(disease)
  # subset diseases based on disease name contained in the disease column
  diseases <- diseases[grep(disease, diseases$disease_name, ignore.case = TRUE),]
  if (is.null(nrow(diseases)))
    return(data.frame("DISEASE"="")[-1,])

  if (vocabulary=="HPO")
    vocabulary <- "hpo_id"

  if (vocabulary=="OMIM")
    vocabulary <- "database_id"

  diseases <- diseases[,vocabulary]
  if (length(diseases)==0)
  {
    diseases <- data.frame("DISEASE"="")[-1,]
    }
  else
    diseases <- data.frame("DISEASE"=unique(diseases))
  return(diseases)
}
