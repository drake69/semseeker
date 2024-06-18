name_cleaning <- function(name)
{
  name <- toupper(name)
  name <- toupper(gsub("[.]", "_", name))
  name <- toupper(gsub(" ", "_", name))
  name <- toupper(gsub("-", "_", name))
  name <- toupper(gsub("'", "_", name))
  name <- toupper(gsub("__", "_", name))
  return(name)
}
