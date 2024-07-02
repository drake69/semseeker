name_cleaning <- function(name)
{
  name <- toupper(name)
  name <- toupper(gsub("[.]", "_", name))
  name <- toupper(gsub(" ", "_", name))
  name <- toupper(gsub("-", "_", name))
  name <- toupper(gsub("'", "_", name))
  name <- toupper(gsub("__", "_", name))
  name <- toupper(gsub("=", "_", name))
  name <- toupper(gsub("'", "_", name))

  # replace any special character with underscore
  name <- gsub("[^a-zA-Z0-9]", "_", name)
  name <- toupper(gsub("__", "_", name))
  # remove any underscore at the end
  name <- sub("_+$", "", name)

  return(name)
}
