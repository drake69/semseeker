string_normalize <- function (input_string)
{
  input_string <- gsub("__","_", input_string)
  input_string <- gsub(" ", "_", input_string)
  input_string <- (gsub("-", "_", input_string))
  input_string <- (gsub(":", "_", input_string))
  input_string <- (gsub("/", "_", input_string))
  input_string <- (gsub("'", "_", input_string))
  return(toupper(input_string))

}
