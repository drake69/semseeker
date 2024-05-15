sig.formula_vars <- function(sig.formula){

  dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
  dependent_variable <- dep_var[[2]]
  independent_variable <- dep_var[[3]]
  # check if contains +
  if(grepl("\\+",independent_variable))
  {
    independent_variable <- strsplit(independent_variable,"\\+")
    independent_variable <- unlist(independent_variable)
    # take the first
    independent_variable <- independent_variable[1]
  }
  return(list(dependent_variable=dependent_variable,independent_variable=independent_variable))

}
