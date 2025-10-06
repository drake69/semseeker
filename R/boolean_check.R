boolean_check <- function(x)
{
  if (is.na(x) || x=="" || is.null(x))
    return(FALSE)

  if (!is.logical(x))
    x = (x == "TRUE" || x == "T" || x == "true" || x == "True" || x == "1" || x == 1)


  return(x)
}
