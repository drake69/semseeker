# Core wrapping function
wrap_it <- function(x, len)
{
  sapply(x, function(y) paste0(strwrap(y, len),
                              collapse = "\n"),
         USE.NAMES = FALSE)
}


