# .pkgglobalenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {
  .pkgglobalenv <- new.env(parent=emptyenv())
  assign("ssEnv", NULL, envir=.pkgglobalenv)
  # assign("myvar", 42, envir=.pkgglobalenv)
}
