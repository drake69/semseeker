.pkgglobalenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {

  assign("ssEnv", list(), envir=.pkgglobalenv)
}
