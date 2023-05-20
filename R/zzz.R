.pkgglobalenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname) {

  if(Sys.info()[['sysname']]=="Windows")
    options(parallelly.fork.enable= FALSE)
  else
    options(parallelly.fork.enable= TRUE)

  # .pkgglobalenv <- new.env(parent=emptyenv())
  assign("ssEnv", NULL, envir=.pkgglobalenv)
  # assign("myvar", 42, envir=.pkgglobalenv)
}
