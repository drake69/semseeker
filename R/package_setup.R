# .onLoad <- function(libname, pkgname) {
#
#   op <- options()
#   op.devtools <- list(
#     devtools.path = "~/R-dev",
#     devtools.install.args = "",
#     devtools.name = "Your name goes here",
#     devtools.desc.author = '"First Last <first.last@example.com> [aut, cre]"',
#     devtools.desc.license = "What license is it under?",
#     devtools.desc.suggests = NULL,
#     devtools.desc = list()
#   )
#   toset <- !(names(op.devtools) %in% names(op))
#   if (any(toset))
#     options(op.devtools[toset])
#
#   devtools::use_package("doParallel")
#   devtools::use_package("checkr")
#   devtools::use_package("dplyr")
#   devtools::use_package("foreach")
#   devtools::use_package("gtools")
#   devtools::use_package("openxlsx")
#   devtools::use_package("parallel")
#   devtools::use_package("plyr")
#   devtools::use_package("reshape2")
#   devtools::use_package("stringi")
#   devtools::use_package("zoo")
#
#   person("Luigi", "Corsaro", email = "l.corsaro@a-company.it", role = c("aut", "cre"))
#   person("Davide", "Gentilini", email = "davide.gentilini@unipv.it", role = c("aut"))
#   person("Luciano", "Calzari", email = "luciano.calza@gmail.com", role = c("ctb"))
#
#   invisible()
# }
#
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Welcome to my package")
#
#   if (!requireNamespace("checkr", quietly = TRUE)) {
#     stop("checkr is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("doParallel", quietly = TRUE)) {
#     stop("doParallel is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("dplyr", quietly = TRUE)) {
#     stop("dplyr is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("foreach", quietly = TRUE)) {
#     stop("foreach is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("gtools", quietly = TRUE)) {
#     stop("gtools is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("openxlsx", quietly = TRUE)) {
#     stop("openxlsx is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("plyr", quietly = TRUE)) {
#     stop("plyr is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("parallel", quietly = TRUE)) {
#     stop("parallel is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("reshape2", quietly = TRUE)) {
#     stop("reshape2 is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("stringi", quietly = TRUE)) {
#     stop("stringi is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#   if (!requireNamespace("zoo", quietly = TRUE)) {
#     stop("zoo is needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#
#
# }
