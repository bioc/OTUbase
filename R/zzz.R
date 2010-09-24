.onAttach <- function(libname, pkgname) {
  version <- packageDescription("OTUbase", field="Version")
  message(.getBar())
  message("Welcome to OTUbase version ", version)
}


