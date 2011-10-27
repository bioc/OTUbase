.onAttach <- function(libname, pkgname) {
  version <- packageDescription("OTUbase", field="Version")
  packageStartupMessage(.getBar())
  packageStartupMessage("Welcome to OTUbase version ", version)
}


