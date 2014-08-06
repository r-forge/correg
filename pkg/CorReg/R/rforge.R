#' Upgrades a package to the lastest version on R-forge
#' @param package name of the packages to upgrade
#' @export
rforge<-function(package="CorReg"){
   detach(unload=TRUE)
   install.packages(package, repos="http://R-Forge.R-project.org")
}

