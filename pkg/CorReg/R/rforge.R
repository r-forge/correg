#' to upgrade a package to the lastest version on R-forge
#' @param package the name of the package to upgrade
#' @export
rforge<-function(package="CorReg"){
   install.packages(package, repos="http://R-Forge.R-project.org")
}