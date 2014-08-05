#' to upgrade a package to the lastest version on R-forge
#' @param package vector of the name of the packages to upgrade
#' @export
rforge<-function(packages=c("CorReg"),repos="http://R-Forge.R-project.org"){
   install.packages(packages, repos=repos)
}