#' Compute the BIC of a given structure
#' @export
#' @param X the dataset
#'@param Z binary adjacency matrix of the structure (size p)
#' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
#' @param Bic_old BIC (vector) associated to Zold
#' @param Zold another structure with some common parts with Z (allows to compute only the differences)
#' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
#' @param star boolean defining wether classical BIC or BIC* is computed
BicZ<-function(X=X,Z=Z,Bic_null_vect=NULL,Bic_old=NULL,methode=1,Zold=NULL,star=FALSE){
   if(is.null(Bic_null_vect)){
      Bic_null_vect=density_estimation(X=X)$BIC_vect 
   }
 if(is.null(Zold)| is.null(Bic_old)){
   Zold=0*Z
   Bic_old=Bic_null_vect
 }
    res=.Call( "BicZ",as.matrix(X),Z,Bic_null_vect,Bic_old,methode,Zold, PACKAGE = "CorReg")
   if(star){
      res$BIC=sum(res$BIC)-ProbaZ(Z,star=TRUE)
   }
    return(res$BIC)

}