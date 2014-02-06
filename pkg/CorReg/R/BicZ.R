#' Compute the BIC of a given structure
#' @export
#' @param X the dataset
#' @param Z the structure
#' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
#' @param BicOld BIC (vector) associated to Zold
#' @param Zold another structure with some common parts with Z (allows to compute only the differences)
#' @param methode parameter for OLS 
#' @param star boolean defining wether classical BIC or BIC* is computed
BicZ<-function(X=X,Z=Z,Bic_null_vect=NULL,BicOld=NULL,methode=1,Zold=NULL,star=TRUE){
   if(is.null(Bic_null_vect)){
      Bic_null_vect=density_estimation(X=X)$BIC_vect 
   }
 if(is.null(Zold)| is.null(BicOld)){
   Zold=0*Z
   BicOld=Bic_null_vect
 }
    res=.Call( "BicZ",X,Z,Bic_null_vect,BicOld,methode,Zold, PACKAGE = "CorReg")
   if(star){
      res$BIC=sum(res$BIC)-ProbaZ(Z,star=TRUE)
   }
    return(res$BIC)

}