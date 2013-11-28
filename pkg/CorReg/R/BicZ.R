#' Compute the BIC of a given structure
#' @export
#' @param X the dataset
#' @param Z the structure
#' @param Bic_vide_vect the BIC of the null hypothesis (used for independent variables)
#' @param BicOld BIC (vector) associated to Zold
#' @param Zold another structure with some common parts with Z (allows to compute only the differences)
#' @param methode parameter for OLS 
BicZ<-function(X=X,Z=Z,Bic_vide_vect=NULL,BicOld=NULL,methode=1,Zold=NULL){
 if(is.null(Bic_vide_vect)){
   val=calcul_BIC_mixmod(X=X)
   Bic_vide_vect=val$BIC_vect
 }
 if(is.null(Zold)| is.null(BicOld)){
   Zold=0*Z
   BicOld=Bic_vide_vect
 }
    res=.Call( "BicZ",X,Z,Bic_vide_vect,BicOld,methode,Zold, PACKAGE = "CorReg")
    return(res$BIC)

}