#' calcul du vecteur BIC de la matrice Z
#' 

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