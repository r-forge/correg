#' Probability of Z without knowing the dataset. It also gives the exact number of binary nilpotent matrices of size p.
#' @param p the number of covariates
#' @param Z the structure
#' @export
#' 
#'
ProbaZ<-function(p=NULL,Z=NULL,proba=FALSE){
   if(is.null(p)){
      if(!is.null(Z)){
         p=ncol(Z)
         
      }else{
         print("missing parameters")
      }
   }
   nb=1#modele vide
   if(p>1){
      #calcul du nombre de modèles
      for (i in 1:(p-1)){
         nb=nb+choose(p,i)*(2^(p-i)-1)^i
      }
   }
   if(proba){
      return(1/nb)
   }else{
      return(nb)
   }
}