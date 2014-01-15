#' Probability of Z without knowing the dataset. It also gives the exact number of binary nilpotent matrices of size p.
#' @param p the number of covariates
#' @param Z the structure
#' @param star gives the log pra under uniform law for p2
#' @param proba gives the proba under the uniform law for Z
#' @export
#' 
#'
ProbaZ<-function(Z=NULL,p=NULL,proba=FALSE,star=TRUE){
   if(star & !is.null(Z)){
      p=ncol(Z)
      p1j=colSums(Z)
      I2=which(p1j!=0)
      p1j=p1j[I2]
      p2=length(I2)
      logproba=0
      logproba=logproba-log(p2)-p2*log(p-p2)-log(choose(p,p2))
      for (j in 1:p2){
         logproba=logproba-log(choose((p-p2),p1j[j]))
      }
      return(logproba)
   }else{
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
}