#' Curve of the BIC for each possible p2 with a fixed Z
#' @export
BicZcurve<-function(X=X,Z=Z,Bic_null_vect=Bic_null_vect,plot=T,star=F,trunc=NULL){
   p2=sum(colSums(Z)!=0)
   curve=sum(BicZ(X=X,Z=0*Z,Bic_null_vect=Bic_null_vect,star=star))
   
   if(p2>0){
      I2=which(colSums(Z)!=0)
      sigmavect=R2Z(Z=Z,X=X,crit="R2",adj=T)
      ordre=order(sigmavect[I2],decreasing=T)
      for (i in 1:p2){
         Zloc=Z;Zloc[,-I2[ordre[1:i]]]=0
         curve=c(curve,sum(BicZ(X=X,Z=Zloc,Bic_null_vect=Bic_null_vect,star=star)))
      }
      if(plot){
         plot(curve[-1])
         abline(h=curve[1])
      }
   }   
   quimin=which.min(curve)
   if(quimin>1){
      quimin=quimin-1
      Zopt=Z;Zopt[,-I2[ordre[1:quimin]]]=0
      if(plot){
         abline(v=quimin)
      }
   }else{
      Zopt=0*Z
   }
   if(!is.null(trunc)){
      trunc=min(trunc,p2)
      trunc=max(0,trunc)
      Zopt=Z;Zopt[,-I2[ordre[1:trunc]]]=0
      if(plot){
         abline(v=trunc,col="red")
      }
   }
   return(list(curve=curve,Zopt=Zopt))
}