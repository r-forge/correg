SEM<-function(nbit_gibbs=1,nbit_SEM=50,warm=10,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
              Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir,compout=TRUE,Xout=FALSE,alphaout=TRUE){
   last=FALSE
   result=list()
   
   #initialisation
   if(is.null(alpha)){
      alpha=hatB(Z = Z,X=X)
   }
   if(is.null(Zc)){
      Zc=colSums(Z)
   }
   if(is.null(Ir)){
      Ir=which(Zc!=0)
   }
   if(is.null(comp_vect)){
      comp_vect=matrix(rep(1,times=p*rowmiss),ncol=p)
   }
   if(is.null(sigma_IR)){
      sigma_IR=rep(0,times=p)#résidus des regressions (fixés avec alpha donc stables pour gibbs)
      for (j in Ir){
         sigma_IR[j]=sd(X[,j]-X%*%alpha[-1,j])
      }
   }
   
   
   for (i in 1:(nbit_SEM+warm)){
      print(i)
      if(i==(nbit_SEM+warm)){last=TRUE}
      #SE step
      resgibbs2=Gibbs(last=last,nbit=nbit_gibbs,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                     Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir)
      comp_vect=resgibbs2$comp_vect
      X=resgibbs2$X
#       print(comp_vect)
      #step M
      resM=Mstep(Z=Z,X=X,sigma_IR=sigma_IR,Ir=Ir)
      alpha=resM$alpha
      sigma_IR=resM$sigma_IR
      if(i>warm){
         if(i==warm+1){
            if(alphaout){result$alpha=alpha/nbit_SEM}
            if(Xout){result$X=X/nbit_SEM}
         }else{
            if(alphaout){result$alpha=result$alpha+alpha/nbit_SEM}#optimiser au format creux
            if(Xout){result$X=result$X+X/nbit_SEM}#optimiser en ne modifiant que les manquants
         }   
      }
   }
   if(compout){result$comp=comp_vect}
   return(result)
}