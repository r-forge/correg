SEM<-function(M=M,nbit_gibbs=1,n=n,nbit_SEM=50,warm=10,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
              Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,loglikout=FALSE,nbclust_vect=nbclust_vect,Ir=Ir,compout=TRUE,Xout=FALSE,alphaout=TRUE,gibbsfin=0){
   last=FALSE
   result=list()
   
   #initialisation
  
   if(is.null(Zc)){
      Zc=colSums(Z)
   }
   if(is.null(Ir)){
      Ir=which(Zc!=0)
   }
   if(is.null(comp_vect)){
      comp_vect=matrix(rep(1,times=p*n),ncol=p)
#       for(i in 1:p){
#          if(nbclust_vect[i]>1){
#             comp_vect[,i]=which(rmultinom(n =missrow, size = 1,prob=mixmod[[i]][,1])==1,arr.ind = TRUE)[,1]
#          }
#       }
      #Imputation des classes####
      for(j in 1:p){
         if(nbclust_vect[j]>1){
            for (i in 1:n){
               if(!is.na(X[i,j])){#valeur observée (ou imputée, peu importe) on estime la classe
                  vect=rmultinom(1, 1, tik(x=X[i,j],nbclust=nbclust_vect[j],mixmod=mixmod[[j]]))
                  comp_vect[i,j]=match(1,vect)
               }else{#manquant donc on impute la classe et la valeur 
                  vect=rmultinom(1, 1, mixmod[[j]][,1])
                  comp_vect[i,j]=match(1,vect)
                  X[i,j]=rnorm(n=1,mean=mixmod[[j]][comp_vect[i,j],2],sd=sqrt(mixmod[[j]][comp_vect[i,j],3]))#on tire dans la classe choisie
               }
            }
         }
      }
   }
if(is.null(alpha)){
   alpha=hatB(Z = Z,X=X)
}
   if(is.null(sigma_IR)){
      sigma_IR=rep(0,times=p)#résidus des regressions (fixés avec alpha donc stables pour gibbs)
      for (j in Ir){
         sigma_IR[j]=sd(X[,j]-X%*%alpha[-1,j])
      }
   }
   
   loglik_bool=FALSE
   for (i in 1:(nbit_SEM+warm)){
      print(i)
      if(i==(nbit_SEM+warm) & gibbsfin<=0){last=TRUE}
      if(i>warm & loglikout){
         loglik_bool=TRUE
      }#on commence à calculer les vraisemblances
     
      #SE step####
      resgibbs2=Gibbs(M=M,last=last,nbit=nbit_gibbs,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                     Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir,loglik_bool=loglik_bool)
      comp_vect=resgibbs2$comp_vect
      X=resgibbs2$X
      
      #step M####
      resM=Mstep(Z=Z,X=X,sigma_IR=sigma_IR,Ir=Ir)
      alpha=resM$alpha
      sigma_IR=resM$sigma_IR
      
      if(i>warm){
         if(loglikout){#on calcule la vraisemblance locale
            loglik=resgibbs2$loglik#scalaire
         }
         if(i==warm+1){
            if(alphaout){result$alpha=alpha/nbit_SEM}
#             if(Xout){result$X=X/nbit_SEM}
            if(loglikout){result$loglik=loglik/nbit_SEM}
         }else{
            if(alphaout){result$alpha=result$alpha+alpha/nbit_SEM}#optimiser au format creux
#              if(Xout){result$X=result$X+X/nbit_SEM}#optimiser en ne modifiant que les manquants
            if(loglikout){result$loglik=result$loglik+loglik/nbit_SEM}
         }   
      }
   }
if(Xout){result$X=X}

   if(compout){result$comp=comp_vect}

if(gibbsfin>0){
   resgibbs2=Gibbs(M=M,last=TRUE,nbit=gibbsfin,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                   Z=Z,Zc=Zc,alpha=result$alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir,loglik_bool=loglik_bool)
   result$X=resgibbs2$X
   result$loglik=resgibbs2$loglik
}
   return(result)
}