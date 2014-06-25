Gibbs<-function(last=FALSE,nbit=1,mixmod=mixmod,X=X,comp_vect=comp_vect,missrow=missrow,quimiss=quimiss,
                Z=Z,Zc=Zc,alpha=alpha,sigma_IR=sigma_IR,nbclust_vect=nbclust_vect,Ir=Ir){
   for(iter in 1:nbit){
      missrow=1
      resmui=muiZ(mixmod = mixmod,components=comp_vect[missrow,],Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR)#on commence ligne 1
      mui=resmui$mui
      sigmai=resmui$sigmai
      Sigma=resmui$Sigma
      for(i in 1:nbmiss){
#          print(paste("i",i))
         miss=quimiss[i,]
         if(miss[1]!=missrow){
            missrow=miss[1]
            resmui=muiZ(mixmod = mixmod,components=comp_vect[missrow,],Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR)#mise à jour
            mui=resmui$mui
            sigmai=resmui$sigmai
            Sigma=resmui$Sigma
         }
#          print(paste("i2",i))
#          print(Ir)
#          print(miss[2])
#          print(miss[2]%in%Ir)
         if(miss[2]%in%Ir){#redundant covariate missing
#             print("left")
            
            X[miss[1],miss[2]]=Gibbs_X_ij_IR(X=X[miss[1],],sigma=sigma_IR[miss[2]],alpha=alpha[,miss[2]])
         }else{#missing right
#             print("right")
            X[miss[1],miss[2]]=Gibbs_X_ij_IF(Z=Z,X=X[miss[1],],mui=mui,sigmai=sigmai,Sigma=Sigma,alpha=alpha,mixmod=mixmod,j=miss[2],components=comp_vect[miss[1],])
         }
#          print(paste("i3",i))
      }
      if(!last | (iter!=nbit)){
         #Imputation des classes
         for(j in 1:p){
#             print(paste("j",j))
            if(nbclust_vect[j]>1){
               for (i in 1:rowmiss){
#                   print("a")
                  vect=rmultinom(1, 1, tik(x=X[i,j],nbclust=nbclust_vect[j],mixmod=mixmod[[j]]))
#                   print("b")
#                   print(paste("vect",vect))
                  comp_vect[i,j]=match(1,vect)
               }
            }
         }
      } 
   }
#    print("c")
#    print(X)
#    print(comp_vect)
   return(list(X=X,comp_vect=comp_vect))
}