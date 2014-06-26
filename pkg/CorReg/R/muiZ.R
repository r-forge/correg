
# ' @param alpha the matrix of coefficients of the sub-regressions
muiZ<-function(mixmod=mixmod,components=components,Z=Z,Zc=Zc,alpha=alpha,Ir=Ir,sigma_IR=sigma_IR){ 
   mui=components#en C on les prendra en entrée pour gagner une initialisation à chaque fois
   sigmai=components
   for(i in (1:p)[-Ir]){#vaiable à droite donc mixmod
         mui[i]=mixmod[[i]][components[i],2]
         sigmai[i]=mixmod[[i]][components[i],3]#attention c'est une variance
   }
   for(j in Ir){
      mui[j]=alpha[1,j]+t(as.matrix(mui))%*%(alpha[-1,j])#intercept plus combi linéaire
      sigmai[j]=t(as.matrix(sigmai))%*%((alpha[-1,j])^2)+sigma_IR[j]^2#attention c'est une variance
   }
   Sigma=Matrix(0,ncol=5,nrow=5)# à intialiser au début en creux symétrique et à passer en argument
 
   diag(Sigma)=sigmai
   for(j in Ir){
      for (i in (1:p)[-j]){
         if(Zc[i]!=0){#Les deux sont à gauche
            #veiller à ne pas stocker des 0 en testant avant d'écrire ?
            Sigma[i,j]=(alpha[-1,][,i]*alpha[-1,][,j])%*%sigmai#attention premier produit hadamard
            Sigma[j,i]=Sigma[i,j]
         }else{#i est à droite
            Sigma[i,j]=alpha[-1,][i,j]*sigmai[i]
            Sigma[j,i]=Sigma[i,j]
         }
      }
      
   }
   return(list(mui=mui,sigmai=sigmai,Sigma=Sigma))
}