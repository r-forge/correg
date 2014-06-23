
# ' @param alpha the matrix of coefficients of the sub-regressions
muiZ<-function(mixmod=mixmod,components=components,Z=Z,Zc=Zc,alpha=alpha){
   
   mui=components#en C on les prendra en entrée pour gagner une initialisation à chaque fois
   sigmai=components
   for(i in 1:p){
      mui[i]=mixmod[[i]][components[i],2]
      sigmai[i]=mixmod[[i]][components[i],3]#attention c'est une variance
   }
   Sigma=0*Z# à intialiser au début en creux symétrique et à passer en argument
   for(i in 1:p){
      for(j in 1:i){
         if(i==j){#variance sur la diagonale
            if(Zc[i]>0){#variable à gauche
               Sigma[i,i]=sigmai[i]+alpha[-1][,j]^2%*%sigmai#les alphas nuls font le tri tout seul
            }else{#variable à droite
               Sigma[i,i]=sigmai[i]
            }
         }else if (Zc[i]!=0){#si une à gauche
            if(Zc[j]!=0){#Les deux sont à gauche
               #veiller à ne pas stocker des 0 en testant avant d'écrire ?
               Sigma[i,j]=(alpha[-1][,i]*alpha[-1][,j])%*%sigmai#attention premier produit hadamard
               Sigma[j,i]=Sigma[i,j]
            }else{#j est à droite
               Sigma[i,j]=alpha[-1][j,i]*sqrt(sigmai[j])
               Sigma[j,i]=Sigma[i,j]
            }
         }#sinon, deux distinctes droite alors covariance nulle
      }
   }
   return(list(mui=mui,sigmai=sigmai,Sigma=Sigma))
}