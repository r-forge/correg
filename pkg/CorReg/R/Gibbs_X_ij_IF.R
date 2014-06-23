#' Gibbs sampler for redundant covariates
#' @param X est la ligne concernée (X[i,]) 
#' @param components est le vecteur des classes pour i (vecteur creux de taille p avec pr elements non nuls)
#' @param mixmod is mixmod$details
Gibbs_X_ij_IF<-function(Z=Z,X=X,p=p,mui=mui,sigmai=sigmai,alpha=alpha,mixmod=mixmod,j=j,components=components){
   Sigma_j_reste=rep(0,times=p-1)
   Sigma_reste_reste=matrix(0,nrow=p-1,ncol=p-1)
   for(i in (1:p)[-j]){
      if(i!=j){
         #variable à droite de la manquante (impossible car elle est à droite) donc Z[i,j]=0
         if(Z[j,i]!=0){#variable à gauche de la manquante 
            Sigma_j_reste[i]=alpha[-1][j,i]*sigmai[i]
         }#sinon indépendance car pas à droite non plus donc rien à faire car initialisation à 0
      
        
      }
   }
   
   mu=mui[j]+Sigma_j_reste%*%solve(Sigma_reste_reste)%*%(X[-j]-mui[-j])
   sigma=0
   return(rnorm(1,mean=mu,sd=sigma))
}