#' Compute the log-likelihood of a gaussian mixture value due to missing values in a regression
#' @param X is the vector of the regressors
#' @param Y is the value of the response variable (scalar)
#' @param B is the vector of regression parameters, including the intercept in first position  (zero if no intercept)
#' @param sigma is the standard deviation of the residual of the subregression
#' @param M is the indicatrice vector (submatrix) of the missing values
#' @param mixmod is result of calcul_BIC_mixmod2.0(X=dataset,nbclustmax=nbclustmax,details=T)
#' @param log boolean to define if you want the likelihood or the log-likelihood
GM_Loglikelihood<-function(Y=Y,X=X,B=B,sigma=sigma,M=NULL,mixmod,log=T,intercept=T){
  if(is.null(M)){
    M=0*X
    M[is.na(M)]=1
    X[is.na(X)]=0
  }
  return(.Call( "GM_likelihood",Y,X,B,sigma,M,mixmodcpp$nbclust,mixmodcpp$details,log,intercept, package = "CorReg"))
#   res=0
#   #on calcule d'abord la partie fixe qui sera commune (ajoutée) à toutes les classes
#   quibon=which(!is.na(X))
#   meanfix=B[1]+sum(X[quibon]*B[-1][quibon])#attention au passage en C++ quand quibon est vide (tester la taille)
#   #si tout est bon, on est sur cette loi là :
#   proptot=1
#   meantot=meanfix
#   vartot=sigma^2
#   #calcul des paramètres####
#   quimank=which(is.na(X))  
#   for(i in quimank){#pour chaque variable manquante
#     proptot=kronecker(proptot,mixmod$details[[i]][,1])
#     meantot=kroneckersum_vect(meantot,mixmod$details[[i]][,2]*B[-1][i])
#     vartot=kroneckersum_vect(vartot,mixmod$details[[i]][,3]*(B[-1][i])^2)    
#   }
#   #calcul de la vraisemblance
#   for(i in 1:length(proptot)){
#     res=res+proptot[i]*dnorm(Y,mean=meantot[i],sd=sqrt(vartot[i]))
#   }
#   if(log){
#     res=log(res)
#   }
  return(res)
}


