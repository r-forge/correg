#' Gibbs sampler for redundant covariates
#' 
Gibbs_X_ij_IR<-function(X=X,sigma=sigma,alpha=alpha){
   mu=alpha[1]+X%*%alpha[-1]
   return(rnorm(1,mean=mu,sd=sigma))
}