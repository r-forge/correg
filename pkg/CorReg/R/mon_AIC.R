mon_AIC<-function(theta=theta,Y=Y,X=X,intercept=TRUE){
   if(intercept){X=cbind(1,X)}
   AIC=-2*log_likelihood(theta=theta,Y=Y,X=X)+2*length(theta[theta!=0])
   return(AIC)
}