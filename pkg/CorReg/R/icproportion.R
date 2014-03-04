icproportion<-function(prop=prop, n=n, alpha=0.05){
   inter=alpha/2#seuil de confiance
   lower =prop + qnorm(inter) * sqrt(prop * (1 -prop)/n)
   upper =prop + qnorm(1-inter) * sqrt(prop * (1 -prop)/n)
   return(c(lower, upper))
}