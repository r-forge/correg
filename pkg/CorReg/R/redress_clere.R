# '
redress_clere<-function(model){
   clus <- clusters(model, threshold = NULL)
   p=model@p
   A=rep(0,times=p+1)
   A[1]=model@intercept#intercept
   A[-1]=model@b[clus]
   return(A)
}