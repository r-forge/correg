# '
A_clere<-function(y=Y,x=X,g=NULL){
   if(is.null(g)){
      g=ncol(x)
   }
   model <- fit.clere(y = y, x = x, g = g, plotit = FALSE)
   clus <- clusters(model, threshold = NULL)
   A=c(model@intercept,model@b[clus])
   return(A)
}