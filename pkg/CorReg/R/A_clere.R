# '
A_clere<-function(y=Y,x=X,g=NULL,analysis="aic"){
   x=1*as.matrix(x)
   if(is.null(g)){
      g=min(5,ncol(x))#ncol(x)
   }
   model <- fit.clere(y = as.numeric(y), x = x, g = g, plotit = FALSE, analysis=analysis)
   A=c(model@intercept,rowMeans(model@Bw))
   return(A)
}