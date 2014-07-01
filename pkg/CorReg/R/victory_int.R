victory_int<-function(x=x,y=y,col=c("black","red","green"),pos=10){
   quimin=apply(y,1,which.min)
   points(x,rep(pos,times=length(x)),col=col[quimin],lwd="2",pch=15)
}