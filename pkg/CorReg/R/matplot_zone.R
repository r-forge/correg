#' @export
matplot_zone<-function(x=x,y=y,col=1:6,alpha=0.2,what=which.min,ylim=NULL){
   matplot(x,y,ylim=ylim,type="l",xlab="R2",ylab="MSE",main=paste('n=',n, ", sigmaY=",sigmaY),lwd=lwd)
   victory_int(x=x,y=y,col=col,what=what)
   matplot(x,y,ylim=ylim,add=TRUE,col=col,type="l",xlab="R2",ylab="MSE",main=paste('n=',n, ", sigmaY=",sigmaY),lwd=lwd)
}