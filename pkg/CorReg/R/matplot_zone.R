#' @export
matplot_zone<-function(x=x,y=y,col=1:6,alpha=0.2,what=which.min,ylim=NULL,type="p",xlab=NULL,ylab="NULL",main=NULL,lwd=lwd){
   matplot(x,y,ylim=ylim,type=type,xlab=xlab,ylab=ylab,main=main,lwd=lwd)
   victory_int(x=x,y=y,col=col,what=what)
   matplot(x,y,ylim=ylim,add=TRUE,col=col,type=type,xlab=xlab,ylab=ylab,main=main,lwd=lwd)
}