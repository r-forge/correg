#' draws matplot with conditionnal background for easier comparison of curves.
#' @param x the abscisses
#' @param y matrix of the curves (columns)
#' @param col list of colors (like in matplot)
#' @param what a function to choose a winner
#' @param alpha parameter for transparency of the background
#' @param ylim vector for vertical limits
#' @param xlim vector for horizontal limits
#' @param type the type of curve (like in matplot)
#' @param xlab (like in matplot)
#' @param ylab (like in matplot)
#' @param lwd  (like in matplot)
#' @param lty  (like in matplot)
#' @param main the main title (like in matplot)
#' @export
matplot_zone<-function(x=x,y=y,col=1:6,alpha=0.2,what=which.min,ylim=NULL,xlim=NULL,type="p",xlab=NULL,ylab="NULL",main=NULL,lwd=lwd,lty = 1:5){
   matplot(x,y,ylim=ylim,type=type,xlab=xlab,ylab=ylab,main=main,lwd=lwd,xlim=xlim)
   victory_int(x=x,y=y,col=col,what=what,alpha=alpha)
   matplot(x,y,ylim=ylim,add=TRUE,col=col,type=type,xlab=xlab,ylab=ylab,main=main,lwd=lwd,xlim=xlim)
}