#' a summary-like function
#' @param A coefficient vector
#' @param labels name of the covariates
#' @param X the dataset (named) if labels is null
#' @param intercept boolean defining wether A contains an intercept or not
#'@export
readY<-function(A=A,labels=NULL,X=NULL,intercept=TRUE){
   if(is.null(labels)){
      labels=names(X)
   }
   if(intercept){labels=c("intercept",labels)}
   interp=cbind(A[A!=0],labels[A!=0])
   return(interp)
}