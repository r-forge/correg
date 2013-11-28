#' simple MSE function
#' @export
#' @param intercept indicates wether A contains an intercept or not
MSE_loc<-function(Y=Y,X=X,A=A,intercept=T){
  if(intercept){
    X=as.matrix(cbind(1,X))
  }
  return(mean((Y-X%*%A)^2))
}