#' simple MSE function
#' @export
#' @param Y the response variable (vector)
#' @param X the dataset (matrix of covariates)
#' @param A the vector of coefficients
#' @param intercept boolean (to add a column of 1 to X if A contains an intercept and X don't)
MSE_loc<-function(Y=Y,X=X,A=A,intercept=T){
  if(intercept){
    X=as.matrix(cbind(1,X))
  }
  return(mean((Y-X%*%A)^2))
}