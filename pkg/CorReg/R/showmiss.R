#'whow the missing values of a dataset
#'
#'@param x the dataset to analyse

showmiss<-function(X){
  M=which(is.na(X),arr.ind=T)
  if(nrow(M)>1){
    plot(M[,c(2,1)],pch=7)
    title("Missing values in the dataset")
    
    
    
  }else{
    print("No missing values")
  }
}