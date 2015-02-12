#' show the missing values of a dataset
#'
#'@param X the dataset to analyse
#'@param what indicates what to plot
#'@param pch for missing, symbol to plot (can set pch="." for large datasets)
#'@export

showdata<-function(X=X,what=c("miss","correl"),pch=7){
   what=what[1]
   if(what=="miss"){
      M=which(is.na(X),arr.ind=T)
      if(nrow(M)>1){
         plot(M[,c(2,1)],pch=pch)
         title("Missing values in the dataset")  
      }else{
         print("No missing values")
      }
   }else if (what=="correl"){
      corrplot(corr=X,addrect=NULL,is.corr=T,method="color",tl.pos="n",diag=F,outline=F)
   }else{
      correl=cor(X[,!is.na(colSums(X)) & apply(X,2,sd)!=0])
      corrplot(corr=correl,addrect=NULL,is.corr=T,method="color",tl.pos="n",diag=F,outline=F)
   }
}
