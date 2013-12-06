#' Computes the MSE on the joint distribution of the dataset
#'@export
MSEZ<-function(X=X,X_appr=NULL,B=NULL,Z=Z,scale=TRUE){
   X=scale(X)
   I1=which(colSums(Z)==0)
   X1=cbind(1,X[,I1])
   X2=X[,-I1]
   I2=(1:ncol(Z))[-I1]
   res=sum(apply(X1,2,var))
   
   if(is.null(B)){
      if(is.null(X_appr)){X_appr=X
      }else{
         X_appr=scale(X_appr)
      }
      B=hatB(Z=Z,X=X_appr)
   }
   X2=X[,-colSums(Z)!=0]
   res=res+sum(apply(X2-X1%*%B[-(I2+1),I2],2,var))
   return(res)
}