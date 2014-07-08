#' Estimates B matrix
#'@param Z binary adjacency matrix of the structure (size p)
#' @param X the dataset
#' @param methode parameter for OLS (matrix inversion) methode_BIC  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
#' @export 
hatB<-function(Z=Z,X=X,methode=1){
  p=ncol(Z)
  B=matrix(0,ncol=p,nrow=p+1)
  quiI2=which(colSums(Z)!=0)
  for(i in quiI2){
    qui=which(Z[,i]!=0)
    Xloc=as.matrix(X[,qui])
    Yloc=as.matrix(X[,i])
    quimank=which(is.na(Xloc),arr.ind=T)[,1]
    quimank=c(quimank,which(is.na(Yloc),arr.ind=T))
    quimank=unique(quimank)
    if(length(quimank)>0){
      Xloc=matrix(Xloc[-quimank,],ncol=length(qui)) #si des valeurs sont manquantes,on supprimes les lignes localement
      Yloc=Yloc[-quimank]
    }
    beta=OLS(X=matrix(as.double(Xloc),ncol=ncol(Xloc)),Y=as.double(Yloc),intercept=T,methode=methode)$beta
    B[c(1,qui+1),i]=beta
    if(any(is.infinite(beta)) | any(is.nan(beta))){
      print(paste("hatB singularity col=",i,"set to 0"))
      B[c(1,qui+1),i][which(is.infinite(beta))]=0
      B[c(1,qui+1),i][which(is.nan(beta))]=0
    }
  }
  return(B)
}




