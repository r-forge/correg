#' read the structure and explain it
#' 
#'@param Z binary adjacency matrix of the structure (size p)
#' @param B is the complete structure (wheighted)
#' @param crit define the criterion to use c("none","R2","F","sigmaX")
#' @param varnames the names of the variables (same order)
#' @param output indicates the content of the output output=c("index","names","all")
#' @param X is a dataframe or matrix containing the dataset
#' @param order Define the order used (0: none, -1: decreasing, 1: growing)
#' @export
readZ<-function(Z=Z,B=NULL,crit=c("none","R2","F","sigmaX"),varnames=NULL,output=c("index","names","all"),X=NULL,order=1){
  p=ncol(Z)
  output=output[1]
  res=list()
  if(output!="index" & is.null(varnames)){#if names needed but unknown
    if(length(names(X))!=p){
      varnames=1:p
    }else{
      varnames=names(X)
    }
  }
  if(is.null(B) & (crit!="none"| output=="all")){#if B needed but unknown
     B=hatB(Z=Z,X=X)
  }
  if(crit!="none"){
    R2vect=R2Z(Z=Z,X=X,adj=TRUE,crit=crit)
    R2vect=R2vect[R2vect!=0]
    quiI2=which(colSums(Z)!=0)
    if(order==1){#increasing order
       neworder=order(R2vect,decreasing=FALSE)
       
    }else if(order==-1){
       neworder=order(R2vect,decreasing=TRUE)
       
    }else{
       neworder=1:length(quiI2)     
    }
    if(output=="index"){
      for(i in 1:length(quiI2)){
        res[[i]]=c(quiI2[neworder[i]],which(Z[,quiI2[neworder[i]]]!=0))
      }
    }else if(output=="names"){
      for(i in 1:length(quiI2)){
        res[[i]]=varnames[c(quiI2[neworder[i]],which(Z[,quiI2[neworder[i]]]!=0))]
      }
    }else{#all
      for(i in 1:length(quiI2)){
        beta=B[,quiI2[neworder[i]]]
        ssreg=data.frame(cbind(c(R2vect[neworder[i]],NA,beta[beta!=0]),c(crit,varnames[quiI2[neworder[i]]],"intercept",varnames[which(Z[,quiI2[neworder[i]]]!=0)])))
        names(ssreg)=c("coefs","var")
        res[[i]]=ssreg
      }
    }
  }else{
    quiI2=which(colSums(Z)!=0)
    if(output=="index"){
      for(i in 1:length(quiI2)){
        res[[i]]=c(quiI2[i],which(Z[,quiI2[i]]!=0))
      }
    }else if(output=="names"){
      for(i in 1:length(quiI2)){
        res[[i]]=varnames[c(quiI2[i],which(Z[,quiI2[i]]!=0))]
      }
    }else{#all
      for(i in 1:length(quiI2)){
        beta=B[,quiI2[i]]
        ssreg=data.frame(cbind(c(NA,beta[beta!=0]),c(varnames[quiI2[i]],"intercept",varnames[which(Z[,quiI2[i]]!=0)])))
        names(ssreg)=c("coefs","var")
        res[[i]]=ssreg
      }
    }
  }
  return(res)
}