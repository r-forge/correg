#' read the structure and explain it
#' 
#' @param Z is the binary structure
#' @param B is the complete structure (if known)
#' @param R2 boolean. Indicates wether you want to order Z by R2 (adjusted)
#' @param varnames the names of the variables (same order)
#' @param output indicates the content of the output
#' @param X is a dataframe containing the dataset
#' @export
readZ<-function(Z=Z,B=NULL,R2=T,varnames=NULL,output=c("index","names","all"),X=NULL,decr=T){
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
  if(is.null(B) & (R2| output=="all")){#if B needed but unknown
     B=hatB(Z=Z,X=X)
  }
  if(R2){
    R2vect=R2Z(Z=Z,X=X,adj=T,crit="R2")
    R2vect=R2vect[R2vect!=0]
    quiI2=which(colSums(Z)!=0)
    neworder=order(R2vect,decreasing=decr)
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
        ssreg=data.frame(cbind(c(R2vect[neworder[i]],NA,beta[beta!=0]),c("R2adj",varnames[quiI2[neworder[i]]],"intercept",varnames[which(Z[,quiI2[neworder[i]]]!=0)])))
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