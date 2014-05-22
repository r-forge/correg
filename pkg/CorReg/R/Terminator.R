#' Destructing values to have missing ones
#' @export
#' @param wrath the ratio of missing values in the output
#' @param target the dataset (matrix or data.frame) in which missing values will be made
#' 
Terminator<-function(target="Sarah Connor", wrath=0.1){
   if(target[1]=="Sarah Connor" ){
      print("I'll be back !")
   }else if (wrath>0){
      wrath=
      target=as.matrix(target)
      n=nrow(target)
      p=ncol(target)
      mankmax=n*p-max(n,p)
      nbmank=floor(wrath*n*p)
      quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
      loc=target
      target=NA*target
      target[quidiag]=loc[quidiag]
      if(nbmank<mankmax){
         candidat=which(is.na(target),arr.ind=TRUE)
         quimank=sample(1:nrow(candidat),size=nbmank)
         target=loc
         target[candidat[quimank,]]=NA
      }
   }
   return(target)
}