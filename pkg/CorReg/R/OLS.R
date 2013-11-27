#' Ordinary Least Square efficiently computed with SEM for missing values
#' 
#' @param X the covariates (double)
#' @param Y the response
#' @param intercept (boolean) is an intercept intended ?
#' @param sigma (boolean) is it necessary to compute the standard deviation of errors ?
#' @param Bic (boolean) is the BIC criterion computation intended ?
#' @param methode how to compute the least squares ?
#' @param miss to indicate wether there are missing values in X (miss=1) or not (miss=0) next version would allow missing values in Y
#' @param nbit number of iteration for SEM
#' @param nbclustmax max number of cluster for mixmod (ignored if miss=0)
#' @param result of calcul_BIC_mixmod2.0(X=X,nbclustmax=nbclustmax,details=T)
#' @param M binary matrix (size of X) with 1 where X is missing
#' @param mixmod Gaussian Mixutre hypothesis if needed. Or result of calcul_BIC_mixmod(X=X,nbclustmax=nbclustmax,bla=F,details=T)
OLS<- function(X=X,Y=Y,M=NULL,intercept=F,sigma=F,Bic=F,methode=1,miss=0,mixmod=NULL,nbit=100,nbclustmax=10){
  if(miss & nbit>0){
    if(is.null(mixmod)){
      mixmod=calcul_BIC_mixmod(X=X,nbclustmax=nbclustmax,bla=F,details=T)
      mixmod=mixmod_adapter(mixmod=mixmod)
    }
    if(is.null(M)){
      M=0*X
      M[is.na(M)]=1
      X[is.na(X)]=0
    }
    quimank=which(M==1,arr.ind=T)#sparse version of M
    
    resmat=matrix(ncol=ncol(X)+intercept+sigma,nrow=nbit)
    for (i in 1:nbit){
      for (j in 1:nrow(quimank)){#optimisable en calculant probas une fois seulement par variable
        probas=mixmod$details[mixmod$details[,4]==quimank[j,2],1]#proportions de chaque classe
#         print(probas);print("ok");print(j);print(quimank[j,2])
        quelclasse=sample(length(probas),size=1,prob=probas)#on tire une classe
        #on impute selon la classe choisie
#         print(quelclasse)
#         print((mixmod$details[,4]==quimank[j,2]))
#         print((mixmod$details[,4]==quimank[j,2])[quelclasse])
        meanloc=mixmod$details[(mixmod$details[,4]==quimank[j,2]),2][quelclasse]
        sdloc=sqrt(mixmod$details[(mixmod$details[,4]==quimank[j,2]),3][quelclasse])
        val=rnorm(1,mean=meanloc,sd=sdloc)
#         print(val)
        X[quimank[j,1],quimank[j,2]]=val
      }
      res=.Call( "OLS",X, as.double(Y),intercept,sigma,F,methode, PACKAGE = "CorReg")
      resmat=rbind(resmat,c(res$sigma,res$beta))
      resmat[i,]=c(res$sigma,res$beta)
    }
    result=list()
    if(sigma){
      result$sigma=mean(resmat[,1])
      result$beta=colMeans(resmat[,-1])
    }else{
      result$beta=colMeans(resmat)
    }
    if(Bic){
      result$BIC=NULL#BicTheta(X,)#version modifiée qui tient compte des mélanges
    }
    return(result)
    
  }else{
    if(any(is.na(X))){
      ou=unique(which(is.na(X),arr.ind=T)[,1])
      X=X[-ou,]
      Y=Y[-ou]
    }
    if(any(is.na(Y))){
      ou=which(is.na(Y))
      Y=Y[-ou]
      X=X[-ou,]
    }
    if(length(Y)>0){    
      res=.Call( "OLS",X, as.double(Y),intercept,sigma,Bic,methode, PACKAGE = "CorReg")
    }else{
      print("not any complete line")
      res=NA
    }
    return (res)
  }
}