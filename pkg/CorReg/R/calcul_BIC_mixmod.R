#' calcul du BIC initial via package externe Mixmod (ou mclust)
#' @export
#' @param X the dataset (matrix)
#' @param nbclsutmax max number of clusters in the gaussian mixtures
#' @param bla verbose or not
#' @param details boolean to give the details of the mixtures found
#' @param max boolean. Use an heuristic to shrink nbclustmax according to the number of individuals in the dataset
#' @param mclust boolean. Use mclust instead of Rmixmod
calcul_BIC_mixmod<-function(X=X,nbclustmax=10,bla=FALSE,details=FALSE,max=TRUE,mclust=TRUE){
  #X est la matrice sans constante
  n=nrow(X)
  if(max){
    nbclustmax=round(min(nbclustmax,1+n^(0.3)))
  }
  p=ncol(X)
  nbclust=c()
  BIC_vect=c()
  if(details){
    detailsmat=list() 
  }
  if(mclust==F){#si on veut utiliser mixmod
    for (i in 1:p){
      vect=X[!is.na(X[,i]),i]#donnees observees seulement
      res=mixmodCluster(data=vect,criterion="BIC",nbCluster=c(1:nbclustmax))["bestResult"]
      if(bla){print(res)}
      nbclust[i]=res[1]
      BIC_vect[i]=res[3]
      if(details){
        prop=res[6][1]#proportions
        meansvect=c(res[6][2])#means
        varvect=unlist(res[6][3])#variances
        detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
      }
    }
  }else{#on utilise mclust
    options(warn=-1)
    for (i in 1:p){
      vect=X[!is.na(X[,i]),i]#donnees observees seulement
      res=Mclust(vect,G=c(1:nbclustmax),modelNames="V")[c("bic","parameters")]
      nbclust[i]=res$parameters$variance$G
      BIC_vect[i]=-res$bic
      if(details){
        prop=res$parameters$pro#proportions
        meansvect=res$parameters$mean#means
        varvect=res$parameters$variance$sigmasq#variances
        detailsmat[[i]]=cbind(prop,meansvect,varvect,i)
      }
    }
    options(warn=1)
  }
  if(details){#boucle à la main pour sortie utilisable
    return(list(BIC_vect=BIC_vect,nbclust=nbclust,BIC=sum(BIC_vect),details=detailsmat))
  }else{
    return(list(BIC_vect=BIC_vect,nbclust=nbclust,BIC=sum(BIC_vect)))
  }
}