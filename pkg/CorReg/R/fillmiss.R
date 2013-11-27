#' Fill the missing values in the dataset
#' 
#' @param X the dataset (matrix) with missing values
#' @param Z the structure associated to X. Can be a matrix of zeros if non structure.
#' @param res_mixmod the best results found by mixmod if already computed
#' @param mixmod boolean to say if the function has to use mixmod hypothesis or just the observed mean.
#' @param nbclustmax the max number of cluster for mixmod.
#' @param X1 booleant to say if dependent variables on the right will be filled based on the structure
#' @param Bt the matrix used for X1 if X1=TRUE

fillmiss<-function(X=X,Z=NULL,mixmod=F,B=NULL,Bt=NULL,res_mixmod=NULL,nbclustmax=10,X1=T){
   quimank=which(is.na(X),arr.ind=T)
   p=ncol(X)
   if(is.null(Z)){
     Z=matrix(0,ncol=p,nrow=p)
   }
   quiou=WhoIs(Z=Z,I1=T,I2=T,I3=T)
   if(is.null(B)){
     B=hatB(Z=Z,X=X)
   } 
   if(X1 & is.null(Bt)){
     Bt=hatB(Z=t(Z),X=X)
   }
   simple=F
   if(!anyDuplicated(quimank[1,])){
      simple=T#maximum un trou par ligne
   }else{
     print("trous multiples")
   }
   for(miss in 1:nrow(quimank)){#pour chaque valeur manquante
      if(anyDuplicated(c(quimank[miss,2],quiou$I2))){#Si La valeur manquante est à gauche
         X[quimank[miss,1],quimank[miss,2]]=B[1,quimank[miss,2]]+as.matrix(X[quimank[miss,1],-quimank[miss,2]])%*%as.matrix(B[-1,quimank[miss,2]][-quimank[miss,2]])
      }else if(X1 & !anyDuplicated(c(quimank[miss,2],quiou$I3))){#si la variable est à droite et qu'on en tient compte
         beta=Bt[,quimank[miss,2]]
         X[quimank[miss,1],quimank[miss,2]]=beta[1]+as.matrix(X[quimank[miss,1],-quimank[miss,2]])%*%beta[-1][-quimank[miss,2]]
      }else if(mixmod==T){#mixmod
        if(is.null(res_mixmod)){#si mixmod n'a pas tourné, on le fait tourner
          n=nrow(X)
          nbclustmax=round(min(nbclustmax,1+n^(0.3)))
          Xloc=X[,quimank[miss,2]]
          Xloc=Xloc[!is.na(Xloc)]
          res_mixmod=mixmodCluster(data=Xloc,criterion="BIC",nbCluster=c(1:nbclustmax))["bestResult"]
        } 
          resmix=list()
          resmix$nbclust=res_mixmod[1]
          resmix$prop=c(res_mixmod[6][1])#vecteur des proportions
          resmix$mean=c(res_mixmod[6][2])#vecteur des moyennes
          resmix$var=as.numeric(res_mixmod[6][3])#vecteur des variances
          quelclust=1+sum(cumsum(resmix$prop)<runif(n=1))#on tire la classe (on peut aussi faire avec sample(,prob=prop))
          X[quimank[miss,1],quimank[miss,2]]=rnorm(n=1,mean=resmix$mean[quelclust],sd=sqrt(resmix$var[quelclust]))   
      }else{#remplissage par la moyenne
        Xloc=X[,quimank[miss,2]]
        Xloc=Xloc[!is.na(Xloc)]
        X[quimank[miss,1],quimank[miss,2]]=mean(Xloc)
      }  
   }
   return(X)
}

