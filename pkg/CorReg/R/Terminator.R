#' Destructing values to have missing ones
#' @export
#' @param wrath the ratio of missing values in the output
#' @param target the dataset (matrix or data.frame) in which missing values will be made
#' @param diag if >0 it is the thickness of the diagonal band of missing values
#' @param Z adjacency matrix to coerce a maximum of 1 missing value per sub-regression for each individual
Terminator<-function(target="Sarah Connor", wrath=0.1,diag=0,Z=NULL){
   if(target[1]=="Sarah Connor" ){
      print("I'll be back !")
   }else if (target[1]=="bender"){
      Bender()
   }else if(diag>0){
      n=nrow(target)
      p=ncol(target)
      quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
      target[quidiag]=NA
      for(j in 2:diag){
         quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(c(j:p,1:(j-1)),length.out=max(n,p)))
         quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
         target[quidiag]=NA
      }
   }else if (wrath>0){
      target=as.matrix(target)
      n=nrow(target)
      p=ncol(target)
      quidiag=cbind(rep(1:n,length.out=max(n,p)),rep(1:p,length.out=max(n,p)))
      nbmank=floor(wrath*n*p)
      loc=target
      target=NA*target
      if(!is.null(Z)){
         Zc=colSums(Z)
         quiZ_vect=which(Z!=0,arr.ind=TRUE)#liste des impliqués en ssreg
         quiZtot=unique(c(quiZ_vect))#liste des impliqués en ssreg
         for (i in 1:n){
            quiZ=quiZtot#liste des impliqués en ssreg
            quiblok=c()
            for(j in 1:length(quiZ)){#maxi pr manquants
               if(length(quiZ)>0){
                  mankloc=sample(quiZ,size=1)#on tue quelqu'un
                  quiZ=quiZ[quiZ!=mankloc]#le mort n'est plus candidat ni bloquable
                  if(Zc[mankloc]>0){#variable à gauche, on retire la regression
                     reste=which(Z[,mankloc]!=0)
                     quiblok=c(quiblok,reste)
                     quiZ=unique(c(reste,quiZ))[-c(1:length(reste))]#on enleve les bloqués
                  }else{#variable à droite, on retire la gauche, les autres à droites, et les autres régressions touchées
                     impact=quiZ_vect[quiZ_vect[1,]==mankloc,2]
                     for(k in impact){
                        if(length(quiZ)>0){
                           ssreg=c(k,which(Z[,k]!=0))
                           ssreg=ssreg[ssreg!=mankloc]
                           quiblok=unique(c(quiblok,ssreg))
                           quiZ=unique(c(quiblok,quiZ))[-c(1:length(quiblok))]#on enleve les bloqués
                        }else{
                           break
                        }
                     }
                  }                  
               }else{
                  break
               }
            }
            quidiag=rbind(quidiag,cbind(i,unique(quiblok)))
         }
      }
      target[quidiag]=loc[quidiag]
      mankmax=n*p-length(quidiag)
      
      if(nbmank<mankmax){
         candidat=which(is.na(target),arr.ind=TRUE)
         quimank=sample(1:nrow(candidat),size=nbmank)
         target=loc
         target[candidat[quimank,]]=NA
      }
   }
   return(target)
}