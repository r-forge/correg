#' Imputation of missing values knowing alpha (E step of the EM)
#' @param X the dataset with missing values
#' @param M binary matrix X-sized (1=missing) 
#' @param Ir the set of indices of the covariates on the left (computed once before EM to be efficient)
#' @param Zc profil (sum) colonne de Z
#' @param mixmod the matrix from mixmod
#' @param Zc p-sized vector to see if Z column is empty or not
#' @param alpha matrix of the coefficients (p+1)xp format Matrix
#'
Estep<-function(X=X,alpha=alpha,M=NULL,Z=NULL,mixmod=mixmod,Zc=Zc,X1=FALSE){
   #X bouge, mais pas alpha ni Z ni M ni mixmod
   require(Matrix)
   quimank=which(is.na(X),arr.ind=TRUE)# apriori on stockera les deux formats (M et quimank en permanence) et on fera pareil pour Z
   
   p=ncol(X)
   n=nrow(X)
   #on liste tout ce qu'il faut (on verra après comment extraire ça pour optimiser)
   nbmank=nrow(quimank)
   M=sparseMatrix(i=quimank[,1],j=quimank[,2],dims=c(n,p))
   Z=as(Matrix(Z,sparse=TRUE),"nsparseMatrix")
   Zc=Matrix(colSums(Z))
   X=fillmiss(X=X,X1=FALSE,mixmod=FALSE)
   
   if(is.null(mixmod)){print("ok")
      mixmod=density_estimation(X=X,detailed=TRUE,matshape=TRUE)$details
   }
   #peut-être entrer quimank et écrire uniquement dans le Cpp M au format Eigen creux (sans contact avec R)
   #on remplit
   for (miss in 1:nbmank){#pour chaque valeur manquante
      miss=quimank[miss,]
      if(Zc[miss[2]]!=0){#trou à gauche
         #imputation par moyenne du mélange gaussien résultant de l'intégration de la sous-regression
         quimanklign=which(M[miss[1],]==1 & Z[,miss[2]]!=0)
         #la constante
         loc=alpha[1,miss[2]]
         #la régression (hors manquants)
         if(length(quimanklign)>0){
            loc=loc+t(as.matrix(X[miss[1],-quimanklign]))%*%as.matrix(alpha[-1,miss[2]][-quimanklign])
         }else{
            loc=loc+t(as.matrix(X[miss[1],]))%*%as.matrix(alpha[-1,miss[2]])
         }
         #les moyennes pondérées par alpha des manquants à droite
         for(i in quimanklign){
            if(i!=miss[2]){#inutile si on a filtré sur les z non-nuls
               mixloc=matrix(mixmod[mixmod[,4]==i,c(1,2)],ncol=2)
               mixloc=sum(unlist(mixloc[,1])*unlist(mixloc[,2]))
               loc=loc+alpha[-1,miss[2]][i]*mixloc
            }
         }
         X[miss[1],miss[2]]=loc
      }else if(X1){#trou à droite
         Irj= which(Z[miss[2],]!=0)#variables expliquées par le manquant
         IrjO=which(Z[miss[2],]!=0 & M[miss[1],]==0)
         if(length(IrjO)<1){#variable dans I3 ou trop de manquants
            #imputation par la moyenne, mais déjà fait à la base
         }else{#imputation conditionnelle
            pi_ij=unlist(mixmod[mixmod[,4]==miss[2],1])
            for(l in IrjO){
               pi_ij=kronecker(pi_ij,unlist(mixmod[mixmod[,4]==l,1]))
            }
            Kij=length(pi_ij)
            loc=0
            for (k in 1:Kij){
               mixmodj=mixmod[mixmod[,4]==miss[2],]
               muloc=unlist(mixmodj[,2])[k%%nrow(mixmodj)]#ne marche pas, d'abord tout dilater...
               #calcul de la covariance
               #produit par l'inverse de la variance
               #produit par la covariance
               #total
               muloc=muloc
               loc=loc+pi_ij[k]*muloc
            }
            X[miss[1],miss[2]]=(1/Kij)*loc
         }
      }
   }
   return(X)
}