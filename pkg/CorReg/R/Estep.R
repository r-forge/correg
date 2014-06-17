#' Imputation of missing values knowing alpha (E step of the EM)
#' 
#' @param Ir the set of indices of the covariates on the left (computed once before EM to be efficient)
#' @param Zc profil (sum) colonne de Z
#' @param alpha matrix of the coefficients (p+1)xp
#' @export
Estep<-function(X=X,alpha=alpha,M=M,Z=Z,mixmod=mixmod,Zc=Zc){
   #X bouge, mais pas alpha ni Z ni M ni mixmod
   #on liste tout ce qu'il faut (on verra après comment extraire ça pour optimiser)
   quimank=which(is.na(X),arr.ind=TRUE)# apriori on stockera les deux formats (M et quimank en permanence) et on fera pareil pour Z
   nbmank=nrow(quimank)
   Zc=colSums(Z)
   #peut-être entrer quimank et écrire uniquement dans le Cpp M au format Eigen creux (sans contact avec R)
   
   #on remplit
   for (miss in 1:nbmank){#pour chaque valeur manquante
      miss=quimank[miss,]
      if(Zc[miss[2]]!=0){#trou à gauche
         #imputation par moyenne du mélange gaussien résultant de l'intégration de la sous-regression
         X[miss[1],miss[2]]=alpha[1,miss[2]]+as.matrix(X[miss[1],-miss[2]])%*%as.matrix(alpha[-1,quimank[miss,2]][-quimank[miss,2]])

      }else{#trou à droite
         
      }
   }
   return(X)
}