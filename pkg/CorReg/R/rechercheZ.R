#' recherche de structure
#'@param X the dataset
#'@param Z Z est une matrice nulle si on ne lui a pas mis de valeur
#'@param Bic_null_vect vecteur BIC de la matrice nulle
#'@param methode_tirage 0:ligne et colonne,-1:seulement la colonne, entier>0:nombre aleatoire de candidats, -2 : tout le monde (hors diagonale), -3 : uniquement les non-nuls
#'@param reject 0:mode relax, 1:mode reject
#'@param methode_BIC 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
#'@param Rmax complexite maximum d'une sous-regression
#'@param p2max nombre maximal de sous-regressions
#'@param Maxiter nombre d'etapes
#'@param plot T:retourne le type de changement, la complexite et le BIC de chaque etapes
#'@param best T:permet d'aller systematiquement au meilleur BIC si il est meilleur que tout les autres deja rencontres 
#'@param better T:permet d'aller systematiquement au meilleur BIC si il est meilleur que l'etape precedente
#'@param random F:permet de s'ameliorer ou de rester sur place
#'@param bla 0:pas de messages, 1:affiche le BIC,le numero d'etape et la complexite de Z quand il y'a un meilleur BIC, 2:affiche le BIC,le numero d'etape,la complexite de Z,le nombre de candidats et le BIC minimum observe parmi les candidats quand il y'a un meilleur BIC, 3: affiche en plus de bla=1 la complexite locale et le BIC local
#'@param nb_opt_max stop criterion defining how many times the chain can stay at the max found
#'@param exact boolean. If exact subregression is found it gives its content.
#'@param nbini Number of initialisations (using Winitial). if NULL and Zini is NULL : only one chain beginning with zero matrix.
#'@param star boolean to compute BIC* instead of BIC (stronger penalization of the complexity)
#'@param ... parameters to be passed (for Winitial).
#'@return etape 0:suppression,1 ajout,2 stationarite
#'@export
rechercheZ<-function(X=X,Z=NULL,Bic_null_vect=NULL,methode_tirage=-1,reject=0,methode_BIC=1,Rmax=5,p2max=NULL,Maxiter=1,plot=FALSE,best=TRUE,better=FALSE,random=TRUE,bla=1,nb_opt_max=NULL,exact=TRUE,nbini=NULL,star=TRUE,...){
  params=match.call()
  if(is.null(p2max)){
    p2max=ncol(X)+1 
  }
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Z)){
    Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
  }
  if(is.null(Bic_null_vect)){
     Bic_null_vect=calcul_BIC_mixmod(X=X,nbclustmax=10,bla=FALSE,details=FALSE,mclust=TRUE)$BIC_vect
  }
  if(is.null(nbini)){
     if(reject==0){#relax mode
        res=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }else{# reject mode
        res=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }
  }else{
     BICnull=sum(Bic_null_vect)
     if(!("W" %in% names(params))){
        W=cor(X)
     }
     res=list()
     if(nbini>1){#first try with zero matrix
        if(reject==0){#relax mode
           resloc=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           resloc=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        if(resloc$bic_opt<=min(res$bic_opt,BICnull)){
           res=resloc
        }
        nbini=nbini-1
     }
     for(i in 1:nbini){
        Z=Winitial(W=W,X=X,Rmax=Rmax,Bic_null_vect=Bic_null_vect,p2max=p2max)
        if(reject==0){#relax mode
           resloc=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           resloc=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,methode_tirage,methode_BIC,Rmax,p2max,Maxiter,plot,best,better,random,bla,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        if(resloc$bic_opt<=min(res$bic_opt,BICnull)){
           res=resloc
        }
     }
     return(res)
  }

}