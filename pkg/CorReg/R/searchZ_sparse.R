#' Sparse structure research
#' @param X the dataset
#' @param Zi indices of the rows of the 1
#' @param Zj indices of the columns of the 1
#' @param Si nombre de 1 dans chaque lignes
#' @param Sj nombre de 1 dans chaque colonnes
#' @param bic_vide_vect vecteur BIC de la matrice nulle
#' @param methode_tirage 0:ligne et colonne,-1:seulement la colonne, entier>0:nombre aleatoire de candidats
#' @param methode_BIC 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
#' @param Rmax complexite maximum d'une sous-regression
#' @param Maxiter nombre d'etapes
#' @param plot T:retourne le type de changement, la complexite et le BIC de chaque etapes
#' @param best T:permet d'aller systematiquement au meilleur BIC si il est meilleur que tout les autres deja rencontres 
#' @param better T:permet d'aller systematiquement au meilleur BIC si il est meilleur que l'etape precedente
#' @param random F:permet de s'ameliorer ou de rester sur place
#' @param bla 0:pas de messages, 1:affiche le BIC,le numero d'etape et la complexite de Z quand il y'a un meilleur BIC, 2:affiche le BIC,le numero d'etape,la complexite de Z,le nombre de candidats et le BIC minimum observe parmi les candidats quand il y'a un meilleur BIC, 3: affiche en plus de bla=1 la complexite locale et le BIC local
#' @param nb_opt_max convergence criterion (number of step at the optimum)
# ' @param Mixmod
#' @export
#' @return etape 0:suppression,1 ajout,2 stationarite
# '
searchZ_sparse<-function(X=X,Zi=NULL,Zj=NULL,Si=NULL,Sj=NULL,bic_vide_vect=bic_vide_vect,methode_tirage=2,methode_BIC=1,Rmax=5,Maxiter=1,plot=F,best=T,better=F,random=T,bla=1,nb_opt_max=NULL){
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Zi) | is.null(Zj) ){
     Zi=as.vector(0)
     Zj=as.vector(0)
     p=ncol(X)
     Si=rep(0,times=p)
     Sj=rep(0,times=p)
  }
  if(is.null(Si) | is.null(Sj)){
     p=ncol(X)
     Si=rep(0,times=p)
     Sj=rep(0,times=p)
     for(i in 1:p){
        Si[i]=length(which(Zi==i))
        Sj[i]=length(which(Zj==i))
     }
  }
  res=.Call( "rechercheZ_sparse_relax",X,Zi,Zj,Si,Sj,bic_vide_vect,methode_tirage,methode_BIC,Rmax,Maxiter,plot,best,better,random,bla,nb_opt_max, PACKAGE = "CorReg")
  return(res)
}