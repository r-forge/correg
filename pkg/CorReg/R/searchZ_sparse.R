#' Sparse structure research
#' @param X the dataset
#' @param Zi indices of the rows of the 1
#' @param Zj indices of the columns of the 1
#' @param Si nombre de 1 dans chaque lignes
#' @param Sj nombre de 1 dans chaque colonnes
#' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
#' @param methode_tirage 0:ligne et colonne,-1:seulement la colonne, entier>0:nombre aleatoire de candidats
#' @param methode_BIC 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
#'@param p1max maximum complexity for a regression
#'@param Maxiter number of steps
#'@param plot TRUE: returns for each step the type of move, complexity and BIC
#'@param best TRUE: systematically jumps to the best BIC seen ever when seen (it is stored even if best=FALSE)
#'@param better TRUE: systematically jumps to the best candidate if better than stationnarity (random wheighted jump otherwise)
#'@param random if FALSE:moves only to improve and only to the best
#'@param verbose 0:none, 1:BIC,step and complexity when best BIC found 2:BIC, step, complexity, nb candidates and best candidate when best BIC found
#'@param nb_opt_max stop criterion defining how many times the chain can walk (or stay) on the max found
# ' @param Mixmod
#' @export
#'@return step 0:delete, 1: add, 2: stationnarity
# '
searchZ_sparse<-function(X=X,Zi=NULL,Zj=NULL,Si=NULL,Sj=NULL,Bic_null_vect=NULL,methode_tirage=2,methode_BIC=1,Rmax=5,Maxiter=1,plot=F,best=T,better=F,random=T,verbose=1,nb_opt_max=NULL){
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Bic_null_vect)){
     Bic_null_vect=density_estimation(X=X)$BIC_vect 
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
  res=.Call( "rechercheZ_sparse_relax",X,Zi,Zj,Si,Sj,Bic_null_vect,methode_tirage,methode_BIC,Rmax,Maxiter,plot,best,better,random,verbose,nb_opt_max, PACKAGE = "CorReg")
  return(res)
}