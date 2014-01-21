#' cleanZ (if BIC improved)
#' @export
#'@param X the dataset
#'@param Z matrice Z a nettoyer
#'@param bic_vide_vect vecteur BIC de la matrice nulle
#'@param methode_BIC 1:utilisation de la fonction householderQr, 2:utilisation de la fonction colPivHouseholderQr
#'@param star to use BIC*
#'@param plot T:retourne le BIC de chaque etapes
#'@param bla 0:pas de messages, 1:affiche le BIC,le numero d'etape et la complexite de Z quand il y'a un meilleur BIC, 2:affiche le BIC,le numero d'etape,la complexite de Z,le nombre de candidats et le BIC minimum observe parmi les candidats quand il y'a un meilleur BIC
#'@return etape 0:suppression,1 ajout,2 stationarite
#'
cleanZ<-function(X=X,Z=Z,bic_vide_vect=bic_vide_vect,methode_BIC=1,plot=F,bla=1,star=FALSE){
   res=.Call( "cleancolZ",X,Z,bic_vide_vect,methode_BIC,plot,bla, PACKAGE = "CorReg")
  res=.Call( "cleanZ",X,res$Z,bic_vide_vect,methode_BIC,plot,bla, PACKAGE = "CorReg")
  return(res)
  
}