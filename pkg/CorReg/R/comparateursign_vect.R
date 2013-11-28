#' compare les signes dans des vecteurs
#' @export
comparateursign_vect<-function(vraiA=vraiA,Aalgo=Aalgo){
  quivrai0=which(vraiA==0)
  nbbon0=length(which(Aalgo[quivrai0]==0))
  nbbon1=length(which(Aalgo[-quivrai0]!=0))
  nbfaux0=length(which(Aalgo[-quivrai0]==0))
  nb0mank=length(quivrai0)-nbbon0
  
  #comptage des signes +
  quivraiplus=which(vraiA>0)
  quivraimoins=which(vraiA<0)
  #attention , on n'a pas une partition (à cause des 0)
  nbbonplus=length(which(Aalgo[quivraiplus]>0))
  nbbonmoins=length(which(Aalgo[quivraimoins]<0))
  nbfauxplus=length(which(Aalgo[quivraimoins]>0))
  nbfauxmoins=length(which(Aalgo[quivraiplus]<0))
 
  return(list(nbbonplus=nbbonplus,nbbonmoins=nbbonmoins,nbfauxplus=nbfauxplus,nbfauxmoins=nbfauxmoins))    
}