#' compare les 0 dans des vecteurs
#' @export
comparateur0_vect<-function(vraiA=vraiA,Aalgo=Aalgo,taux=F){
  quivrai0=which(vraiA==0)
  if(length(quivrai0)>0){
    nbbon0=length(which(Aalgo[quivrai0]==0))
    nbbon1=length(which(Aalgo[-quivrai0]!=0))
    nbfaux0=length(which(Aalgo[-quivrai0]==0))
    nb0mank=length(quivrai0)-nbbon0
  }else{
    nbbon0=0
    nbbon1=length(which(Aalgo!=0))
    nbfaux0=length(which(Aalgo==0))
    nb0mank=0
  }
  if(taux){
    tauxbon0=nbbon0/length(quivrai0)
    tauxfaux0=nbfaux0/length(which(Aalgo==0))
    return(list(nbbon0=nbbon0,nbbon1=nbbon1,nbfaux0=nbfaux0,nb0mank=nb0mank,tauxbon0=tauxbon0,tauxfaux0=tauxfaux0))
  }else{
    return(list(nbbon0=nbbon0,nbbon1=nbbon1,nbfaux0=nbfaux0,nb0mank=nb0mank))    
  }
}