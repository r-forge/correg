#' compare zero in two coefficient vector (same size)
#'@export

comparateur0_vect<-function(vraiA=vraiA,Aalgo=Aalgo,taux=F){
   quivrai0=which(vraiA==0)
   nbbon0=length(which(Aalgo[quivrai0]==0))
   nbfaux0=length(which(Aalgo[-quivrai0]==0))
   if(taux){
      tauxbon0=nbbon0/length(quivrai0)
      tauxfaux0=nbfaux0/length(which(Aalgo==0))
      return(list(nbbon0=nbbon0,nbfaux0=nbfaux0,tauxbon0,tauxfaux0))
   }else{
      return(list(nbbon0=nbbon0,nbfaux0=nbfaux0))    
   }
}