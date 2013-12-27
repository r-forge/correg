
#' @export
newhatB<-function(X,Z,A,Atilde,Bold,intercept=TRUE){
   qui=WhoIs(Z=Z,I3=F,I2=T,I1=T)
   I2=qui$I2
   X=cbind(1,X)
   if(!intercept){
      A=c(0,A)
      Atilde=c(0,Atilde)
   }
   B=Bold
   Atilde=matrix(Atilde,ncol=1);A=matrix(A,ncol=1)
   for(j in I2){
      if(as.numeric(A[j+intercept])!=0){#si A2 est nul,  on ne peut utiliser la formule et donc on garde Bold
         I1loc=c(1,which(Z[,j]!=0)+1)
         if(sum(Atilde[I1loc]!=A[I1loc])!=0){#A et Atilde sont égaux donc A2 est nul, on ne peut utiliser la formule et donc on garde Bold
            Xloc=X[,I1loc]
            B[I1loc,j]=(1/as.numeric(A[j+intercept])(Atilde[I1loc]-A[I1loc]))              
            #print("calcul effectif")
         }
      }
   }
   return(B)
}
