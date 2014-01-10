# 'B matrice p+1 x p   
# ' Sigma vecteur de taille p2 des sigma ssreg
MakeJ<-function(X=X,Z=Z,B=B,Sigma=Sigma,A=A){
   X=cbind(1,X)
   I2=which(colSums(Z)!=0)
   p2=length(I2)
   p1=ncol(X)-p2#prendre donc en compte la constante
   Z=rbind(0,Z)
   Z[1,I2]=1#on ajoute une constante à chaque ssreg
   Z=cbind(0,Z)
   I2=I2+1
   pz=sum(Z!=0)
   n=nrow(X)
   J=matrix(0,ncol=(p1+p2+pz),nrow=(p1+p2+pz))
   barZ=which(Z!=0,arr.ind=T)
   for(j in 1:p2){
      I1j=barZ[barZ[,2]==I2[j],1]
      J[pz+p1+j,pz+p1+j]=2*Sigma[j]#bloc J9
      debcolj=nrow(barZ[barZ[,2]<I2[j],])
      colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
      J[pz+p1+j,colonne]=(2/n)*t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]])%*%X[,I1j]#bloc J7
      J[pz+I1j,which(barZ[,2]==I2[j])]=A[j+1]#attention on compte l'intercept #blocJ4
      J[colonne,pz+p1+j]=(-2/(Sigma[j]^3))*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]]) #bloc J3 
      J[which(barZ[,2]==I2[j]),pz+I1j]=A[j+1]#attention on compte l'intercept #bloc J2
      diag(J[I1j,I1j])=(-1/(Sigma[j]^2))*t(X[,I1j])%*%(X[,I1j])#bloc J1
   }    
   return(J)
}