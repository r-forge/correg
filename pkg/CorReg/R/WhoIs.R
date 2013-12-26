# ' qui est qui  ?
# ' 
# ' @param Z the structure (square matrix)
# ' @param I1 (and others) wanted output
# ' 
#'@export 
WhoIs<-function(Z=Z,I3=F,I2=T,I1=T){
  res=list()
  if(I2){
    res$I2=which(colSums(Z)!=0)# Qui est ? gauche
  }
  if(I1){
    res$I1=which(colSums(Z)==0)#qui n'est pas ? gauche
  }
  if(I3){
    res$I3=which(colSums(Z)==0 & rowSums(Z)==0)# qui est totalement isol?
  }
  return(res)
}