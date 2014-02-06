#' To compare structures (Z)
#' @export
#' @param trueZ first structure
#' @param Zalgo second structure
#' @param all boolean. Also compute the ratio for each stat.
#' @param mode how to modify the structures before comparison. mode=c("NULL","hybrid","clique","sym")
compare_struct<-function(trueZ=trueZ,Zalgo=Zalgo,all=TRUE,mode="NULL"){
  if(mode=="hybrid"){
    trueZ=trueZ+t(trueZ)+t(trueZ)%*%trueZ#attention ? bien multiplier par la transpos?e ? gauche
    trueZ[trueZ>1]=1     
    Zalgo=Zalgo+t(Zalgo)+t(Zalgo)%*%Zalgo#attention ? bien multiplier par la transpos?e ? gauche
    Zalgo[Zalgo>1]=1      
  }else if(mode=="clique"){#on fait des cliques
    trueZ=cliquefaction(trueZ)
    Zalgo=cliquefaction(Zalgo)
  }else if (mode=="sym"){#on fait juste la d?sorientation
    trueZ=trueZ+t(trueZ)
    Zalgo=Zalgo+t(Zalgo)
  }else{
    #on garde les structures d'origine
  }
  trueZ=as.matrix(trueZ)
  Zalgo=as.matrix(Zalgo)
  res=as.matrix(trueZ-Zalgo)
  nbbon1=sum(Zalgo*trueZ)#produit de hadamard
  nbtrop=-sum(res[res==-1])#faux 1
  nbmank=sum(res[res==1]) #faux 0
  nbbon0=ncol(trueZ)^2-nbbon1-nbtrop-nbmank #vrai 0
  #dist=nbtrop+nbmank
  #on passe aux pourcentages
  taux_bon1=nbbon1/sum(trueZ)#taux de v?rit? d?couverte
  taux_bon0=nbbon0/(ncol(trueZ)^2-sum(trueZ))
  taux_faux1=nbtrop/sum(Zalgo) #taux d'ajout dans ce qui est dit doit tendre vers 0 (max si on ne dit que des aneries)
  taux_faux0=1-taux_bon0 #taux d'oublis doit tendre vers 0 (max si on n'a rien dit de vrai)
  vraissreg=which(colSums(trueZ)>0) 
  ssregalgo=which(colSums(Zalgo)>0)
  deltap2=length(vraissreg)-length(ssregalgo)
  bon_gauche=sum(duplicated(c(vraissreg,ssregalgo)))
  faux_gauche=length(ssregalgo)-bon_gauche
  if(all){
    return(list(taux_bon1=taux_bon1,taux_bon0=taux_bon0,taux_faux1=taux_faux1,taux_faux0=taux_faux0,nbbon1=nbbon1,nbtrop=nbtrop,nbmank=nbmank,deltap2=deltap2,bon_gauche=bon_gauche,faux_gauche=faux_gauche))
  }else{
    return(list(verite_decouverte=taux_bon1,taux_faux1=taux_faux1))
  }
}