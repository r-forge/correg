
# ' Zc colSums(Z)
missing_penalty<-function(nbclust_vect=nbclust_vect,Z=Z,M=M,n=n,p=p,Zc=Zc){
   penalty=0
   penalty_vect=nbclust_vect*3-1#3 paramètres par classe mixmod mais somme des proportions à 1
   compl_vect=Zc+2#coef + constante+bruit
   used=rep(0,times=p)
   for (i in 1:n){
      used=0*used
      for(j in 1:p){
         if(M[i,j]==0){#variable observée
            if(Zc[j]==0){#variable à droite
               penalty=penalty+penalty_vect[j]
            }else{#variable à gauche
               penalty=penalty+compl_vect[j]#on compte la structure
            }
         }
      }
      quimankdroit=which(Z%*%(-(M[i,]-1))*M[i,]>0)
      penalty=penalty+penalty_vect[quimankdroit]
   }
   #variable à droite observée : mixmod compte à droite mais pas à gauche car sachant
   #variable à droite manquante : mixmod ne compte pas à droite car pas observée mais intervient dans la loi à gauche si la gauche est observée
   
   #vecteur binaire d'utilisation mis à jour à chaque ligne'
   penalty=log(n)*penalty/n#nombre moyen
   return(penalty)
}