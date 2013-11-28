# ' generateur de param?tres selon la structure
# ' @param pb kind of problem : 0=none, 1=simple, 2=strong
generateurA_ou<-function (Z = Z, tp1 = 1, tp2 = 1, tp3 = 1, positif = 0.6, lambdapois = 5,pb=2,Amax=NULL,B=NULL) 
{
  Z = as.matrix(Z)
  p = ncol(Z)
  quip1 = which(rowSums(Z) != 0)
  quip2 = which(colSums(Z) != 0)
  quip3 = (1:p)[-c(quip1, quip2)]
  A = 1+rpois(p + 1, lambdapois) * (rep(-1, p + 1) + 2 * rbinom(p +1, 1, positif))
  if(pb>0 & !is.null(B)){
    Apb=as.vector((B-diag(diag(as.matrix(B))))%*%A)#on fait une combinaison linéaire des sous-régressions pour mettre le lasso en défaut 
    qui0=which(Apb==0)
    Apb[qui0]=A[qui0]
    A=round(Apb)#on impose des valeurs entières pour plus de facilité lors de la comparaison des résultats
  }
  A[-1][quip2]=A[-1][quip2]*rbinom(length(quip2),1,tp2)#on considère le taux comme une proba
  A[-1][quip3]=A[-1][quip3]*rbinom(length(quip3),1,tp3)#on considère le taux comme une proba
  
  if(pb<2){
      A[-1][quip1]=A[-1][quip1]*rbinom(length(quip1),1,tp1)#on considère le taux comme une proba
  }else{#on force le probleme 
    quip1order=unique(which(Z!=0,arr.ind=T)[,1])#dans l'ordre d'apparition
    nb0p1=round((1 - tp1) * length(quip1))
    if(nb0p1>0){#si on doit supprimer des points
      A[-1][quip1order][1:nb0p1] = 0
    }
  }
  A1=length(which(A!=0))
  if(!is.null(Amax) & Amax<A1){#si on doit mettre une borne à Amax
    nb0=A1-Amax#nombre de 0 à placer
    if(pb<2){
      A[A!=0][sample(A1,size=nb0)]=0
    }else{#on conserve les pb
      #mise à jour des indices des survivants
      priorite=1+c(quip3,0,quip2,quip1order)
      priorite=priorite[A[priorite]!=0]
      A[priorite][1:nb0]=0
    } 
  }
  return(A)
}