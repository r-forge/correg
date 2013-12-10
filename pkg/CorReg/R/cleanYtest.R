#' Selection method based on p-values (coefficients)
#' @param Y the response variable
#' @param X the dataset of the covariates (without the response)
#' @param pvalmin the maximal bound for p-value associated to remaining coefficients
#' @param bonferroni boolean defining wether a Bonferroni correction is applied or not
#' @param A optional vector of coefficients to coerce some zeros
#' @export
cleanYtest<-function (Y = Y, X = X, pvalmin = 0.05, bonferroni=F,A=NULL) 
{
  p=ncol(X)
  qui=NULL
  if(!is.null(A)){
    qui=which(A[-1]!=0)
    X=X[,qui]
  }
  #on regarde chaque coef donc on boucle jusqu'à stabilité
  pvalminini=pvalmin
  change=T#changement potentiel
  quinonzero=1:(ncol(X))
  loc=length(quinonzero)
  A=rep.int(0, times=ncol(X)+1)
  while(change & ncol(X)>0){
    if(bonferroni){pvalmin=pvalminini/(ncol(X))}
      lmloc=lm(Y~.,data=data.frame(X))
      summar=summary(lmloc)
      coefs_pval=coef(summar)[,4]#p-values des coefficients
      quivarzero=which(coefs_pval[-1]>pvalmin)
    if(length(quivarzero)>0){#on élague juste les coefs pourris
      quinonzero=quinonzero[-quivarzero]
      X=as.matrix(X[,quinonzero])
      loc=length(quinonzero)
    }else{#on n'a rien changé
      change=F
    }
  }
  #on regarde la constante et on l'enlève si besoin
  if(coefs_pval[1]>pvalmin){
    Aloc=lm(Y~0+.,data=data.frame(X))$coefficients
    quinonzero=quinonzero[-1] 
    A[quinonzero+1]=Aloc
  }else{
    Aloc=lmloc$coefficients
  }
  A[c(1,quinonzero)]=Aloc
  if(!is.null(qui)){#on avait un A plus grand
    Along=rep(0,times=(p+1))
    Along[c(1,qui+1)]=A
    return(Along)
  }else{
    return(A)
  }
}

