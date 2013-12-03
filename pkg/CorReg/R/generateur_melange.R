#' generateur de donnees en sous-regression par modele de melange gaussien
#utilise Matrix (matrices creuses) et generateurZ
#' @param n the number of individuals in the learning dataset
#' @param p the number of covariates (without the response)
#' @param ratio the ratio of explained covariates (dependent)
#' @param max_compl the number of covariates in each subregression
#' @param valid the size of the validation sample
#' @param positif the ratio of positive coefficients in both the regression and the subregressions
#' @param sigma_Y standard deviation for the noise of the regression
#' @param sigma_sousreg standard deviation for the noise of the subregression (all). ignored if gamma=T
#' @param gamma boolean to generate a p-sized vector sigma_sousreg gamma-distributed
#' @param gammashape shape parameter of the gamma distribution (if needed)
#' @param gammascale scale parameter of the gamma distribution (if needed)
#' @param meanvar vector of means for the covariates.
#' @param sigmavar standard deviation of the covariates.
#' @param lambda paramater of the law that define the number of components in gaussian mixture models
#' @param Amax the maximum number of covariates with non-zero coefficients in the regression
#' @param tp1 the ratio of right-side covariates allowed to have a non-zero coefficient in the regression
#' @param tp2 the ratio of left-side covariates allowed to have a non-zero coefficient in the regression
#' @param tp3 the ratio of strictly independent covariates allowed to have a non-zero coefficient in the regression
#' @param lambdapois parameter used to generate the coefficient in the subregressions. poisson distribution.
#' @param pb Defines an heuristic to generate Y in a way that will give some issues with correlations.
#' @export
generateur_melange<-function(n=130,
                                p=100,
                                ratio=0.4,
                                max_compl=3,
                                valid=1000,
                                positif=0.6,
                                sigma_Y=10,
                                sigma_sousreg=0.25,
                                meanvar=NULL,
                                sigmavar=NULL,
                                lambda=5,#pour l enombre de composantes des m?langes gaussiens
                                Amax=15,
                                lambdapois=5,#pour les valeurs des coefs
                                gamma=T,
                                gammashape=1,
                                gammascale=0.5,tp1=1,tp2=1,tp3=1,pb=2
){
  max_compl=min(max_compl,p)
  Amax=min(p+1,Amax) # min entre p+1 et Amax  why?
  R=round(ratio*p) # R : entier nombre de personne a gauche
  if(R==0){pb=0}
  if(is.null(Amax) | Amax>p){Amax=p}
  #lmabda param?tre le nombre de composantes des m?langes gaussiens   
  valid=max(valid,n)
  B=Matrix(0,nrow=(p+1),ncol=(p+1)) #B matrice creuse
  taille=n+valid
  qui=2:(p+1)  # vecteur de taille p qui contient les valeurs allant de 2 à p+1 avec un pas de 1
  if(R>0){
    list_X2=sample(qui,R) # melange R individus pri au hasard 
    B[1,list_X2]=rpois(R,lambdapois)*(rep(-1,R)+2*rbinom(R,1,positif)) # 1ere ligne de B = ?????
    for(j in list_X2){
      B[sample(qui[-c(list_X2-1)],size=max_compl),j]=(1/max_compl)*rpois(max_compl,5)*(rep(-1,max_compl)+2*rbinom(max_compl,1,positif))
    }
    #ajout de G
    G=Diagonal(p+1)
    G[,list_X2]=0
    B=B+G  
  }
  vraiZ=B
  vraiZ[vraiZ!=0]=1
  vraiZ=as.matrix(vraiZ)
  vraiZ=vraiZ-diag(diag(vraiZ))
  vraiZ=vraiZ[-1,-1]
  if(sum(c(tp1,tp2,tp3))==0){
    A=rpois(p+1,lambdapois)*(rep(-1,p+1)+2*rbinom(p+1,1,positif)) 
    #on vient maintenant tailler dans A en fonction des param?tres
    A[-sample(p+1,size=Amax)]=0
  }else{
    A=generateurA_ou(Z=vraiZ,tp1=tp1,tp2=tp2,tp3=tp3,positif=positif,lambdapois=lambdapois,pb=pb,Amax=Amax,B=B)
  }

  composantes=rpois(p-R,lambda=lambda)
  composantes[composantes>n]=n #pas plus de composantes que de points
  composantes[composantes==0]=1#au moins une composante
  ploc=sum(composantes)
  meanvar=NULL
  sigmavar=NULL
  if(is.null(meanvar)){
    meanvar=rpois(ploc,ploc)*(rep(-1,ploc)+2*rbinom(ploc,1,positif))
  }
  if(is.null(sigmavar)){
    sigmavar=5
  }
  prop=runif(ploc)   
  comp_cum=cumsum(composantes)
  comp_cum=c(0,comp_cum)  
  X=matrix(0,ncol=p+1,nrow=taille)
  epsX=X#matrice nulle
  X1g=matrix(rnorm(taille*ploc,mean=meanvar,sd=sigmavar),ncol=ploc,nrow=taille,byrow=T)
  dim(X1g)
  X1=cbind(rep(1,times=taille),matrix(0,ncol=p-R,nrow=taille))  
  for(i in 1:(p-R)){
    prop[(comp_cum[i]+1):comp_cum[i+1]]=prop[(comp_cum[i]+1):comp_cum[i+1]]/sum(prop[(comp_cum[i]+1):comp_cum[i+1]])
    qui=sample(n)
    quiv=sample((n+1):taille)
    combien=rep(1,times=composantes[i])+floor((n-composantes[i])*prop[(comp_cum[i]+1):comp_cum[i+1]])
    if(sum(combien)<n){
      quiloc=sample(composantes[i],size=(n-sum(combien)))
      combien[quiloc]=combien[quiloc]+1
    }
    combien=c(0,cumsum(combien))
    #idem pour la validation
    combienv=rep(1,times=composantes[i])+floor((valid-composantes[i])*prop[(comp_cum[i]+1):comp_cum[i+1]])
    if(sum(combienv)<valid){
      quiloc=sample(composantes[i],size=(valid-sum(combienv)))
      combienv[quiloc]=combienv[quiloc]+1
    }
    combienv=c(0,cumsum(combienv))
    for (j in 1:composantes[i]){
      X1[qui[(combien[j]+1):combien[j+1]],1+i]=X1g[(combien[j]+1):combien[j+1],comp_cum[i]+j]
      X1[quiv[(combienv[j]+1):combienv[j+1]],1+i]=X1g[(combien[j+1]+combienv[j]+1):(combien[j+1]+combienv[j+1]),comp_cum[i]+j]
    }
  }
  rm(X1g)

  if(R>0){ 
    X[,-list_X2]=X1
    if(gamma){
      sigma_sousreg=rgamma(R,shape=gammashape,scale=gammascale)
    }
    epsX[,list_X2]=matrix(rnorm(taille*R,mean=rep(0,times=R),sd=sigma_sousreg),ncol=R,nrow=taille,byrow=T)
    X=X%*%B+epsX
  }else{
    X=X1
  }  
  X=as.matrix(X)
  #names(X)=c("cste",paste("X_",1:p, sep=""))

  Y=X%*%A+rnorm(taille,mean=0,sd=sigma_Y)
  X_appr=as.matrix(X[1:n,-1])
  X_test=as.matrix(X[(n+1):taille,-1])
  Y_appr=as.matrix(Y[1:n])
  Y_test=as.matrix(Y[(n+1):taille]) 
  return(list(X_appr=X_appr,Y_appr=Y_appr,A=A,B=B,Z=vraiZ,
              X_test=data.frame(X_test),Y_test=Y_test,sigma_sousreg=sigma_sousreg,mixmod=composantes))   
}


