\name{CorReg-package}
\alias{CorReg-package}
\docType{package}
\title{Algorithms for regression with correlated covariates}


\description{
Sequential linear regression based on a structural equation model
    (explicit correlations). It permits to face highly correlated datasets.
}



\author{

Maintainer: Clement THERY <clement.thery@arcelormittal.com>
}

%\references{

%}

\keyword{ package }



\examples{
require(CorReg)
base=generateur_melange(n=25,p=10,ratio=0.4,sigma_sousreg=0,tp1=1,tp2=0,tp3=0)
X_appr=base$X_appr
Y_appr=base$Y_appr
Y_test=base$Y_test
X_test=base$X_test
TrueZ=base$Z#True generative structure (binary matrix)
Trueclust=base$mixmod#classes of the gaussian mixtures (orthogonal covariates only)

mixmod=calcul_BIC_mixmod(X=X_appr,nbclustmax=10,details=TRUE,mclust=TRUE)#mclust=TRUE is faster
BIC_empty_vect=mixmod$BIC_vect
BIC_empty=sum(BIC_empty_vect)

TrueBIC_vect=BicZ(X=X_appr,Z=TrueZ,Bic_null_vect=BIC_empty_vect)
TrueBIC=sum(TrueBIC_vect)

correl=cor(X_appr)#Correlations matrix will weight initialisation

Zini=Winitial(W=correl,X=X_appr)#generates initial Z for the Markov chain
BIC_ini_vect=BicZ(X=X_appr,Z=Zini,Bic_null_vect=BIC_empty_vect)
Bic_ini=sum(BIC_ini_vect)

#one chain
res=rechercheZ(X=X_appr,Z=Zini,Maxiter=100,Bic_null_vect=BIC_empty_vect,methode_tirage=-1)
#multiple chain
res=rechercheZ(X=X_appr,Z=NULL,nbini=20,Maxiter=100,Bic_null_vect=BIC_empty_vect,methode_tirage=-1)

hatZ=res$Zopt
hatBic=res$bic_opt
sum(BicZ(X=X_appr,Z=hatZ,Bic_null_vect=BIC_empty_vect))#same value

readZ(Z=hatZ,crit="R2",X=X_appr,output="all",order=1)# <NA>line : name of subregressed covariate

compZ=comparaison_struct(vraiZ=TrueZ,Zalgo=hatZ)#qualitative comparison
select="lar"#LASSO with the Least Angle Regression algorithm to find the best model
resY=correg(X=X_appr,Y=Y_appr,Z=hatZ,compl=TRUE,expl=TRUE,pred=TRUE,
               select=select,K=10,final=FALSE,retour=TRUE)

MSE_complete=MSE_loc(Y=Y_test,X=X_test,intercept=TRUE,A=resY$compl$A)
MSE_explicative=MSE_loc(Y=Y_test,X=X_test,intercept=TRUE,A=resY$expl$A)
MSE_predictive=MSE_loc(Y=Y_test,X=X_test,intercept=TRUE,A=resY$pred$A)

MSE=data.frame(MSE_complete,MSE_explicative,MSE_predictive)
#MSE
#barplot(as.matrix(MSE),main="MSE on validation dataset", sub=paste("select=",select))
}