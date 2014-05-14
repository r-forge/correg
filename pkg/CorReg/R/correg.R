#' Estimates the response variable using a structure
#'  @import Rcpp
#'  @import Rmixmod
#' @import lars
#' @import elasticnet
#' @import Matrix
#' @import mclust 
#' @import ridge
#' @import MASS
#' @import parcor
#' @import  clere
#' @import spikeslab
#' @import  rpart
#' @import corrplot
#' @useDynLib CorReg
#' @export
#' @param B the (p+1)xp matrix associated to Z and that contains the parameters of the sub-regressions
#' @param lambda parameter for elasticnet (quadratic penalty)
#' @param X the data matrix (covariates) without the intercept
#' @param Y The response variable vector
#' @param Z The structure (adjacency matrix) between the covariates
#' @param compl boolean to decide if the complete modele is computed
#' @param expl boolean to decide if the explicative model is in the output
#' @param pred boolean to decide if the predictive model is computed
#' @param select selection method in ("lar","lasso","forward.stagewise","stepwise", "elasticnet", "NULL","ridge","adalasso","clere","spikeslab")
#' @param criterion the criterion used to compare the models
#' @param K the number of clusters for cross-validation
#' @param groupe a vector to define the groups used for cross-validation (to obtain a reproductible result)
#' @param Amax the maximum number of covariates in the final model
#' @param returning boolean : second predictive step (selection on I1 knowing I2 coefficients)
#' @param final boolean : recompute estimators without selection on the remaining parameters of the predictive model
#' @param X_test validation sample
#' @param Y_test response for the validation sample
#' @param intercept boolean. If FALSE intercept will be set to 0 in each model.
#' @param alpha Coefficients of the explicative model to coerce the predictive step. if not NULL explicative step is not computed.
#' @param prednew alternate optimisation for predictive
#' @param nbalter number of alternance for prednew
#' @param deltamin criterion to stop alternance 
#' @param g number of group of variables for clere
#' @param explnew select the number of sub-regression to take into account (by AIC on the corresponding final model)
#' 
correg<-function (X = X, Y = Y, Z = NULL, B = NULL, compl = TRUE, expl = FALSE, explnew=FALSE,
                pred = FALSE,prednew=FALSE,
                select = "lar",
                criterion = c("MSE", "BIC"),
                X_test = NULL, Y_test = NULL, intercept = TRUE, 
                K = 10, groupe = NULL, Amax = NULL, lambda = 1,returning=FALSE,final=FALSE,
                nbalter=10,deltamin=0.01,alpha=NULL,g=5) 
{
  res = list()
  X = 1*as.matrix(X)
  K = abs(K)
  K = min(K, nrow(X))
  Y=1*as.matrix(Y)
  if (is.null(groupe)) {
    groupe = rep(0:(K - 1), length.out = nrow(as.matrix(X)))
    groupe = sample(groupe)
  }
  select = select[1]
  if(select=="NULL"){
    returning=FALSE
  }
  if(final){pred=TRUE}
if(explnew){compl=TRUE;expl=TRUE}
  criterion = criterion[1]
  if (is.null(Amax)) {
    Amax = ncol(X) + 1
  }
  if (sum(Z) == 0) {
    compl = TRUE
  }
  if(!is.null(alpha)){
     expl=TRUE
  }
  if (compl) {
    if (select == "NULL") {
      res$compl$A = c(OLS(X = X, Y = Y, intercept = intercept)$beta)
    }else if (select != "elasticnet" & select != "ridge" & select != "adalasso" & select!="clere" & select!="spikeslab") {
      lars_compl = lars(x = X, y = Y, type = select, intercept = intercept)
      res$compl = meilleur_lars(lars = lars_compl, X = X, 
                                Y = Y, mode = criterion, intercept = intercept, 
                                K = K, groupe = groupe, Amax = Amax)
    }else if (select=="elasticnet"){
      lars_compl = renet(x = X, y = Y, intercept = intercept, 
                        lambda = lambda)
      names(lars_compl)[4] = "coefficients"
      res$compl = meilleur_lars(lars = lars_compl, X = X, 
                                Y = Y, mode = criterion, intercept = intercept, 
                                K = K, groupe = groupe, Amax = Amax)
    }else if(select=="adalasso"){
       resada=adalasso(X=X,y=Y,k=K)    
       if(intercept){
          if(is.null(resada$intercept.adalasso)){
             resada$intercept.adalasso=0
          }
          res$compl$A=c(resada$intercept.adalasso,resada$coefficients.adalasso)
       }else{
          res$compl$A=c(resada$coefficients.adalasso)
       }
       Xloc=X[,resada$coefficients.adalasso!=0]
       res$compl$A[res$compl$A!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)
       res$compl$A=c(res$compl$A)
    }else if(select=="clere"){
       res$compl$A=A_clere(y=as.numeric(Y),x=X,g=g)
    }else if(select=="spikeslab"){
       respike=spikeslab(x=X,y=Y)
       res$compl$A=rep(0,times=ncol(X)+intercept)
       if(intercept){
          res$compl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }else{
          res$compl$A[respike$gnet.scale!=0]=OLS(X=X[,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }
   }else{#ridge
      res_ridge = linearRidge(Y~.,data=data.frame(X))
      res$compl$A=coef(res_ridge)
    }
    res$compl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                             beta = res$compl$A)
    res$compl$AIC=mon_AIC(theta=res$compl$A,Y=Y,X=X,intercept=intercept) 
  }
#explicatif
  if (sum(Z) != 0 & (expl | pred)) {
    qui = WhoIs(Z = Z)
    I1 = qui$I1
    I2 = qui$I2
    if(!is.null(alpha)){
       res$expl$A=alpha
    }else if (select == "NULL") {
      res$expl$A = OLS(X = as.matrix(X[, I1]), Y = Y, intercept = intercept)$beta
    }else if (select != "elasticnet" & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
      lars_expl = lars(x = as.matrix(X[, I1]), y = Y, type = select, 
                       intercept = intercept)
      res$expl = meilleur_lars(lars = lars_expl, X = as.matrix(X[, I1]), Y = Y, 
                               mode = criterion, intercept = intercept, 
                               K = K, groupe = groupe, Amax = Amax)
    }else if (select=="elasticnet"){
      lars_expl = renet(x = as.matrix(X[, I1]), y = Y, intercept = intercept, 
                       lambda = lambda)
      names(lars_expl)[4] = "coefficients"
      res$expl = meilleur_lars(lars = lars_expl, X = as.matrix(X[,I1]), Y = Y, 
                               mode = criterion, intercept = intercept, 
                               K = K, groupe = groupe, Amax = Amax)
    }else if(select=="adalasso"){
       resada=adalasso(X=X[,I1],y=Y,k=K)    
       if(intercept){
          if(is.null(resada$intercept.adalasso)){
             resada$intercept.adalasso=0
          }
          res$expl$A=c(resada$intercept.adalasso,resada$coefficients.adalasso)
       }else{
          res$expl$A=c(resada$coefficients.adalasso)
       }
       Xloc=X[,I1][,resada$coefficients.adalasso!=0]
       res$expl$A[res$expl$A!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)       
    }else if(select=="clere"){
       res$expl$A=A_clere(y=as.numeric(Y),x=X[,I1],g=g)
    }else if(select=="spikeslab"){
       respike=spikeslab(x=X[,I1],y=Y)
       res$expl$A=rep(0,times=ncol(X[,I1])+intercept)
       if(intercept){
          res$expl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }else{
          res$expl$A[respike$gnet.scale!=0]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
       }
    }else{#ridge
      lars_expl = linearRidge(Y~.,data=data.frame(X[,I1]))
      res$expl$A=coef(lars_expl)
    }
    if(is.null(alpha)){
       A_expl = rep(0, times = ncol(X) + intercept)
       if(intercept){
          A_expl[c(intercept, I1 + intercept)] = res$expl$A
       }else{
          A_expl[I1] = res$expl$A
       }
       res$expl$A = A_expl
    }else{
       A_expl=alpha
    }
    res$expl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, beta = A_expl)
    res$expl$AIC=mon_AIC(theta=res$expl$A,Y=Y,X=X,intercept=intercept) 
    
    #nouvel explicatif
    if(explnew){#selection sur explicatif####
      qui = WhoIs(Z = Z)
      I1 = qui$I1
      I2 = qui$I2
      R2_vect=R2Z(Z=Z,X=X,adj=TRUE)#evaluation des R2
      R2_vect=R2_vect[R2_vect!=0]
      #tri par R2 => nouvel I2 ordonné
      I2=I2[order(R2_vect,decreasing=FALSE)]#on commence par garder ce qui est mal expliqué par le reste (donc forte perte et peu de corrélations)
      #initialisation du résultat final (par le modèle complet déjà calculé, + AIC associé)
      res$expl2$A=res$expl$A
      resopt=res$expl$A
      AICopt=mon_AIC(theta=resopt,Y=Y,X=X,intercept=intercept)
      for (iexplnew in 1:length(I2)){
         I1=c(qui$I1,I2[1:iexplnew])#c'est ici que tout se joue
          if(!is.null(alpha)){
             res$expl2$A=alpha
          }else if (select == "NULL") {
             res$expl2$A = OLS(X = as.matrix(X[, I1]), Y = Y, intercept = intercept)$beta
          }else if (select != "elasticnet" & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
             lars_expl = lars(x = as.matrix(X[, I1]), y = Y, type = select, 
                              intercept = intercept)
             res$expl2 = meilleur_lars(lars = lars_expl, X = as.matrix(X[, I1]), Y = Y, 
                                      mode = criterion, intercept = intercept, 
                                      K = K, groupe = groupe, Amax = Amax)
          }else if (select=="elasticnet"){
             lars_expl = renet(x = as.matrix(X[, I1]), y = Y, intercept = intercept, 
                               lambda = lambda)
             names(lars_expl)[4] = "coefficients"
             res$expl2 = meilleur_lars(lars = lars_expl, X = as.matrix(X[,I1]), Y = Y, 
                                      mode = criterion, intercept = intercept, 
                                      K = K, groupe = groupe, Amax = Amax)
          }else if(select=="adalasso"){
             resada=adalasso(X=X[,I1],y=Y,k=K)    
             if(intercept){
                if(is.null(resada$intercept.adalasso)){
                   resada$intercept.adalasso=0
                }
                res$expl2$A=c(resada$intercept.adalasso,resada$coefficients.adalasso)
             }else{
                res$expl2$A=c(resada$coefficients.adalasso)
             }
             Xloc=X[,I1][,resada$coefficients.adalasso!=0]
             res$expl2$A[res$expl2$A!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)       
          }else if(select=="clere"){
             res$expl2$A=A_clere(y=as.numeric(Y),x=X[,I1],g=g)
          }else if(select=="spikeslab"){
             respike=spikeslab(x=X[,I1],y=Y)
             res$expl2$A=rep(0,times=ncol(X[,I1])+intercept)
             if(intercept){
                res$expl2$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
             }else{
                res$expl2$A[respike$gnet.scale!=0]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
             }
          }else{#ridge
             lars_expl = linearRidge(Y~.,data=data.frame(X[,I1]))
             res$expl2$A=coef(lars_expl)
          }
          if(is.null(alpha)){
             A_expl = rep(0, times = ncol(X) + intercept)
             if(intercept){
                A_expl[c(intercept, I1 + intercept)] = res$expl2$A
             }else{
                A_expl[I1] = res$expl2$A
             }
             res$expl2$A = A_expl
          }else{
             A_expl=alpha
          }
          res$expl2$BIC = BicTheta(X = X, Y = Y, intercept = intercept, beta = A_expl)
         #calcul du AIC
         AICloc=mon_AIC(theta=res$expl2$A,Y=Y,X=X,intercept=intercept)
         if(AICloc<AICopt){
            resopt=res$expl2$A
            AICopt=AICloc
         }
      }
      #on remet les choses à leur place
      res$expl2$A=resopt
      res$expl2$AIC=AICopt
      res$expl2$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                              beta = res$expl2$A)
      qui = WhoIs(Z = Z)
      I1 = qui$I1
      I2 = qui$I2
    }#fin du nouvel explicatif
    if (pred & length(I2)>0) {
      if(length(I2)<2 & select!="NULL"){
         select="lar"
      }
      
      if (is.null(B)) {
        B = hatB(Z = Z, X = X)
      }
      Xtilde = X[, I2] - cbind(rep(1, times = nrow(X)), X[, I1]) %*% B[c(1, I1 + 1), I2]
      Xtilde = as.matrix(Xtilde)
      if(intercept){
        Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[-1][I1] - A_expl[1]
      }else{
        Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[I1]
      }
      if (select == "NULL") {
        A_inj = OLS(X = Xtilde, Y = Ytilde, intercept = F)$beta
      }else if (select != "elasticnet"  & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
        lars_inj = lars(x = Xtilde, y = Ytilde, type = select, 
                        intercept = F)
        A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                              Y = Ytilde, mode = criterion, intercept = F, 
                              K = K, groupe = groupe)$A
      }else if (select=="elasticnet") {
        lars_inj = renet(x = Xtilde, y = Ytilde, intercept = F, 
                        lambda = lambda)
        names(lars_inj)[4] = "coefficients"
        A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                              Y = Ytilde, mode = criterion, intercept = F, 
                              K = K, groupe = groupe)$A
      }else if(select=="adalasso"){
         resada=adalasso(X=Xtilde,y=Ytilde,k=K)    
         A_inj=c(resada$coefficients.adalasso)
         if(length(A_inj[A_inj!=0])>0){
            Xloc=Xtilde[,A_inj!=0]
            A_inj[A_inj!=0]=c(OLS(X=Xloc,Y=Ytilde,intercept=F)$beta)
         }
      }else if(select=="clere"){
         A_inj=A_clere(y=as.numeric(Ytilde),x=Xtilde,g=g)
         A_inj=A_inj[-1]#vraiment pas propre
      }else if(select=="spikeslab"){
         respike=spikeslab(x=X[,I1],y=Ytilde)
         A_inj=rep(0,times=ncol(Xtilde))
         A_inj[which(respike$gnet.scale!=0)]=OLS(X=Xtilde[,respike$gnet.scale!=0],Y=as.numeric(Ytilde),intercept=intercept)$beta 
      }else{#ridge
        ridge_pred = linearRidge(Ytilde~0+.,data=data.frame(Xtilde))
        A_inj=coef(ridge_pred)
      }
      A_pred = rep(0, times = ncol(X) + intercept)
      A_pred[I2 + intercept] = A_inj
      if(returning){
        Ytildebis=Y-as.matrix(X[,I2])%*%A_pred[I2 + intercept]
        Ytildebis=as.matrix(Ytildebis)
        if (select != "elasticnet" & select != "ridge" & select != "adalasso"  & select!="clere" & select!="spikeslab") {
          lars_retour=lars(x = as.matrix(X[,I1]), y = Ytildebis, type = select, 
                         intercept = intercept)
          A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                   Y = Ytildebis, mode = criterion, intercept = intercept, 
                                   K = K, groupe = groupe)$A
        }else if (select=="elasticnet"){
          lars_retour= renet(x = as.matrix(X[,I1]), y = Ytildebis, intercept =intercept, 
                  lambda = lambda)
          names(lars_retour)[4] = "coefficients"
          A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                   Y = Ytildebis, mode = criterion, intercept = intercept, 
                                   K = K, groupe = groupe)$A
        }else if(select=="adalasso"){
           resada=adalasso(X=X[,I1],y=Ytildebis,k=K)    
           if(intercept){
              A_retour=c(resada$intercept.adalasso,resada$coefficients.adalasso)
           }else{
              A_retour=c(resada$coefficients.adalasso)
           }
           Xloc=X[,I1][,resada$coefficients.adalasso!=0]
           A_retour[A_retour!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)   
         }else if(select=="clere"){
            A_retour=A_clere(y=as.numeric(Y),x=X[,I1],g=g)
         }else if(select=="spikeslab"){
            respike=spikeslab(x=X[,I1],y=Y)
            res$compl$A=rep(0,times=ncol(X[,I1])+intercept)
            if(intercept){
               res$compl$A[c(1,1+which(respike$gnet.scale!=0))]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
            }else{
               res$compl$A[respike$gnet.scale!=0]=OLS(X=X[,I1][,respike$gnet.scale!=0],Y=as.numeric(Y),intercept=intercept)$beta
            }
         }else{#ridge
          ridge_pred = linearRidge(Ytildebis~.,data=data.frame(X[,I1]))
          A_retour=coef(ridge_pred)
        }
        
        if(intercept){
          A_pred[c(intercept, I1 + intercept)] =A_retour
        }else{
          A_pred[ I1 ] = A_retour
        }
      }else{
        if(intercept){
          A_pred[c(intercept, I1 + intercept)] = A_expl[c(intercept,I1 + intercept)] - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }else{
          A_pred[ I1 ] = A_expl[I1] - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }
        
      }
         
      res$pred$A = A_pred
      res$pred$CVMSE = CVMSE(X = X[, which(A_pred[-1] != 0)], Y = Y, intercept = intercept, K = K, groupe = groupe)
      res$pred$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                              beta = A_pred)
      res$pred$AIC=mon_AIC(theta=res$pred$A,Y=Y,X=X,intercept=intercept) 
      
    }
    #nouveau prédictif####
    if (prednew) {
       if(length(I2)<2 & select=="adalasso"){
          select="lar"
       }
       if (is.null(B)) {
          B = hatB(Z = Z, X = X)
       }
       nbit=0
       deltaobs=deltamin+1
       A_old=rep(0,times=(ncol(Z)+intercept))
       while(nbit<nbalter & deltaobs>deltamin){
          nbit=nbit+1
          #on fait tout comme d'hab
             Xtilde = X[, I2] - cbind(rep(1, times = nrow(X)), X[, I1]) %*% B[c(1, I1 + 1), I2]
             Xtilde = as.matrix(Xtilde)
             if(intercept){
                Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[-1][I1] - A_expl[1]
             }else{
                Ytilde = Y - as.matrix(X[, I1]) %*% A_expl[I1]
             }
             if (select == "NULL") {
                A_inj = OLS(X = Xtilde, Y = Ytilde, intercept = F)$beta
             }else if (select != "elasticnet"  & select != "ridge" & select != "adalasso" & select!="clere" & select!="spikeslab") {
                lars_inj = lars(x = Xtilde, y = Ytilde, type = select, 
                                intercept = F)
                A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                                      Y = Ytilde, mode = criterion, intercept = F, 
                                      K = K, groupe = groupe)$A
             }else if (select=="elasticnet") {
                lars_inj = renet(x = Xtilde, y = Ytilde, intercept = F, 
                                 lambda = lambda)
                names(lars_inj)[4] = "coefficients"
                A_inj = meilleur_lars(lars = lars_inj, X = Xtilde, 
                                      Y = Ytilde, mode = criterion, intercept = F, 
                                      K = K, groupe = groupe)$A
             }else if(select=="adalasso"){
                resada=adalasso(X=Xtilde,y=Ytilde,k=K)    
                A_inj=c(resada$coefficients.adalasso)
                if(length(A_inj[A_inj!=0])>0){
                   Xloc=Xtilde[,A_inj!=0]
                   A_inj[A_inj!=0]=c(OLS(X=Xloc,Y=Y,intercept=F)$beta)
                }  
             }else{#ridge
                ridge_pred = linearRidge(Ytilde~0+.,data=data.frame(Xtilde))
                A_inj=coef(ridge_pred)
             }
             A_pred = rep(0, times = ncol(X) + intercept)
             A_pred[I2 + intercept] = A_inj
             if(returning){
                Ytildebis=Y-as.matrix(X[,I2])%*%A_pred[I2 + intercept]
                Ytildebis=as.matrix(Ytildebis)
                if (select != "elasticnet" & select != "ridge" & select != "adalasso" & select!="clere" & select!="spikeslab") {
                   lars_retour=lars(x = as.matrix(X[,I1]), y = Ytildebis, type = select, 
                                    intercept = intercept)
                   A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                            Y = Ytildebis, mode = criterion, intercept = intercept, 
                                            K = K, groupe = groupe)$A
                }else if (select=="elasticnet"){
                   lars_retour= renet(x = as.matrix(X[,I1]), y = Ytildebis, intercept =intercept, 
                                      lambda = lambda)
                   names(lars_retour)[4] = "coefficients"
                   A_retour = meilleur_lars(lars = lars_retour, X = as.matrix(X[,I1]), 
                                            Y = Ytildebis, mode = criterion, intercept = intercept, 
                                            K = K, groupe = groupe)$A
                }else if(select=="adalasso"){
                   resada=adalasso(X=X[,I1],y=Ytildebis,k=K)    
                   if(intercept){
                      A_retour=c(resada$intercept.adalasso,resada$coefficients.adalasso)
                   }else{
                      A_retour=c(resada$coefficients.adalasso)
                   }
                   Xloc=X[,I1][,resada$coefficients.adalasso!=0]
                   A_retour[A_retour!=0]=c(OLS(X=Xloc,Y=Y,intercept=intercept)$beta)   
                }else{#ridge
                   ridge_pred = linearRidge(Ytildebis~.,data=data.frame(X[,I1]))
                   A_retour=coef(ridge_pred)
                }
                
                if(intercept){
                   A_pred[c(intercept, I1 + intercept)] =A_retour
                }else{
                   A_pred[ I1 ] = A_retour
                }
             }else{
                if(intercept){
                   A_pred[c(intercept, I1 + intercept)] = A_expl[c(intercept,I1 + intercept)]- B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
                }else{
                   A_pred[ I1 ] = A_expl[I1]-B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
                }
                
             }
             
             res$prednew$A = A_pred
             res$prednew$CVMSE = CVMSE(X = X[, which(A_pred[-1] != 0)], Y = Y, intercept = intercept, K = K, groupe = groupe)
             res$prednew$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                                     beta = A_pred)
#              if(final){
#                 if(intercept){
#                    quifinal=which(A_pred[-1]!=0)
#                 }else{
#                    quifinal=which(A_pred!=0)
#                 }
#                 Zfinal=as.matrix(Z[quifinal,quifinal])
#                 A_final=correg(X=X[,quifinal],Y=Y,Z=Zfinal,B=as.matrix(B[c(1,quifinal+1),quifinal]),returning=F,final=F,groupe=groupe,K=K,intercept=intercept,criterion=criterion,select="NULL")$pred$A
#                 
#                 res$finalnew$A=res$prednew$A
#                 res$finalnew$A[res$prednew$A!=0]=A_final
#                 res$finalnew$CVMSE = CVMSE(X = as.matrix(X[,quifinal]), Y = Y, intercept = intercept, K = K, groupe = groupe)
#                 res$finalnew$BIC = BicTheta(X = as.matrix(X[,quifinal]), Y = Y, intercept = intercept, 
#                                          beta = A_final)
#              }
          #une fois A calculé, on regarde si on doit continuer
          deltaobs=sum(abs(A_old-res$prednew$A))
          print(deltaobs)
          if(!is.nan(deltaobs)){
             if(deltaobs>deltamin & nbit<nbalter){#on continue
                #on estime un nouveau B
                A_old=res$prednew$A
                newB=newhatB(Z=Z,A=res$prednew$A,Atilde=res$expl$A,Bold=B,intercept=intercept)
                if(prod(B==newB)){#rien n'a change dans B
                   nbit=nbalter+1
                   print("same B : stop")
                }else if(any(newB==Inf)){
                   nbit=nbalter+1
                   print("numerical explosion of B: stop")
                   res$prednew$A=A_old
                }else{
                   B=newB;print(paste("new B",nbit))
                }
             }
          }else{
             nbit=nbalter+1
             print("numerical explosion : stop")
             res$prednew$A=A_old
          }
          
          
       } #fin de la boucle d'optimisation alternee 
       res$B=B
    }
  }else if (sum(Z) == 0 & (expl | pred)) {
    res$expl$A = res$compl$A
    if (pred) {
      res$pred$A = res$compl$A
    }
    if (explnew){
       res$expl2$A = res$compl$A
    }
  }
   if(pred & final){
      res$final$A=res$pred$A
      A_pred=res$pred$A
      if(intercept){
         quifinal=which(A_pred[-1]!=0)
      }else{
         quifinal=which(A_pred!=0)
      }
      if(length(quifinal)==0){
         res$final=res$pred
      }else{
         resA_final=correg(X=X[,quifinal],Y=Y,returning=F,final=F,groupe=groupe,K=K,intercept=intercept,criterion=criterion,select=select,compl=TRUE,expl=FALSE,pred=FALSE,explnew=FALSE)
         quifinal=which(A_pred!=0)
         res$final$A[quifinal]=resA_final$compl$A
         res$final$BIC = resA_final$compl$BIC
      }
      res$final$AIC=mon_AIC(theta=res$final$A,Y=Y,X=X,intercept=intercept) 
      
   }
  return(res)
}