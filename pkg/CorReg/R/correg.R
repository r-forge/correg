#' Estimates the response variable using a structure
#' @useDynLib CorReg
#' @export
#' @param B the (p+1)xp matrix assiated to Z and that contains the parameters of the sub-regressions
#' @param lambda parameter for elasticnet (quadratic penalty)
#' @param X the data matrix (covariates) without the intercept
#' @param Y The response variable vector
#' @param Z The structure (adjacency matrix) between the covariates
#' @param compl boolean to decide if the complete modele is computed
#' @param expl boolean to decide if the explicative model is in the output
#' @param pred boolean to decide if the predictive model is computed
#' @param select selection method in ("lar","lasso","forward.stagewise","stepwise", "elasticnet", "NULL","ridge")
#' @param criterion the criterion used to compare the models
#' @param K the number of clusters for cross-validation
#' @param groupe a vector to define the groupes used for cross-validation (to obtain a reproductible result)
#' @param Amax the maximum number of covariates in the final model
#' @param retour boolean : second predictive step
#' @param final boolean : recompute estimators without selection on the remaining parameters of the predictive model
#' @param X_test validation sample
#' @param Y_test response for the validation sample
#' @param intercept boolean. If FALSE intercept will be set to 0 in each model.
#' 
correg<-function (X = X, Y = Y, Z = NULL, B = NULL, compl = T, expl = T, 
                pred = T, 
                select = "lar",
                criterion = c("MSE", "BIC"),
                X_test = NULL, Y_test = NULL, intercept = T, 
                K = 10, groupe = NULL, Amax = NULL, lambda = 1,retour=T,final=T) 
{
  res = list()
  X = as.matrix(X)
  K = abs(K)
  K = min(K, nrow(X))
  Y=as.matrix(Y)
  if (is.null(groupe)) {
    groupe = rep(0:(K - 1), length.out = nrow(as.matrix(X)))
    groupe = sample(groupe)
  }
  select = select[1]
  if(select=="NULL"){
    retour=F
    final=F
  }
  criterion = criterion[1]
  if (is.null(Amax)) {
    Amax = ncol(X) + 1
  }
  if (sum(Z) == 0) {
    compl = T
  }
  if (compl) {
    if (select == "NULL") {
      res$compl$A = c(OLS(X = X, Y = Y, intercept = intercept)$beta)
    }else if (select != "elasticnet" & select != "ridge") {
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
    }else{#ridge
      res_ridge = linearRidge(Y~.,data=data.frame(X))
      res$compl$A=coef(res_ridge)
    }
    res$compl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                             beta = res$compl$A)
  }
  if (sum(Z) != 0 & (expl | pred)) {
    qui = WhoIs(Z = Z)
    I1 = qui$I1
    I2 = qui$I2
    if (select == "NULL") {
      res$expl$A = OLS(X = as.matrix(X[, I1]), Y = Y, intercept = intercept)$beta
    }else if (select != "elasticnet" & select != "ridge" ) {
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
    }else{#ridge
      lars_expl = linearRidge(Y~.,data=data.frame(X[,I1]))
      res$expl$A=coef(lars_expl)
    }
    A_expl = rep(0, times = ncol(X) + intercept)
    if(intercept){
      A_expl[c(intercept, I1 + intercept)] = res$expl$A
    }else{
      A_expl[I1] = res$expl$A
    }
    res$expl$A = A_expl
    res$expl$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                            beta = A_expl)
    if (pred) {
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
      }else if (select != "elasticnet"  & select != "ridge") {
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
      }else{#ridge
        ridge_pred = linearRidge(Ytilde~0+.,data=data.frame(Xtilde))
        A_inj=coef(ridge_pred)
      }
      A_pred = rep(0, times = ncol(X) + intercept)
      A_pred[I2 + intercept] = A_inj
      if(retour){
        Ytildebis=Y-as.matrix(X[,I2])%*%A_pred[I2 + intercept]
        Ytildebis=as.matrix(Ytildebis)
        if (select != "elasticnet" & select != "ridge") {
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
          A_pred[c(intercept, I1 + intercept)] = A_expl[c(intercept,I1 + intercept)] 
          - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }else{
          A_pred[ I1 ] = A_expl[I1] 
          - B[c(1, I1 + 1), I2] %*% as.matrix(A_inj)
        }
        
      }
      
      res$pred$A = A_pred
      res$pred$CVMSE = CVMSE(X = X[, which(A_pred[-1] != 
                                             0)], Y = Y, intercept = intercept, K = K, groupe = groupe)
      res$pred$BIC = BicTheta(X = X, Y = Y, intercept = intercept, 
                              beta = A_pred)
      if(final){
        if(intercept){
          quifinal=which(A_pred[-1]!=0)
        }else{
          quifinal=which(A_pred!=0)
        }
        Zfinal=as.matrix(Z[quifinal,quifinal])
        A_final=correg(X=X[,quifinal],Y=Y,Z=Zfinal,B=as.matrix(B[c(1,quifinal+1),quifinal]),retour=F,final=F,groupe=groupe,K=K,intercept=intercept,criterion=criterion,select="NULL")$pred$A
        
        res$final$A=res$pred$A
        res$final$A[res$pred$A!=0]=A_final
        res$final$CVMSE = CVMSE(X = as.matrix(X[,quifinal]), Y = Y, intercept = intercept, K = K, groupe = groupe)
        res$final$BIC = BicTheta(X = as.matrix(X[,quifinal]), Y = Y, intercept = intercept, 
                                 beta = A_final)
        
      }
      
    }
    
  }else if (sum(Z) == 0 & (expl | pred)) {
    res$expl$A = res$compl$A
    if (pred) {
      res$pred$A = res$compl$A
    }
  }
  return(res)
}