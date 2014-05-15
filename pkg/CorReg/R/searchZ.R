#' MCMC algo to find a structure between the covariates
#'@param X the dataset
#'@param Z binary adjacency matrix of size p. if NULL zero matrix is used
#' @param Bic_null_vect the BIC of the null hypothesis (used for independent variables)
#'@param candidates 0:row and column,-1:column only, int>0:random int candidates, -2 : all (but the diag), -3 : non-zeros
#'@param reject 0: constraint relaxation, 1: reject mode
#' @param methode  parameter for OLS (matrix inversion) 1:householderQr, 2:colPivHouseholderQr
#'@param p1max maximum complexity for a regression
#'@param p2max maximum number of regressions 
#'@param Maxiter number of steps
#'@param plot TRUE: returns for each step the type of move, complexity and BIC
#'@param best TRUE: systematically jumps to the best BIC seen ever when seen (it is stored even if best=FALSE)
#'@param better TRUE: systematically jumps to the best candidate if better than stationnarity (random wheighted jump otherwise)
#'@param random if FALSE:moves only to improve and only to the best
#'@param verbose 0:none, 1:BIC,step and complexity when best BIC found 2:BIC, step, complexity, nb candidates and best candidate when best BIC found
#'@param nb_opt_max stop criterion defining how many times the chain can walk (or stay) on the max found
#'@param exact boolean. If exact subregression is found it gives its content (another verbose mode).
#'@param nbini Number of initialisations (using Winitial if Z is NULL). if NULL and Zini is NULL : only one chain beginning with zero matrix.
#'@param star boolean to compute BIC* instead of BIC (stronger penalization of the complexity). WARNING : star=TRUE implies p2max<=p/2
#'@param clean cleaning steps at the end of the walk (testing each remainging 1 for removal)
#'@param ... parameters to be passed (for Winitial).
#'@return step 0:delete, 1: add, 2: stationnarity
#'@export
searchZ<-function(X=X,Z=NULL,Bic_null_vect=NULL,candidates=-1,reject=0,methode=1,p1max=5,p2max=NULL,Maxiter=1,plot=FALSE,best=TRUE,better=FALSE,random=TRUE,verbose=1,nb_opt_max=NULL,exact=TRUE,nbini=NULL,star=TRUE,clean=TRUE,...){
  params=match.call()
  Wini=FALSE
  X=1*as.matrix(X)
  if(is.null(p2max)){
    p2max=ncol(X)+1 
  }
  if(star){
     p2max=floor(min(p2max,ncol(X)/2))
     p1max=floor(min(p1max,ncol(X)/2))
  }
  if(is.null(nb_opt_max)){
    nb_opt_max=Maxiter
  }
  if(is.null(Bic_null_vect)){
     Bic_null_vect=density_estimation(X=X,nbclustmax=10,verbose=FALSE,detailed=FALSE,mclust=TRUE)$BIC_vect
  }
  if(!is.null(nbini)){
     if (nbini<1){nbini=NULL}
  }
  if(is.null(nbini)){
     if(is.null(Z)){
        Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
     }
     if(reject==0){#relax mode
        res=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }else{# reject mode
        res=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        return(res)
     }
  }else{
     BICnull=sum(Bic_null_vect)
     if(!("W" %in% names(params))){
        W=cor(X)
     }
     if(is.null(Z)){
        Z=matrix(0,ncol=ncol(X),nrow=ncol(X))
        Wini=TRUE
     }
     res=list()
     if(nbini>1){#first try with provided Z matrix (or null if not provided)
        if(reject==0){#relax mode
           res=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           res=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        nbini=nbini-1#to finally have the exact number of tries
        if(clean){
           resclean=cleanZ(X=X,Z=res$Z_opt,Bic_null_vect=Bic_null_vect,star=star,verbose=verbose)#nettoyage colonnes puis ponctuel
           res$Z_opt=resclean$Z_opt
           res$bic_opt=resclean$bic_opt
        }
     }
    
     for(i in 1:nbini){
        if(Wini){#only if no Z provided
           Z=Winitial(W=W,X=X,p1max=p1max,Bic_null_vect=Bic_null_vect,p2max=p2max)
        }
        if(reject==0){#relax mode
           resloc=.Call( "rechercheZ_relax",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }else{# reject mode
           resloc=.Call( "rechercheZ_rejet",X,Z,Bic_null_vect,candidates,methode,p1max,p2max,Maxiter,plot,best,better,random,verbose,nb_opt_max,exact,star, PACKAGE = "CorReg")
        }
        if(length(res)==0){res=resloc}
        if(resloc$bic_opt<=min(res$bic_opt,BICnull)){
           res=resloc
        }
     }
     if(clean){
        resclean=cleanZ(X=X,Z=res$Z_opt,Bic_null_vect=Bic_null_vect,star=star,verbose=verbose)#nettoyage colonnes puis ponctuel
        res$Z_opt=resclean$Z_opt
        res$bic_opt=resclean$bic_opt
     }
     return(res)
  }

}