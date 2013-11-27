#include "CVMSE.h"
#include "OLS_cpp.h"
#include <math.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
using namespace Rcpp ;
using namespace std;
using Rcpp::as;
using namespace Eigen;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


SEXP CVMSE(SEXP RX,SEXP RY,SEXP RK,SEXP Rintercept,SEXP Rmethode,SEXP Rgroupe)
{
  BEGIN_RCPP
  //déclaration des paramètres
  const Map<MatrixXd> X(as<Map<MatrixXd> >(RX));//on utilise le format EIgen pour la matrice pour bénéficier des outils matriciels de base
  const Map<MatrixXd> Y(as<Map<MatrixXd> >(RY));
  const Map<VectorXi> groupe(as<Map<VectorXi> >(Rgroupe));
  
  int K = Rcpp::as<int>(RK);     // length vector 
  int methode = Rcpp::as<int>(Rmethode);     // length vector 
  bool intercept = Rcpp::as<bool>(Rintercept);     // length vector 

  //initialisations hors paramètres
  int n = X.rows();
  int p = X.cols(); 
  VectorXi compteur(K);//vecteur des groupes
  compteur <<VectorXi::Zero(K);
 // srand (time(NULL)); //initialisation de la graine aléatoire
  
  for(int i=0;i<n;i++){//recensement des groupes
    compteur[groupe[i]]++;
  }
  int inter=0;
  if(intercept){
    inter++;
  }
  MatrixXd beta(p+inter,1);// déclaration du beta qui ne change jamais de taille
  beta << MatrixXd::Zero(p+inter,1);
  double MSE=0;
  double MSEloc=0;
  int compt_test=0;
  int compt_appr=0;
  for(int j=0;j<K;j++){//pour chaque groupe
    //déclarations locales
    MatrixXd Xappr(n-compteur[j],p);
    MatrixXd Xtest(compteur[j],p+inter);
    MatrixXd Yappr(n-compteur[j],1);
    MatrixXd Ytest(compteur[j],1);
    MatrixXd residus(compteur[j],1);
   
    for(int i=0;i<n;i++){// on remplit les matrices
      if(groupe[i]==j){//echantillon test
        Ytest(compt_test,0)=Y(i,0);
        if(inter){
           Xtest(compt_test,0)=1;
        }
        for(int k=0;k<p;k++){//on remplit toute la ligne de Xtest
          Xtest(compt_test,k+inter)=X(i,k);
        }
        compt_test++;
      }else{// apprentissage
        Yappr(compt_appr,0)=Y(i,0);
        for(int k=0;k<p;k++){//on remplit toute la ligne de Xtest
          Xappr(compt_appr,k)=X(i,k);
        }
        compt_appr++;
      }
    }
    compt_test=0;
    compt_appr=0;
    beta=OLS_cpp(Xappr,Yappr, intercept, methode) ;
    residus=Ytest-Xtest*beta;
    for(int i=0;i<compteur[groupe[j]];i++){//on fait la somme et on passe au carré en même temps
      MSEloc=MSEloc+residus(i,0)*residus(i,0);
    }
    MSEloc=MSEloc/ (double) compteur[groupe[j]];
    MSE=MSE+MSEloc;
    MSEloc=0;//réinitialisation
  }
  MSE=MSE/ (double) K;
  return List::create(
      Named("MSE")=  MSE
    );
  END_RCPP
}
