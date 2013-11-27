#include "cleancolZ.h"
#include "BicZ_cpp2.h"
#include "BicZ_cpp.h"
#include <iostream>
#include <string>
#include <math.h>
using namespace Rcpp ;
using namespace Eigen;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

SEXP cleancolZ(SEXP X,SEXP Z,SEXP bic_vide_vect,SEXP methode_BIC,SEXP plot,SEXP bla)
{
 
  BEGIN_RCPP
  //déclaration des varibles
  const Map<MatrixXd> matZ(as<Map<MatrixXd> >(Z));//Z
  const Map<MatrixXd> matX(as<Map<MatrixXd> >(X));//X
  const Map<VectorXd> Bic_vide_vect(as<Map<VectorXd> >(bic_vide_vect));//bic_vide_vect
  
  Rcpp::NumericVector met_BIC(methode_BIC),Plot(plot),Bla(bla);
  typedef Rcpp::NumericVector::iterator vec_iterator;
  vec_iterator imet_BIC=met_BIC.begin(),iplot=Plot.begin(),ibla=Bla.begin();
  int p=matZ.cols();//nombre de colonne de la matrice Z
  Eigen::MatrixXd Zopt;//meilleur Z obtenu
  double Bicbest;//BIC associé au meilleur modèle
  Eigen::VectorXd bicvect;//vecteur BIC des matrices Z retenues
  Eigen::MatrixXd newZ;//matrice Z modifié à chaque étapes
  Eigen::VectorXd list_cand;//vecteur qui contient les indices des colonnes candidates (liste candidats)
  int compte=0;//permet de cree la matrice liste en désignant le numéro du candidat (liste candidats)
  int nbcand;//nombre de candidats
  int numcand;//numero du candidat
  int i_loc;//premiere coordonnée du candidat (modification Z)
  Eigen::MatrixXd Zcand;//matrice Z du candidat (modification Z)
  Eigen::ArrayXXd SumCol(1,p);//k est un vecteur qui contient la somme de chaque colonne de Zcand (modification Z)
  Eigen::VectorXd BIC_cand;//vecteur qui contient les BIC de chaque colonne des la matrice Zcand (calcul du BIC)
  double Sum_BIC_cand;//somme des BIC de BIC_cand (calcul BIC)
  Eigen::VectorXd stock_BIC;//vecteur qui contient le BIC de tout les candidats (stock)
  double sumbic;//BIC de chaque etapes
  double BIC_min;
  bool minimum;
  int i;
  int iret;
  Eigen::VectorXd bic_etape;//vecteur qui stock le BIC de chaque etapes
  Eigen::VectorXd complexite_etape;//vecteur qui stock la compléxité de Z a chaque etapes
  Eigen::VectorXd BIC_step;//vecteur qui stock le BIC de chaque etapes
  Eigen::VectorXd complexite_step ;//vecteur qui stock la compléxité de Z a chaque etapes
  //initialisation
  bicvect=BicZ_cpp2(matX,matZ,Bic_vide_vect,imet_BIC[0]);
//somme a la main
  sumbic=bicvect.sum();
  if (ibla[0]>0)
  {
    Rcout<<sumbic<<"\n";
  }
  Bicbest =sumbic;
  Zopt=matZ;
  newZ=matZ;
  
  int step =0;
  SumCol=matZ.colwise().sum();//nombre d'éléments dans chaque colonne
  nbcand=(SumCol>0).count();//nombre de colonnes candidates (quand il y'a au moins 1 élément dans une colonne)
  list_cand.resize(nbcand);//le nombre de candidats est la compléxité du modèle
  bic_etape.resize(nbcand);
  complexite_etape.resize(nbcand);
  while(step<nbcand)  
  {
    compte=0;//initialisation du vecteur liste désigne le numéro du candidat
    
  //liste candidats (couples [i,j])
    for(int n=0;n<p;n++)//parcours les colonnes
    {
      if(newZ.col(n).sum()>0)//on cherche les candidats
      {
        list_cand(compte)=n;//stock la colonne du candidat
        compte=compte+1;//passage au candidat suivant
      }
    }
    
    stock_BIC.resize(compte);
    //pour chaque candidat
    for (numcand=0;numcand<compte;numcand++)
    {
    //modification (calcul du Z)    (Zcand avec methode (rejet ou relax)   Z,i,j,methode="relax", p2max=inf,rmax=inf)
      
      Zcand=newZ;
      for(i_loc=0;i_loc<p;i_loc++)//parcours des lignes
      {
        Zcand(i_loc,list_cand(numcand))=0;//modification de Z par rapport a la colonne candidate
      }
      
    //calcul du bic (du nouveau Z généré)
      BIC_cand=BicZ_cpp(matX,Zcand,Bic_vide_vect,bicvect,imet_BIC[0],newZ);//calcul du vecteur BIC du candidat
      Sum_BIC_cand=BIC_cand.sum();
    //stockage du BIC
      
      stock_BIC(numcand)=Sum_BIC_cand;
    }
    //on regarde la meilleur candidat
    BIC_min=Bicbest;//initialisation BIC_min
    minimum=false;
    for(i=0;i<compte;i++)
    {
      if(stock_BIC(i)<Bicbest)//si on a trouvé un minimum
      {
        iret=i;//on retient sa position
        minimum=true;
        Bicbest=stock_BIC(i);//on modifie Bicbest
      }
    }
    if (minimum==true)//si il y a un minimum, on y va
    {
      Zopt=newZ;//initialisation à la stationarité
      for(int i_opt=0;i_opt<p;i_opt++)//on met la colonne à 0
      {
        Zopt(i_opt,list_cand(iret))=0;
      }
      bicvect=BicZ_cpp(matX,Zopt,Bic_vide_vect,bicvect,imet_BIC[0],newZ);//calcul du vecteur BIC du nouveau Zopt

      newZ=Zopt;
      
      //partie commentaires
      if (ibla[0]>0)
      {
        Rcout<<"etape : "<<step<<"\n";
        Rcout<<" Bicbest: "<<Bicbest<<" complexite: "<<newZ.sum()<<"\n";
        if (ibla[0]==2)
        {
          Rcout<<"nb_cand "<<compte<<"BIC min "<<BIC_min<<"\n";
        }
      }
      
      //partie graphiques
      if(iplot[0]==1)//si on veut des graphiques on mets à jour les vecteurs
      {
        bic_etape(step)=Bicbest;
        complexite_etape(step)=newZ.sum();
      }
      
    step=step+1;//passage a l'etape suivante s'il y a eu un minimum
    }
    else//pas de minimum
    {
      if (ibla[0]>0){
        Rcout<<"nettoyage fini";
      }      
      step=nbcand+1;//arret des étapes
    }
  }//fin des étapes
  
  if(iplot[0]==1)//si on veut des graphiques on mets à jour les vecteurs
  {
    BIC_step.resize(nbcand-compte);//on redimensionne pour les plot
    complexite_step.resize(nbcand-compte);
    if(nbcand-compte>0)//si il y'a eu au moins une étape
    {
      for(int k=0;k<nbcand-compte;k++)
      {
        BIC_step(k)=bic_etape(k);
        complexite_step(k)=complexite_etape(k);
      }
    }
    return List::create(
      Named("Z")=  newZ,
      Named("Zopt")=  Zopt,
      Named("bic_opt")=  Bicbest,
      Named("bic_etape")=  BIC_step,
      Named("complexite_etape")=  complexite_step
    );
  }
  else
  {
  return List::create(
      Named("Z")=  newZ,
      Named("Zopt")=  Zopt,
      Named("bic_opt")=  Bicbest
    );
  }
  END_RCPP
}  
