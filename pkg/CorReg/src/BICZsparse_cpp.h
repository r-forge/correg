#ifndef _CorReg_BICZsparse_cpp_H
#define _CorReg_BICZsparse_cpp_H

#include <RcppEigen.h>
using namespace Rcpp ;
RcppExport Eigen::VectorXd BICZsparse_cpp(Eigen::MatrixXd X,NumericVector Z_Zi,Eigen::VectorXd Z_Sj,Eigen::VectorXd Bic_vide_vect,Eigen::VectorXd BicOld,int methode,Eigen::VectorXd Zold_Sj,int nrow) ;

#endif
