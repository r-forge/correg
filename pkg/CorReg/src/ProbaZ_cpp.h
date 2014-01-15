#ifndef _CorReg_BicLoc_cpp_H
#define _CorReg_BicLoc_cpp_H

#include <RcppEigen.h>
using namespace Rcpp ;
RcppExport double BicLoc_cpp(Eigen::MatrixXd X,Eigen::MatrixXd Y,bool intercept,int methode) ;

#endif
