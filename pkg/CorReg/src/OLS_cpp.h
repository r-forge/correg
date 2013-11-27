#ifndef _CorReg_OLS_cpp_H
#define _CorReg_OLS_cpp_H

#include <RcppEigen.h>
using namespace Rcpp ;
using namespace Eigen;

RcppExport MatrixXd OLS_cpp(MatrixXd X, MatrixXd Y , bool intercept,int methode) ;

#endif
