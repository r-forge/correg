\name{BicZ}
\alias{BicZ}
\title{Compute the BIC of a given structure}
\usage{
  BicZ(X = X, Z = Z, Bic_null_vect = NULL, BicOld = NULL,
    methode = 1, Zold = NULL)
}
\arguments{
  \item{X}{the dataset}

  \item{Z}{the structure}

  \item{Bic_null_vect}{the BIC of the null hypothesis (used
  for independent variables)}

  \item{BicOld}{BIC (vector) associated to Zold}

  \item{Zold}{another structure with some common parts with
  Z (allows to compute only the differences)}

  \item{methode}{parameter for OLS}
}
\description{
  Compute the BIC of a given structure
}

