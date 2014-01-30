\name{BicZcurve}
\alias{BicZcurve}
\title{Curve of the BIC for each possible p2 with a fixed Z and truncature of Z}
\usage{
BicZcurve(X = X, Z = Z, Bic_null_vect = Bic_null_vect, plot = T,
  star = F, trunc = NULL)
}
\arguments{
  \item{X}{matrix containing the dataset}

  \item{Z}{adjacency matrix (binary) describing the
  structure between the variables}

  \item{Bic_null_vect}{vector of the BIC for each variable.
  used when the variable is independant}

  \item{plot}{boolean to plot or not the curve}

  \item{star}{boolean to use BIC* (hierarchical uniform law
  on the structure)}

  \item{trunc}{number of sub-regression to keep (best R2).
  if NULL the min of BIC is kept}
}
\description{
Curve of the BIC for each possible p2 with a fixed Z and
truncature of Z
}

