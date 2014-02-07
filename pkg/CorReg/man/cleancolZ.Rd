\name{cleancolZ}
\alias{cleancolZ}
\title{clean Z columns (if BIC improved)}
\usage{
cleancolZ(X = X, Z = Z, Bic_null_vect = bic_vide_vect, methode_BIC = 1,
  plot = F, verbose = 1, star = FALSE)
}
\arguments{
  \item{X}{the dataset}

  \item{Z}{binary adjacency matrix of the structure (size
  p)}

  \item{Bic_null_vect}{vector of the BIC for each
  covariate}

  \item{methode_BIC}{parameter for OLS (matrix inversion)
  1:householderQr, 2:colPivHouseholderQr}

  \item{star}{to use BIC*}

  \item{plot}{if TRUE returns the vector of BIC for each
  step}

  \item{verbose}{0:none, 1:BIC,step and complexity when
  best BIC found 2:BIC, step, complexity, nb candidates and
  best candidate when best BIC found}
}
\description{
clean Z columns (if BIC improved)
}

