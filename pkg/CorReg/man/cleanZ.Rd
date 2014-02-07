\name{cleanZ}
\alias{cleanZ}
\title{clean Z (if BIC improved)}
\usage{
cleanZ(X = X, Z = Z, Bic_null_vect = Bic_null_vect, methode = 1,
  plot = F, verbose = 1, star = FALSE)
}
\arguments{
  \item{X}{the dataset}

  \item{Z}{binary adjacency matrix of the structure (size
  p)}

  \item{Bic_null_vect}{the BIC of the null hypothesis (used
  for independent variables)}

  \item{methode}{parameter for OLS (matrix inversion)
  methode_BIC parameter for OLS (matrix inversion)
  1:householderQr, 2:colPivHouseholderQr}

  \item{star}{boolean defining wether classical BIC or BIC*
  is computed}

  \item{plot}{if TRUE returns the vector of BIC for each
  step}

  \item{verbose}{0:none, 1:BIC,step and complexity when
  best BIC found 2:BIC, step, complexity, nb candidates and
  best candidate when best BIC found}
}
\description{
clean Z (if BIC improved)
}

