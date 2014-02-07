\name{R2Z}
\alias{R2Z}
\title{Estimates R2 of each subression}
\usage{
R2Z(Z = Z, X = X, methode = 1, adj = F, crit = c("R2", "F", "sigmaX"))
}
\arguments{
  \item{Z}{binary adjacency matrix of the structure (size
  p)}

  \item{X}{the dataset}

  \item{methode}{parameter for OLS (matrix inversion)
  methode_BIC parameter for OLS (matrix inversion)
  1:householderQr, 2:colPivHouseholderQr}

  \item{adj}{boolean to choose between adjusted R-squared
  and classical one}

  \item{crit}{to choose between the R-squared and the F
  statistic (p-value)}
}
\description{
Estimates R2 of each subression
}

