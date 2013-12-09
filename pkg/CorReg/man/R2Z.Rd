\name{R2Z}
\alias{R2Z}
\title{Estimates R2 of each subression}
\usage{
R2Z(Z = Z, X = X, methode = 1, adj = F, crit = c("R2", "F"))
}
\arguments{
  \item{Z}{is the structure (binary matrix)}

  \item{X}{the dataset}

  \item{methode}{parameter for Ordinary Least square (don't
  change it)}

  \item{adj}{boolean to choose between adjusted R-squared
  and classical one}

  \item{crit}{to choose between the R-squared and the F
  statistic (p-value)}
}
\description{
Estimates R2 of each subression
}

