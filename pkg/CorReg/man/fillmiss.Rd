\name{fillmiss}
\alias{fillmiss}
\title{Fill the missing values in the dataset}
\usage{
  fillmiss(X = X, Z = NULL, mixmod = F, B = NULL,
    Bt = NULL, res_mixmod = NULL, nbclustmax = 10, X1 = T)
}
\arguments{
  \item{X}{the dataset (matrix) with missing values}

  \item{Z}{the structure associated to X. Can be a matrix
  of zeros if non structure.}

  \item{res_mixmod}{the best results found by mixmod if
  already computed}

  \item{mixmod}{boolean to say if the function has to use
  mixmod hypothesis or just the observed mean.}

  \item{nbclustmax}{the max number of cluster for mixmod.}

  \item{X1}{booleant to say if dependent variables on the
  right will be filled based on the structure}

  \item{Bt}{the matrix used for X1 if X1=TRUE}
}
\description{
  Fill the missing values in the dataset
}

