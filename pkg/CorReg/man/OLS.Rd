\name{OLS}
\alias{OLS}
\title{Ordinary Least Square efficiently computed with SEM for missing values}
\usage{
OLS(X = X, Y = Y, M = NULL, intercept = F, sigma = F, Bic = F,
  methode = 1, miss = 0, mixmod = NULL, nbit = 100, nbclustmax = 10)
}
\arguments{
  \item{X}{the covariates (double)}

  \item{Y}{the response}

  \item{intercept}{(boolean) is an intercept intended ?}

  \item{sigma}{(boolean) is it necessary to compute the
  standard deviation of errors ?}

  \item{Bic}{(boolean) is the BIC criterion computation
  intended ?}

  \item{methode}{how to compute the least squares ?}

  \item{miss}{to indicate wether there are missing values
  in X (miss=1) or not (miss=0) next version would allow
  missing values in Y}

  \item{nbit}{number of iteration for SEM}

  \item{nbclustmax}{max number of cluster for mixmod
  (ignored if miss=0)}

  \item{M}{binary matrix (size of X) with 1 where X is
  missing}

  \item{mixmod}{Gaussian Mixture hypothesis if needed. Or
  result of
  calcul_BIC_mixmod(X=X,nbclustmax=nbclustmax,bla=F,details=T)}
}
\description{
Ordinary Least Square efficiently computed with SEM for
missing values
}

