\name{calcul_BIC_mixmod}
\alias{calcul_BIC_mixmod}
\title{BIC of the marginal densities}
\usage{
calcul_BIC_mixmod(X = X, nbclustmax = 10, bla = FALSE, details = FALSE,
  max = TRUE, mclust = TRUE, nbini = 50)
}
\arguments{
  \item{X}{the dataset (matrix)}

  \item{nbclustmax}{max number of clusters in the gaussian
  mixtures}

  \item{bla}{verbose or not}

  \item{details}{boolean to give the details of the
  mixtures found}

  \item{max}{boolean. Use an heuristic to shrink nbclustmax
  according to the number of individuals in the dataset}

  \item{mclust}{boolean. Use mclust instead of Rmixmod}

  \item{nbini}{number of initial points for Rmixmod}
}
\description{
Estimates the density of each covariates with gaussian
mixture models and then gives the associated BIC.
}

