\name{calcul_BIC_mixmod}
\alias{calcul_BIC_mixmod}
\title{calcul du BIC initial via package externe Mixmod (ou mclust)}
\usage{
  calcul_BIC_mixmod(X = X, nbclustmax = 10, bla = FALSE,
    details = FALSE, max = TRUE, mclust = TRUE)
}
\arguments{
  \item{X}{the dataset (matrix)}

  \item{nbclsutmax}{max number of clusters in the gaussian
  mixtures}

  \item{bla}{verbose or not}

  \item{details}{boolean to give the details of the
  mixtures found}

  \item{max}{boolean. Use an heuristic to shrink nbclustmax
  according to the number of individuals in the dataset}

  \item{mclust}{boolean. Use mclust instead of Rmixmod}
}
\description{
  calcul du BIC initial via package externe Mixmod (ou
  mclust)
}

