\name{Winitial}
\alias{Winitial}
\title{initialization based on a wheight matrix (correlation or other)}
\usage{
Winitial(W = W, X = NULL, p1max = NULL, BIC = F, Bic_null_vect = NULL,
  relax = T, random = T, nbclustmax = 5, sorted = T, mode = c("all",
  "sorted", "multinom"), p2max = NULL)
}
\arguments{
  \item{W}{Weight matrix with values between 0 and 1 ( e.g.
  : abs(cor(X)) )}

  \item{X}{the dataset (if BIC=TRUE)}

  \item{Bic_null_vect}{the BIC of the null hypothesis (used
  for independent variables)}

  \item{relax}{TRUE: constraint relaxation, FALSE: reject
  mode}

  \item{p1max}{maximum complexity for a regression}

  \item{random}{if FALSE W is not used.}

  \item{mode}{candidates strategy}

  \item{p2max}{maximum number of regressions}

  \item{sorted}{boolean to sort the candidates based on W
  (best first)}

  \item{BIC}{boolean to use BIC for probabilities.}

  \item{nbclustmax}{parameter for calcul_BIC_mixmod if
  needed}
}
\description{
initialization based on a wheight matrix (correlation or
other)
}

