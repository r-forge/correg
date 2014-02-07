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

  \item{reject}{0: constraint relaxation, 1: reject mode}

  \item{p1max}{est le nombre max de 1 sur une colonne de Z}

  \item{random}{indique si on fait le tirage au sort ou pas
  (si on utilise BIc on a moins besoin du tirage). sans
  tirage = sans tenir compte de W}

  \item{mode}{indique la strategie pour le choix des
  candidats a la mise a 1}

  \item{p2max}{est le nombre max de sous-reg}

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

