\name{Winitial}
\alias{Winitial}
\title{initialisation based on a wheight matrix (correlation or other)}
\usage{
  Winitial(W = W, X = NULL, Rmax = NULL, BIC = F,
    Bic_null_vect = NULL, relax = T, random = T,
    nbclustmax = 5, sorted = T,
    mode = c("all", "sorted", "multinom"), p2max = NULL)
}
\arguments{
  \item{W}{est la matrice des poids (entre 0 et 1)}

  \item{X}{est la matrice des donnees (ne sert que si
  BIC=T)}

  \item{Bic_null_vect}{est le vecteur bic vide (issu de
  mixmod par exemple)}

  \item{relax}{indique si on fait la relaxation ou le rejet
  (booleen)}

  \item{Rmax}{est le nombre max de 1 sur une colonne de Z}

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
  initialisation based on a wheight matrix (correlation or
  other)
}

