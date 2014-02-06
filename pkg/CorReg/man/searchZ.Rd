\name{searchZ}
\alias{searchZ}
\title{MCMC algo to find a structure between the covariates}
\usage{
searchZ(X = X, Z = NULL, Bic_null_vect = NULL, methode_tirage = -1,
  reject = 1, methode_BIC = 1, Rmax = 5, p2max = NULL, Maxiter = 1,
  plot = FALSE, best = TRUE, better = FALSE, random = TRUE, bla = 1,
  nb_opt_max = NULL, exact = TRUE, nbini = NULL, star = TRUE,
  nett = TRUE, ...)
}
\arguments{
  \item{X}{the dataset}

  \item{Z}{binary adjacency matrix of size p. if NULL zero
  matrix is used}

  \item{Bic_null_vect}{vector of the BIC for each
  covariate}

  \item{methode_tirage}{0:row and column,-1:column only,
  int>0:random int candidates, -2 : all (but the diag), -3
  : non-zeros}

  \item{reject}{0:mode relax, 1:mode reject}

  \item{methode_BIC}{1:utilisation de la fonction
  householderQr, 2:utilisation de la fonction
  colPivHouseholderQr}

  \item{Rmax}{complexite maximum d'une sous-regression}

  \item{p2max}{nombre maximal de sous-regressions}

  \item{Maxiter}{nombre d'etapes}

  \item{plot}{T:retourne le type de changement, la
  complexite et le BIC de chaque etapes}

  \item{best}{T:permet d'aller systematiquement au meilleur
  BIC si il est meilleur que tout les autres deja
  rencontres}

  \item{better}{T:permet d'aller systematiquement au
  meilleur BIC si il est meilleur que l'etape precedente}

  \item{random}{F:permet de s'ameliorer ou de rester sur
  place}

  \item{bla}{0:pas de messages, 1:affiche le BIC,le numero
  d'etape et la complexite de Z quand il y'a un meilleur
  BIC, 2:affiche le BIC,le numero d'etape,la complexite de
  Z,le nombre de candidats et le BIC minimum observe parmi
  les candidats quand il y'a un meilleur BIC, 3: affiche en
  plus de bla=1 la complexite locale et le BIC local}

  \item{nb_opt_max}{stop criterion defining how many times
  the chain can stay at the max found}

  \item{exact}{boolean. If exact subregression is found it
  gives its content.}

  \item{nbini}{Number of initialisations (using Winitial).
  if NULL and Zini is NULL : only one chain beginning with
  zero matrix.}

  \item{star}{boolean to compute BIC* instead of BIC
  (stronger penalization of the complexity)}

  \item{nett}{cleaning steps at the ends}

  \item{...}{parameters to be passed (for Winitial).}
}
\value{
etape 0:suppression,1 ajout,2 stationarite
}
\description{
MCMC algo to find a structure between the covariates
}

