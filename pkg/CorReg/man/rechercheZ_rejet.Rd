\name{rechercheZ_rejet}
\alias{rechercheZ_rejet}
\title{recherche de structure (mode rejet)}
\usage{
  rechercheZ_rejet(X = X, Z = NULL,
    bic_vide_vect = bic_vide_vect, methode_tirage = 0,
    methode_BIC = 1, Rmax = 5, p2max = NULL,
    Maxiter = Maxiter, plot = F, best = T, better = F,
    random = T, bla = 1, nb_opt_max = NULL, Mixmod = T,
    exact = F)
}
\arguments{
  \item{X}{}

  \item{Z}{Z est une matrice nulle si on ne lui a pas mis
  de valeur}

  \item{bic_vide_vect}{vecteur BIC de la matrice nulle}

  \item{methode_tirage}{0:ligne et colonne,-1:seulement la
  colonne,-2:tout les candidats, entier>0:nombre aleatoire
  de candidats}

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

  \item{nb_opt_max}{nombre maximum de fois ou on reste
  stationnaire}

  \item{Mixmod}{}
}
\value{
  etape 0:suppression,1 ajout,2 stationarite
}
\description{
  recherche de structure (mode rejet)
}

