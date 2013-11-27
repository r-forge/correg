\name{rechercheZ_sparse_relax}
\alias{rechercheZ_sparse_relax}
\title{recherche de structure au format creux et en mode relax}
\usage{
  rechercheZ_sparse_relax(X = X, Zi = Zi, Zj = Zj, Si = Si,
    Sj = Sj, bic_vide_vect = bic_vide_vect,
    methode_tirage = 0, methode_BIC = 1, Rmax = 5,
    Maxiter = Maxiter, plot = F, best = T, better = F,
    random = T, bla = 1, nb_opt_max = NULL)
}
\arguments{
  \item{X}{}

  \item{Zi}{numeros des lignes de chaque 1}

  \item{Zj}{numeros des colonnes de chaque 1}

  \item{Sj}{nombre de 1 dans chaque lignes}

  \item{Sj}{nombre de 1 dans chaque colonnes}

  \item{bic_vide_vect}{vecteur BIC de la matrice nulle}

  \item{methode_tirage}{0:ligne et colonne,-1:seulement la
  colonne, entier>0:nombre aleatoire de candidats}

  \item{methode_BIC}{1:utilisation de la fonction
  householderQr, 2:utilisation de la fonction
  colPivHouseholderQr}

  \item{Rmax}{complexite maximum d'une sous-regression}

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

  \item{nb_opt_max}{}

  \item{Mixmod}{}
}
\value{
  etape 0:suppression,1 ajout,2 stationarite
}
\description{
  recherche de structure au format creux et en mode relax
}

