\name{cleancolZ}
\alias{cleancolZ}
\title{clean Z columns (if BIC improved)}
\usage{
cleancolZ(X = X, Z = Z, bic_vide_vect = bic_vide_vect, methode_BIC = 1,
  plot = F, bla = 1, star = FALSE)
}
\arguments{
  \item{X}{the dataset}

  \item{Z}{matrice Z a nettoyer}

  \item{bic_vide_vect}{vecteur BIC de la matrice nulle}

  \item{methode_BIC}{1:utilisation de la fonction
  householderQr, 2:utilisation de la fonction
  colPivHouseholderQr}

  \item{star}{to use BIC*}

  \item{plot}{T:retourne le BIC de chaque etapes}

  \item{bla}{0:pas de messages, 1:affiche le BIC,le numero
  d'etape et la complexite de Z quand il y'a un meilleur
  BIC, 2:affiche le BIC,le numero d'etape,la complexite de
  Z,le nombre de candidats et le BIC minimum observe parmi
  les candidats quand il y'a un meilleur BIC}
}
\value{
etape 0:suppression,1 ajout,2 stationarite
}
\description{
clean Z columns (if BIC improved)
}

