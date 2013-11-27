\name{BICZsparse}
\alias{BICZsparse}
\title{calcul du BIC en format creux}
\usage{
  BICZsparse(X = X, Z_Zi = Z_Zi, Z_Sj = Z_Sj,
    Bic_vide_vect = Bic_vide_vect, BicOld = BicOld,
    methode = methode, Zold_Sj = Zold_Sj, nrow = nrow)
}
\arguments{
  \item{Z_Zi}{numero des lignes des 1 dans la matrice Z}

  \item{vecteur}{du nombre de 1 dans chaque colonne}

  \item{Bic_vide_vect}{BIC vectoriel de la matrice nulle}

  \item{BicOld}{BIC vectoriel de la matrice precedente}

  \item{methode}{1:utilisation de la fonction
  householderQr, 2:utilisation de la fonction
  colPivHouseholderQr}

  \item{Zold_Sj}{nombre de 1 dans la colonne de la matrice
  precedente}

  \item{nrow}{nombre de ligne de la matrice x}
}
\description{
  calcul du BIC en format creux
}

