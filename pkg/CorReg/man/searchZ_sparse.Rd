\name{searchZ_sparse}
\alias{searchZ_sparse}
\title{Sparse structure research}
\usage{
searchZ_sparse(X = X, Zi = NULL, Zj = NULL, Si = NULL, Sj = NULL,
  Bic_null_vect = NULL, methode_tirage = 2, methode_BIC = 1, Rmax = 5,
  Maxiter = 1, plot = F, best = T, better = F, random = T,
  verbose = 1, nb_opt_max = NULL)
}
\arguments{
  \item{X}{the dataset}

  \item{Zi}{indices of the rows of the 1}

  \item{Zj}{indices of the columns of the 1}

  \item{Si}{nombre de 1 dans chaque lignes}

  \item{Sj}{nombre de 1 dans chaque colonnes}

  \item{Bic_null_vect}{the BIC of the null hypothesis (used
  for independent variables)}

  \item{methode_tirage}{0:ligne et colonne,-1:seulement la
  colonne, entier>0:nombre aleatoire de candidats}

  \item{methode_BIC}{1:utilisation de la fonction
  householderQr, 2:utilisation de la fonction
  colPivHouseholderQr}

  \item{p1max}{maximum complexity for a regression}

  \item{Maxiter}{number of steps}

  \item{plot}{TRUE: returns for each step the type of move,
  complexity and BIC}

  \item{best}{TRUE: systematically jumps to the best BIC
  seen ever when seen (it is stored even if best=FALSE)}

  \item{better}{TRUE: systematically jumps to the best
  candidate if better than stationnarity (random wheighted
  jump otherwise)}

  \item{random}{if FALSE:moves only to improve and only to
  the best}

  \item{verbose}{0:none, 1:BIC,step and complexity when
  best BIC found 2:BIC, step, complexity, nb candidates and
  best candidate when best BIC found}

  \item{nb_opt_max}{stop criterion defining how many times
  the chain can walk (or stay) on the max found}
}
\value{
step 0:delete, 1: add, 2: stationnarity
}
\description{
Sparse structure research
}

