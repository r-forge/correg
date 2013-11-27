\name{generateur_melange}
\alias{generateur_melange}
\title{generateur de donnees en sous-regression par modele de melange gaussien}
\usage{
  generateur_melange(n = 130, p = 100, ratio = 0.4,
    max_compl = 3, valid = 1000, positif = 0.6,
    sigma_Y = 10, sigma_sousreg = 0.25, meanvar = NULL,
    sigmavar = NULL, lambda = 5, Amax = 15, lambdapois = 5,
    gamma = T, gammashape = 1, gammascale = 0.5, tp1 = 1,
    tp2 = 1, tp3 = 1, pb = 2)
}
\arguments{
  \item{lambda}{paramater of the law that define the number
  of components in gaussian mixture models}
}
\description{
  generateur de donnees en sous-regression par modele de
  melange gaussien
}

