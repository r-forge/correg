\name{generateurA_ou}
\alias{generateurA_ou}
\title{generateur de parametres selon la structure}
\usage{
generateurA_ou(Z = Z, tp1 = 1, tp2 = 1, tp3 = 1, positif = 0.6,
  lambdapois = 5, pb = 2, Amax = NULL, B = NULL)
}
\arguments{
  \item{pb}{kind of problem : 0=none, 1=simple, 2=strong}

  \item{Z}{the structure adjacency matrix}

  \item{tp1}{ratio of the right-side covariates in A}

  \item{tp2}{ratio of the left-side covariates in A}

  \item{tp3}{ratio of the totally independent covariates
  (if exists) in A}

  \item{positif}{ratio of positive coefficients in A}

  \item{lambdapois}{parameter of the poisson law to
  generate the coefficients}

  \item{Amax}{max number of non-zero coefficients}

  \item{B}{weighted structure used to coerce some problems
  in the complete model}
}
\description{
generateur de parametres selon la structure
}

