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
  \item{n}{the number of individuals in the learning
  dataset}

  \item{p}{the number of covariates (without the response)}

  \item{ratio}{the ratio of explained covariates
  (dependent)}

  \item{max_compl}{the number of covariates in each
  subregression}

  \item{valid}{the size of the validation sample}

  \item{positif}{the ratio of positive coefficients in both
  the regression and the subregressions}

  \item{sigma_Y}{standard deviation for the noise of the
  regression}

  \item{sigma_sousreg}{standard deviation for the noise of
  the subregression (all). ignored if gamma=T}

  \item{gamma}{boolean to generate a p-sized vector
  sigma_sousreg gamma-distributed}

  \item{gammashape}{shape parameter of the gamma
  distribution (if needed)}

  \item{gammascale}{scale parameter of the gamma
  distribution (if needed)}

  \item{meanvar}{vector of means for the covariates.}

  \item{sigmavar}{standard deviation of the covariates.}

  \item{lambda}{paramater of the law that define the number
  of components in gaussian mixture models}

  \item{Amax}{the maximum number of covariates with
  non-zero coefficients in the regression}

  \item{tp1}{the ratio of right-side covariates allowed to
  have a non-zero coefficient in the regression}

  \item{tp2}{the ratio of left-side covariates allowed to
  have a non-zero coefficient in the regression}

  \item{tp3}{the ratio of strictly independent covariates
  allowed to have a non-zero coefficient in the regression}

  \item{lambdapois}{parameter used to generate the
  coefficient in the subregressions. poisson distribution.}

  \item{pb}{Defines an heuristic to generate Y in a way
  that will give some issues with correlations.}
}
\description{
  generateur de donnees en sous-regression par modele de
  melange gaussien
}

