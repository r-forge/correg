\name{GM_Loglikelihood}
\alias{GM_Loglikelihood}
\title{Compute the log-likelihood of a gaussian mixture value due to missing values in a regression}
\usage{
  GM_Loglikelihood(Y = Y, X = X, B = B, sigma = sigma,
    M = NULL, mixmod, log = T, intercept = T)
}
\arguments{
  \item{X}{is the vector of the regressors}

  \item{Y}{is the value of the response variable (scalar)}

  \item{B}{is the vector of regression parameters,
  including the intercept in first position (zero if no
  intercept)}

  \item{sigma}{is the standard deviation of the residual of
  the subregression}

  \item{M}{is the indicatrice vector (submatrix) of the
  missing values}

  \item{mixmod}{is result of
  calcul_BIC_mixmod2.0(X=dataset,nbclustmax=nbclustmax,details=T)}

  \item{log}{boolean to define if you want the likelihood
  or the log-likelihood}

  \item{intercept}{boolean to define if there is an
  intercept in the model (default=T)}
}
\description{
  Compute the log-likelihood of a gaussian mixture value
  due to missing values in a regression
}

