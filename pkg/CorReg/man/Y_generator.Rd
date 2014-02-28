\name{Y_generator}
\alias{Y_generator}
\title{Response variable generator with a linear model}
\usage{
Y_generator(X = X, Amax = NULL, sigma_Y = 10, positive = 0.6,
  Z = NULL, tp1 = 1, tp2 = 2, tp3 = 3)
}
\arguments{
  \item{X}{the dataset}

  \item{Amax}{maximum number of non-zero coefficients}

  \item{sigma_Y}{the standard deviation of the noise}

  \item{positive}{the ratio of positive coefficients}
}
\description{
Response variable generator with a linear model
}

