\name{confint_coef}
\alias{confint_coef}
\title{plot and give confidence intervals on the coefficients estimated in a model or for proportions}
\usage{
confint_coef(modele = NULL, n = NULL, prop = NULL, mean = NULL,
  alpha = 0.05, labels = NULL)
}
\arguments{
  \item{modele}{a model from lm on whoch to compute the
  confidence intervals}

  \item{n}{a vector of quantities associated to prop}

  \item{prop}{a vector of proportions (between 0 and 1)}

  \item{mean}{a mean value to plot}

  \item{alpha}{the risk (confidence 1-alpha)}

  \item{labels}{a vector of names to put below the bars.}
}
\description{
plot and give confidence intervals on the coefficients
estimated in a model or for proportions
}

