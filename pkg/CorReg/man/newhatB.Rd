\name{newhatB}
\alias{newhatB}
\title{Computes B based on A and Atilde under the constraint of the subregressions}
\usage{
newhatB(Z = Z, A = A, Atilde = Atilde, Bold = B, intercept = TRUE)
}
\arguments{
  \item{Z}{the structure p x p binary matrix}

  \item{A}{the vector of coefficients}

  \item{Atilde}{the explicative model}

  \item{Bold}{the old B}

  \item{intercept}{boolean to define an intercept}
}
\description{
Computes B based on A and Atilde under the constraint of
the subregressions
}

