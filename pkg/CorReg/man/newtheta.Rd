\name{newtheta}
\alias{newtheta}
\title{newtheta for Newton}
\usage{
newtheta(X = X, Z = Z, B = B, Sigma = Sigma, A = A, lambda = NULL,
  Atilde = Atilde, nbit = 1)
}
\arguments{
  \item{X}{dataset without intercept}

  \item{Z}{the structure (square without intercept)}

  \item{B}{the wheighted structure (p+1)xp}

  \item{Atilde}{p-sized vector of explicative coeffcients}

  \item{Sigma}{vector of the standard deviations of the
  subregressions (p2 long)}

  \item{A}{vector of the coefficient for the main
  regression}

  \item{lambda}{lagrange multiplicators}

  \item{nbit}{number of iteration for Newton}
}
\description{
newtheta for Newton
}

