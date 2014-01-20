\name{MakeJ}
\alias{MakeJ}
\title{Constructs the jacobian for the constrained likelihood}
\usage{
MakeJ(X = X, Z = Z, B = B, Sigma = Sigma, A = A)
}
\arguments{
  \item{B}{matrice p+1 x p}

  \item{X}{dataset n x p}

  \item{A}{Fixed p+1 coefficient vector for the main
  regression}

  \item{Z}{p squared binary matrix describing the structure
  (adjacency matrix)}

  \item{Sigma}{p2 vector of standard deviations for each
  subregression}
}
\description{
Constructs the jacobian for the constrained likelihood
}

