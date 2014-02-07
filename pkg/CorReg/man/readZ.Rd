\name{readZ}
\alias{readZ}
\title{read the structure and explain it}
\usage{
readZ(Z = Z, B = NULL, crit = c("none", "R2", "F", "sigmaX"),
  varnames = NULL, output = c("index", "names", "all"), X = NULL,
  order = 1)
}
\arguments{
  \item{Z}{binary adjacency matrix of the structure (size
  p)}

  \item{B}{is the complete structure (wheighted)}

  \item{crit}{define the criterion to use
  c("none","R2","F","sigmaX")}

  \item{varnames}{the names of the variables (same order)}

  \item{output}{indicates the content of the output
  output=c("index","names","all")}

  \item{X}{is a dataframe or matrix containing the dataset}

  \item{order}{Define the order used (0: none, -1:
  decreasing, 1: growing)}
}
\description{
read the structure and explain it
}

