\name{readZ}
\alias{readZ}
\title{read the structure and explain it}
\usage{
readZ(Z = Z, B = NULL, R2 = TRUE, varnames = NULL, output = c("index",
  "names", "all"), X = NULL, order = 1)
}
\arguments{
  \item{Z}{is the binary structure}

  \item{B}{is the complete structure (if known)}

  \item{R2}{boolean. Indicates wether you want to order Z
  by R2 (adjusted)}

  \item{varnames}{the names of the variables (same order)}

  \item{output}{indicates the content of the output}

  \item{X}{is a dataframe containing the dataset}

  \item{order}{Define the order used (0: none, -1:
  decreasing, 1: growing)}
}
\description{
read the structure and explain it
}

