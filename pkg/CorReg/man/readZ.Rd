\name{readZ}
\alias{readZ}
\title{read the structure and explain it}
\usage{
  readZ(Z = Z, B = NULL, R2 = T, varnames = NULL,
    output = c("index", "names", "all"), X = NULL,
    decr = TRUE)
}
\arguments{
  \item{Z}{is the binary structure}

  \item{B}{is the complete structure (if known)}

  \item{R2}{boolean. Indicates wether you want to order Z
  by R2 (adjusted)}

  \item{varnames}{the names of the variables (same order)}

  \item{output}{indicates the content of the output}

  \item{X}{is a dataframe containing the dataset}

  \item{decr}{boolean. Define the order used (don't have to
  scroll to see the best subregresisons if many)}
}
\description{
  read the structure and explain it
}

