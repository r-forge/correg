\name{recursive_tree}
\alias{recursive_tree}
\title{decision tree in a recursive way}
\usage{
recursive_tree(data = data, Y = "Y", modele = NULL, kill = NULL,
  index = NULL, print = TRUE, plot = TRUE)
}
\arguments{
  \item{data}{the dataset including the response}

  \item{Y}{the name of the response}

  \item{modele}{vector of the names used for the last tree
  search (not necessarily choosen)}

  \item{kill}{vector of the names to kill}

  \item{index}{to give a number to the plot}

  \item{print}{boolean to print the tree parameters}

  \item{plot}{boolean to plot the tree}
}
\description{
decision tree in a recursive way
}

