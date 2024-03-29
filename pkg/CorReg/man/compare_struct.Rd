% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_struct.R
\name{compare_struct}
\alias{compare_struct}
\title{To compare sub-regression structures}
\usage{
compare_struct(trueZ = trueZ, Zalgo = Zalgo, all = TRUE,
  mode = "NULL")
}
\arguments{
\item{trueZ}{first structure (binary adjacency matrix)}

\item{Zalgo}{second structure (binary adjacency matrix)}

\item{all}{(boolean) Also compute the ratio for each stat.}

\item{mode}{how to modify the structures before comparison. mode=c("NULL","hybrid","clique","sym")
It allows to compare groups instead of exact sub-regressions. Does nothing by default.}
}
\value{
\item{true1}{ Number of links that exist in both matrices }
\item{false1}{ Number of links that exist only in Zalgo }
\item{false0}{ Number of links that exist only in trueZ }
\item{deltadr}{Number of sub-regressions in trueZ -  Number of sub-regressions in Zalgo (i.e. : negative if too much sub-regressions in Zalgo)}
\item{true_left}{ Number of variables redundant in both matrices }
\item{false_left}{ Number of variables redundant in Zalgo but not in trueZ  }
\item{ratio_true1}{ ratio of links in trueZ that exist also in Zalgo  }
\item{ratio_true0}{ ratio of links not in trueZ that don't exist in  Zalgo.}
}
\description{
Compares two sub-regression structures, considering one of them as the "true one".
}
