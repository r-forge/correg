\name{log_manager}
\alias{log_manager}
\title{managing CRAN logs.}
\usage{
log_manager(start = as.Date("2012-10-01"), today = format(Sys.Date(),
  "\%Y-\%m-\%d"), logdirectory = "CRANlogs", global = FALSE,
  packages = NULL)
}
\arguments{
  \item{start}{starting date}

  \item{today}{last date}

  \item{logdirectory}{the path where files will be placed}

  \item{global}{boolean to define wether a global file will
  be created or not}

  \item{packages}{vector of the names of the packages to
  keep in the global file(NULL will keep every packages)}
}
\description{
code adapted from
http://www.nicebread.de/finally-tracking-cran-packages-downloads/
(original code is under the FreeBSD license, by Felix
Schönbrodt)
}

