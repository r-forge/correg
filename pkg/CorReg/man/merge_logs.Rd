\name{merge_logs}
\alias{merge_logs}
\title{merge Rstudio Cran log files in a dataset}
\usage{
merge_logs(logdirectory = "CRANlogs", packages = "CorReg",
  destfile = NULL)
}
\arguments{
  \item{logdirectory}{the path where files are placed}

  \item{packages}{vector of the names of the packages to
  keep in the global file(NULL will keep every packages)}

  \item{destfile}{the name of the created merged file}
}
\description{
merge Rstudio Cran log files in a dataset
}

