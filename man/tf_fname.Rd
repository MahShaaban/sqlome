% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tf_funs.R
\name{tf_fname}
\alias{tf_fname}
\title{Construct a file name for a Cistrome Cancer file}
\usage{
tf_fname(tf, dir, all = TRUE)
}
\arguments{
\item{tf}{A character string of the official symbol of the transcription factor}

\item{dir}{A character string for a directory to download files.}

\item{all}{A logicall of whether to get the correlation files of all genes (default) or only targets.}
}
\value{
A character string of the cistrome cancer file name for a transcription factor
}
\description{
Construct a file name for a Cistrome Cancer file
}
\examples{
tf_fname('AFF4', dir = 'tmp/tf/')

}
