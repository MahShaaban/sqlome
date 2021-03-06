% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tf_funs.R
\name{tf_url}
\alias{tf_url}
\title{Construct a url for Cistrome Cancer file}
\usage{
tf_url(tf, all = TRUE)
}
\arguments{
\item{tf}{A character string of the official symbol of the transcription factor}

\item{all}{A logicall of whether to get the correlation files of all genes (default) or only targets.}
}
\value{
A character string of the cistrome cancer url for a transcription factor
url <- tf_url('AFF4')
}
\description{
This function constructs a url to downlaod the tarnscription factors correlation data from Cistrome Cancer. Briefly, the function constructs a url
for the file using a predefined handle and a suffix.
}
