% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/w.lcor.R
\name{w.lcor}
\alias{w.lcor}
\title{Calculate the local wavelet correlation.}
\usage{
w.lcor(data, w, scales)
}
\arguments{
\item{data}{the dataset of species and locations.}

\item{w}{the result of species dataset calculate by 'specwt'.}

\item{scales}{the scale selected to analysis.}
}
\value{

}
\description{
Calculate the local wavelet correlation.
}
\examples{
w.lcor(data = Fauna.data[,-12], w = specwt(1:119, Fauna.data[,-12]), scales = 57)
}
