% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothing.R
\name{ksmooth_bi}
\alias{ksmooth_bi}
\title{Bi-variate kernel smoothing for functional data on common design}
\usage{
ksmooth_bi(Y_list, bw_vec, tobs, tout)
}
\arguments{
\item{Y_list}{List, containing the following elements:
\itemize{
\item \strong{$t} Vector of sampling points.
\item \strong{$X} Matrix of observed points, measured on the bi-dimensional grid containing
cartesian product of \verb{$t} with itself.
}}

\item{bw_vec}{Vector of bandwidths, corresponding to the diagonal elements of
the bandwidth matrix.}

\item{xout}{Matrix of evaluation points in two dimensions.}
}
\value{
List, corresponding to the smooth curves.
}
\description{
Multiplicative kernels are used, with a diagonal bandwidth matrix.
Warning! If smoothing is desired after a change of basis from considering
the directional regularity, remember to also change the evaluation points
in addition to the input list of sampling and observed points.
}
