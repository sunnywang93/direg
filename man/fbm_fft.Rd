% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{fbm_fft}
\alias{fbm_fft}
\title{Simulate 1D fractional brownian motion}
\usage{
fbm_fft(H, n, grid_max)
}
\arguments{
\item{H}{Numeric, the Hurst index of the fractional
brownian motion.}

\item{n}{Numeric, the number of evaluation points.}

\item{grid_max}{Numeric, the right-end point of the sampling interval.
The observed points are simulated on the interval \verb{[0, grid_max]}, using
the self-similar property of fbms.}
}
\value{
Vector, the observed points of the fbm.
}
\description{
Performs simulation of one-dimensional fractional
brownian motion on an equally spaced grid using fast-fourier transforms.
}
