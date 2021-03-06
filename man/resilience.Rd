% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resilience.R
\name{resilience}
\alias{resilience}
\title{Resilience index}
\usage{
resilience(
  d,
  sites,
  surveys,
  Rt.quantile = NULL,
  Rt.value = NULL,
  resistance = FALSE
)
}
\arguments{
\item{d}{a symmetric \code{matrix} or an object of class \code{dist} cointaining
the dissimilarity values between pairs of community states (pre-disturbance,
disturbed and post-disturbance states). \code{d} values must be in the range [0, 1].}

\item{sites}{A vector indicating the site of each community state.}

\item{surveys}{A vector indicating the survey of each community state.}

\item{Rt.quantile}{A numeric value of probability in [0, 1] to take the
corresponding \code{\link{quantile}} as threshold of the resistance index below
which the resilience index is calculated.}

\item{Rt.value}{A numeric value in [0, 1] representing the threshold of the
resistance index below which the resilience index is calculated.}

\item{resistance}{A logical flag to indicate if the resistance index is calculated.}
}
\value{
A data frame cointaining the resilience value (Rs) of each site. If
more than three surveys are provided, the resilience index is calculated for
each post-disturbance survey (Rs1, Rs2...). The resistance index values are also
returned if \code{resistance = TRUE}.
}
\description{
The function \code{resilience} measures the recover capacity of
natural communities in relation to the impact of a given disturbance.
}
\details{
The resilience index is a measure of the system elasticity. It considers
the dissimilarity between the pre-disturbance and post-disturbance states of a given
community in relation to the disturbance impact (i.e. the dissimilarity between
the pre-disturbance and disturbed states). The resilience index attains a
maximum value (Rs = 1) when the pre-disturbance and the target post-disturbance state
are identical. Negative values represent systems that move away from the
pre-disturbance state.
For that, the dissimilarity measure used to generate \code{d} must be bounded between 0
and 1. Otherwise, it must be previously transformed.
At least three surveys are needed to calculate the resilience index: Pre-disturbance
(i.e. survey before disturbance), disturbed (i.e. survey during or immediately
after disturbance), and post-disturbance (i.e. survey after time enough for community
to recover to the pre-disturbance state). Additional post-disturbance states allow
calculating resilience at subsequent temporal points.
}
\references{
Sánchez-Pinillos, M., Leduc, A., Ameztegui, A. Kneeshaw, D., Lloret, F.,
Coll, L. 2019. Resistance, resilience or change: post-disturbance dynamics of boreal forests
after insect outbreaks. Ecosystems. https://doi.org/10.1007/s10021-019-00378-6
}
\seealso{
\code{\link{resistance}}
}
\author{
Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
\email{martina.sanchezpinillos@gmail.com}
}
