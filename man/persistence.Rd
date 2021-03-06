% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/persistence.R
\name{persistence}
\alias{persistence}
\title{Persistence Index}
\usage{
persistence(
  abundance,
  RTdata,
  id.species = NULL,
  w = NULL,
  index = NULL,
  na.values = NULL
)
}
\arguments{
\item{abundance}{A data frame containing species relative abundance in target
communities. Sites are in rows and species are in columns.}

\item{RTdata}{A data frame containing numeric values of species response traits.
Species are in rows and response traits are in columns. RTdata must include response
trait values for species in abundance. Values in RTdata must be either quantitative
or binary, although quantitative data need to be standardized to the [0,1] interval.}

\item{id.species}{The name of the column in RTdata that contains species identifier.
If \code{id.species = NULL} row names of RTdata will be considered. Species identifier
(or row names) must match column names of abundance.}

\item{w}{A numerical vector listing weights for traits in RTdata. Weights in \code{w}
must follow the order in RTdata columns. If \code{w = NULL},
all traits will be considered equally relevant.}

\item{index}{A character vector listing the index or indices that the resulting data
frame should contain: \code{"PI"}, \code{"RTR"}, \code{"RTA"}}

\item{na.values}{A numerical value containing the value by which \code{NA} in RTdata
should be replaced.}
}
\value{
A data frame cointaing the value of PI, RTR, and/or RTA for each community.
If \code{index = NULL} the value of PI and its components will be returned. Otherwise,
only the specified index value (PI, RTR, or RTA) will be returned.
}
\description{
The function \code{persistence} calculates the persistence index (PI)
and its components: the response trait richness (RTR) and the response trait
abundance (RTA).
}
\details{
The persistence index is based on the assumption that an ecological community
is more resistant and/or resilient the more abundant the species presenting a given
set of desirable response traits. PI is calculated as the product of two components:
the response trait richness (RTR) quantifies the proportion of desirable traits
present in the community; and the response trait abundance (RTA) is a measure of the
relative abundance of response traits present in the species of a given community.
}
\examples{

# Load "traits" data
data("traits")

# Abundance matrix
head(traits$abundance[, 1:6])

# Response trait matrix
head(traits$RTdata)

# Calculate RTR, and RTA assumming all response traits are equally relevant
persistence(abundance = traits$abundance,
            RTdata = traits$RTdata,
            index = c("RTR", "RTA"))

# Calculate PI, RTR, and RTA considering different trait weights
persistence(abundance = traits$abundance,
            RTdata = traits$RTdata,
            w = c(1, 1, 1/3, 1/4, 1/5))

}
\references{
Sánchez-Pinillos, M., Coll, L., De Cáceres, M., Ameztegui, A. 2016.
Assessing the persistence capacity of communities facing natural disturbances on the
basis of species response traits. Ecological Indicators 66: 76-85
}
\author{
Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
\email{martina.sanchezpinillos@gmail.com}
}
