#' @title Resistance index
#'
#' @description The function \code{resistance} measures how stable a community is
#' during a given disturbance in terms of a pre-defined dissimilarity metric.
#'
#' @param d a symmetric \code{matrix} or an object of class \code{dist} cointaining
#' the dissimilarity values between pairs of community states. \code{d} values
#' must be in the range [0, 1].
#' @param sites A vector indicating the site of each community state.
#' @param surveys A vector indicating the survey of each community state.
#'
#' @details The resistance index is calculated as the one-complement of a
#' dissimilarity measure calculated to compare two community states: the
#' pre-disturbance state (i.e. before disturbance)and the disturbed state
#' (i.e. during or immediately after disturbance). The resistance index
#' ranges between 0 (maximum dissimilarity) and 1 (identical pre-disturbance
#' and disturbed states).
#' For that, the dissimilarity measure used to generate \code{d} must be
#' bounded between 0 and 1. Otherwise, it must be previously transformed.
#' If more than two surveys are provided for a given site, only the
#' first and second surveys are considered as the pre-disturbance and disturbed
#' states, respectively, to calculate resistance.
#'
#' @return A data frame cointaining the resistance value (Rt) of each site.
#'
#' @author Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
#' \email{martina.sanchezpinillos@@gmail.com}
#'
#' @references Sánchez-Pinillos, M., Leduc, A., Ameztegui, A. Kneeshaw, D., Lloret, F.,
#' Coll, L. 2019. Resistance, resilience or change: post-disturbance dynamics of boreal forests
#' after insect outbreaks. Ecosystems. https://doi.org/10.1007/s10021-019-00378-6
#'
#' @seealso \code{\link{resilience}}
#'
#' @examples
#' # Load and define data
#' data("post_disturbance")
#' abundance <- post_disturbance$abundance"
#' sites <- post_disturbance$sites
#' surveys <- post_disturbance$surveys
#'
#' # Calculate dissimilarities between states (Bray-Curtis index)
#' library(vegan)
#' dissim <- as.matrix(vegdist(x = abundance,
#'                             method = "bray"))
#'
#' # Calculate resistance
#' resistance(d = dissim,
#'           sites = sites,
#'           surveys = surveys)
#'
#' @export


resistance <- function(d, sites, surveys){

  # Convert 'd' into a matrix
  dmat <- as.matrix(d)
  if(max(dmat) > 1 | min(dmat) < 0){
    stop("'d' values must be in the range [0, 1]")
  }

  # Metrics on sites and surveys
  siteID <- unique(sites)
  nsites <- length(siteID)
  nsurv <- numeric(length(siteID))
  for(i in 1:nsites){
    nsurv[i] <- length(which(sites == siteID[1]))
  }
  max.nsurv <- max(nsurv)

  # Warnings
  if(length(sites) != nrow(dmat)){
    stop("Length of 'sites' must be equal to the number of rows/columns of 'd'")
  }

  if(max.nsurv < 2){
    stop("All sites must be surveyed at least twice")
  }

  if(any(dmat < 0 | dmat > 1)){
    stop("'d' outside [0, 1]")
  }

  if(max.nsurv > 2){
    warning("One or more sites are associated with more than two surveys. Resistance calculated through dissimilarities between the two first surveys.")
  }

  # Data frame with dissimilrity values between pairs of states
  dist.df <- matrix(NA, nrow = nsites, ncol = (max.nsurv + 1))

  for(i in 1:nsites){
    isite <- which(sites == siteID[i])

    dmat_i = dmat[isite, isite]
    ordsurv = order(surveys[isite])
    dmat_iord = dmat_i[ordsurv, ordsurv]

    dist.df[i, ] <- c(siteID[i], dmat_iord[, 1])

  }

  dist.df <- data.frame(dist.df)
  names(dist.df)[1] <- "sites"
  for(i in 2:ncol(dist.df)){
    names(dist.df)[i] <- paste0("d1", (i-1))
  }

  # Resistance
  indices <- data.frame(sites = siteID,
                        Rt = 1 - dist.df$d12)

  return(indices)

}



