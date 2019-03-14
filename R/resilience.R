#' @title Resilience index
#'
#' @description The function \code{resilience} measures the recover capacity of
#' natural communities in relation to the impact of a given disturbance.
#'
#' @param d a symmetric \code{matrix} or an object of class \code{dist} cointaining
#' the dissimilarity values between pairs of community states (pre-disturbance,
#' disturbed and post-disturbance states). \code{d} values must be in the range [0, 1].
#' @param sites A vector indicating the site of each community state.
#' @param surveys A vector indicating the survey of each community state.
#' @param Rt.quantile A numeric value of probability in [0, 1] to take the
#' corresponding \code{\link{quantile}} as threshold of the resistance index below
#' which the resilience index is calculated.
#' @param Rt.value A numeric value in [0, 1] representing the threshold of the
#' resistance index below which the resilience index is calculated.
#' @param resistance A logical flag to indicate if the resistance index is calculated.
#'
#' @details The resilience index is a measure of the system elasticity. It considers
#' the dissimilarity between the pre-disturbance and post-disturbance states of a given
#' community in relation to the disturbance impact (i.e. the dissimilarity between
#' the pre-disturbance and disturbed states). The resilience index attains a
#' maximum value (Rs = 1) when the pre-disturbance and the target post-disturbance state
#' are identical. Negative values represent systems that move away from the
#' pre-disturbance state.
#' For that, the dissimilarity measure used to generate \code{d} must be bounded between 0
#' and 1. Otherwise, it must be previously transformed.
#' At least three surveys are needed to calculate the resilience index: Pre-disturbance
#' (i.e. survey before disturbance), disturbed (i.e. survey during or immediately
#' after disturbance), and post-disturbance (i.e. survey after time enough for community
#' to recover to the pre-disturbance state). Additional post-disturbance states allow
#' calculating resilience at subsequent temporal points.
#'
#' @return A data frame cointaining the resilience value (Rs) of each site. If
#' more than three surveys are provided, the resilience index is calculated for
#' each post-disturbance survey (Rs1, Rs2...). The resistance index values are also
#' returned if \code{resistance = TRUE}.
#'
#' @author Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
#' \email{martina.sanchezpinillos@@gmail.com}
#'
#' @references Sánchez-Pinillos, M., Leduc, A., Ameztegui, A. Kneeshaw, D., Lloret, F.,
#' Coll, L. Resistance, resilience or change: post-disturbance dynamics of boreal forests
#' after insect outbreaks. (Submitted to Ecosystems)
#'
#' @seealso \code{\link{resistance}}
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
#' # Calculate resilience
#' resilience(d = dissim,
#'           sites = sites,
#'           surveys = surveys)
#'
#' # Calculate resilience for sites with a resistance value lower than the median
#' resilience(d = dissim,
#'            sites = sites,
#'            surveys = surveys,
#'            Rt.quantile = 0.5)
#'
#' # Calculate resilience and resistance for sites with a resistance value lower than 0.7
#' resilience(d = dissim,
#'            sites = sites,
#'            surveys = surveys,
#'            Rt.value = 0.7,
#'            resistance = TRUE)
#'
#' @export

resilience <- function(d, sites, surveys, Rt.quantile = NULL, Rt.value = NULL,
                       resistance = FALSE){

  # Convert 'd' into a matrix
  dmat <- as.matrix(d)

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

  if(max.nsurv < 3){
    stop("All sites must be surveyed at least three times")
  }

  if(any(dmat < 0 | dmat > 1)){
    stop("'d' outside [0, 1]")
  }

  # Data frame with dissimilrity values between pairs of states
  dist.df <- matrix(NA, nrow = nsites, ncol = (max.nsurv + 1))

  for(i in 1:nsites){
    # Select sites
    isite <- which(sites == siteID[i])
    dmat_i = dmat[isite, isite]

    # Surveys could not be ordered
    ordsurv = order(surveys[isite])
    dmat_iord = dmat_i[ordsurv, ordsurv]

    # Dissimilarities for each site
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

  # Resilience
  times_Rs <- ncol(dist.df) - 3
  for(c in 3:(2 + times_Rs)){
    indices[, c] <- 1 - dist.df[, (c+1)] / dist.df$d12
  }

  if(times_Rs == 1){
    names(indices)[3] <- "Rs"
  } else {paste0("Rs", 1:times_Rs)}

  # Subset by resistance values
  if(!is.null(Rt.quantile) & !is.null(Rt.value)){
    stop("Provide only one value for either 'Rt.quantile' or 'Rt.value'.")
  }

  if(!is.null(Rt.quantile)){
    if(Rt.quantile < 0 | Rt.quantile > 1){
      stop("'Rt.quantile' outside [0, 1]")
    }
  }

  if(!is.null(Rt.quantile)){
    Rt.value <- stats::quantile(indices$Rt, Rt.quantile)
  }

  if(!is.null(Rt.value) & (Rt.value < 0 | Rt.value > 1)){
    stop("'Rt.value' outside [0, 1]")
  }
  if(!is.null(Rt.value)){
    indices <- indices[which(indices$Rt < Rt.value), ]
  }

  if(resistance == FALSE){
    indices <- indices[, -2]
  }

  return(indices)
}



