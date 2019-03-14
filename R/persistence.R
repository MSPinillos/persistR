#' @title Persistence Index
#'
#' @description The function \code{persistence} calculates the persistence index (PI)
#' and its components: the response trait richness (RTR) and the response trait
#' abundance (RTA).
#'
#' @param abundance A data frame containing species relative abundance in target
#' communities. Sites are in rows and species are in columns.
#' @param RTdata A data frame containing numeric values of species response traits.
#' Species are in rows and response traits are in columns. RTdata must include response
#' trait values for species in abundance. Values in RTdata must be either quantitative
#' or binary, although quantitative data need to be standardized to the [0,1] interval.
#' @param id.species The name of the column in RTdata that contains species identifier.
#' If \code{id.species = NULL} row names of RTdata will be considered. Species identifier
#' (or row names) must match column names of abundance.
#' @param w A numerical vector listing weights for traits in RTdata. Weights in \code{w}
#' must follow the order in RTdata columns. If \code{w = NULL},
#' all traits will be considered equally relevant.
#' @param index A character vector listing the index or indices that the resulting data
#' frame should contain: \code{"PI"}, \code{"RTR"}, \code{"RTA"}
#' @param na.values A numerical value containing the value by which \code{NA} in RTdata
#' should be replaced.
#'
#' @details The persistence index is based on the assumption that an ecological community
#' is more resistant and/or resilient the more abundant the species presenting a given
#' set of desirable response traits. PI is calculated as the product of two components:
#' the response trait richness (RTR) quantifies the proportion of desirable traits
#' present in the community; and the response trait abundance (RTA) is a measure of the
#' relative abundance of response traits present in the species of a given community.
#'
#' @return A data frame cointaing the value of PI, RTR, and/or RTA for each community.
#' If \code{index = NULL} the value of PI and its components will be returned. Otherwise,
#' only the specified index value (PI, RTR, or RTA) will be returned.
#'
#' @author Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
#' \email{martina.sanchezpinillos@@gmail.com}
#'
#' @references Sánchez-Pinillos, M., Coll, L., De Cáceres, M., Ameztegui, A. 2016.
#' Assessing the persistence capacity of communities facing natural disturbances on the
#' basis of species response traits. Ecological Indicators 66: 76-85
#'
#' @examples
#'
#' # Load "traits" data
#' data("traits")
#'
#' # Abundance matrix
#' head(traits$abundance[, 1:6])
#'
#' # Response trait matrix
#' head(traits$RTdata)
#'
#' # Calculate RTR, and RTA assumming all response traits are equally relevant
#' persistence(abundance = traits$abundance,
#'             RTdata = traits$RTdata,
#'             index = c("RTR", "RTA"))
#'
#' # Calculate PI, RTR, and RTA considering different trait weights
#' persistence(abundance = traits$abundance,
#'             RTdata = traits$RTdata,
#'             w = c(1, 1, 1/3, 1/4, 1/5))
#'
#' @export

#----------------------------------------------------------------------------------------------------------

#### CALCULATION OF THE PERSISTENCE INDEX ####

persistence <- function(abundance, RTdata, id.species = NULL, w = NULL,
                        index = NULL, na.values = NULL){

  # Warning messages for abundance matrix
  if(any(rowSums(abundance) != 1)){
    warning("Total abundance must be equal to 1: defaulting to relative abundances")
    abundance <- as.data.frame(apply(abundance, 2, function(x){x = x/rowSums(abundance)}))
  }

  species <- names(abundance)[(colSums(abundance) > 0)]

  # Set site in abundance
  abundance$site <- rownames(abundance)
  sites <- abundance$site

  # Set species column in RTdata
  if(is.null(id.species)){
    RTdata$species <- rownames(RTdata)
    id.species <- "species"
  }

  traits <- names(RTdata)[names(RTdata) != id.species]

  # Warning messages for RTdata
  set.differences <- setdiff(species, RTdata[, id.species])
  if(length(set.differences) > 0){
    stop("Species names in 'abundance' and 'RTdata' should match")
  }

  if(any(is.na(RTdata)) & is.null(na.values)){
    warning("One or more values in 'RTdata' is NA: defaulting to zero")
    na.values <- 0
  }

  if(!is.null(w) & length(w) != length(traits)){
    stop("'w' lenght must be equal to the number of considered traits")
  }

  # Set parameters (M: number of traits; C: number of communities)
  M <- length(traits)
  C <- nrow(abundance)
  if(is.null(w))
    w <- rep(1, M)

  # Arrange abundance matrix to "long" data
  abundance <- reshape2::melt(abundance, id.vars = "site", variable.name = id.species,
                              value.name = "xi")
  abundance <- abundance[abundance$xi > 0, ]

  # Assigning "na.values" to missing values
  RTdata[, traits][is.na(RTdata[, traits])] <- na.values

  if(any(RTdata[, traits] < 0 | RTdata[, traits] > 1)){
    stop("'RTdata' values must be in the range [0, 1]")
  }

  # Merge species abundances and response traits by species code:
  PIdata <- merge(abundance, RTdata, by = id.species)

  #-----------------------------------------------------------------------------------

  ## RESPONSE TRAIT RICHNESS (RTR) ----

  # Empty data frame to calculate RTR

  # For each community and trait: wt*vt_max
  RTRdata <- data.frame(matrix(nrow = C, ncol = M + 1))
  for(t in 1:M){
    RTRdata[, t] <- w[t]*tapply(PIdata[, traits[t]], PIdata[, "site"], max)
  }

  # Specify ID plot and names of variables
  RTRdata[, ncol(RTRdata)] <- names(tapply(PIdata[, traits[1]], PIdata[, "site"], max))
  names(RTRdata)[1:M] <- traits
  names(RTRdata)[ncol(RTRdata)] <- "site"

  # Calculate RTR
  if(length(traits) > 1){
    RTRdata$RTR <- rowSums(RTRdata[, traits])/sum(w)
  } else {RTRdata$RTR <- RTRdata[, traits]/sum(w)}

  #-----------------------------------------------------------------------------------

  ## RESPONSE TRAIT ABUNDANCE (RTA) ----

  # For each species and trait: w*v*x matrix
  wvx <- PIdata
  RTAdata <- data.frame(matrix(nrow = length(unique(abundance[, "site"])),
                               ncol = length(traits) + 1))
  for(t in 1:M){
    wvx[, traits[t]] <- w[t]*PIdata[, traits[t]]*PIdata$xi
    RTAdata[, t] <- tapply(wvx[, traits[t]], wvx[, "site"], sum)
  }

  # Specify ID plot and names of variables
  RTAdata[, ncol(RTAdata)] <- names(tapply(wvx[, traits[1]], wvx[, "site"], sum))
  names(RTAdata)[1:M] <- traits
  names(RTAdata)[ncol(RTAdata)] <- "site"

  # Calculate RTA
  if(length(traits) > 1){
    RTAdata$RTA <- rowSums(RTAdata[, traits])/sum(w)
  } else {RTAdata$RTA <- RTAdata[, traits]/sum(w)}

  #-----------------------------------------------------------------------------------

  #### PERSISTENCE INDEX (PI) ----

  PI <- merge(RTRdata[, c("site", "RTR")], RTAdata[, c("site", "RTA")])
  PI$PI <- PI$RTR * PI$RTA

  PI <- PI[order(match(PI$site, sites)), ]
  rownames(PI) <- 1:nrow(PI)

  #-----------------------------------------------------------------------------------

  if(is.null(index)){
    return(PI)
  } else {return(PI[, c("site", index)])}

}

