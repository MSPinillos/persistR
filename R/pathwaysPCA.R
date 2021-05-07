#' @title Individual pathways of post-disturbance dynamics (PCA)
#'
#' @description The function \code{pathwaysPCA} represents post-disturbance dynamics of a set
#' of disturbed communities as pathways displayed in an ordination space resulting from
#' applying principal component analyses (PCA) on the abundance matrix.
#'
#' @param abundance A data frame containing species abundance of target communities for
#' three states relative to the occurrence of a disturbance (pre-disturbance, disturbed,
#' and final states). The three surveys of each site are in rows and species
#' are in columns.
#' @param sites A vector indicating the site of each community state.
#' @param surveys A vector indicating the survey of each community state.
#' @param type A character value indicating the type of pathways to be plotted.
#' \code{"mixed"} circular dots representing pre-disturbance states are linked by lines
#' to final states represented by triangles. \code{"points"} pre-disturbance,
#' disturbed and post-disturbance states are represented by circular, squared and triangular
#' points, respectively. \code{"arrows"} post-disturbance dynamics are represented by arrows.
#' @param onlyDraw A vector indicating a subset of sites to be plotted. Sites' identifier
#' must match to those in \code{sites}.
#' @param cluster A numeric or a vector of colors indicating the cluster to which each
#' community state belongs or the color used for each state. This option is only available
#' for \code{type = "mixed"} and \code{type = "points"}
#' @param axes The pair of the principal coordinates to be used in the plot.
#' @param species A logical flag to draw species arrows.
#' @param sp.fac Factor by which species arrows will be multiplied (only for graphical purposes).
#' @param zoom A logical flag to zoom the plot to the community state coordinates.
#' @param add.legend A logical flag to add the plot legend.
#' @param cluster.name A character vector indicating the name of the cluster to which each
#' community state belongs. Values are specified in the legend if \code{add.legend = TRUE}.
#' @param ... Additional parameters to \code{plot}
#'
#' @return Individual post-disturbance dynamics of a given set of communities are displayed
#' as pathways in an ordination space (PCA) by linking the pre-disturbance, disturbed, and
#' final states in chronological order.
#'
#' @author Martina Sanchez-Pinillos (Forest Sciences Centre of Catalonia)
#' \email{martina.sanchezpinillos@@gmail.com}
#'
#' @references Sánchez-Pinillos, M., Leduc, A., Ameztegui, A. Kneeshaw, D., Lloret, F.,
#' Coll, L. 2019. Resistance, resilience or change: post-disturbance dynamics of boreal forests
#' after insect outbreaks. Ecosystems. https://doi.org/10.1007/s10021-019-00378-6
#'
#' @seealso \code{\link{avg.pathwaysPCA}}
#'
#' @examples
#'
#' data("post_disturbance")
#'
#' # Display pre-disturbance, disturbed, and final states as orange,
#' # blue, and red points, respectively, for a set of selected sites.
#' color = rep(c("orange", "blue", "red"),
#'             length(unique(post_disturbance$sites)))
#' selected.sites = c(3, 17, 66, 86, 132, 157)
#' pathwaysPCA(abundance = post_disturbance$abundance,
#'          sites = post_disturbance$sites,
#'          surveys = post_disturbance$surveys,
#'          cluster = color,
#'          type = "points", main = 'type = "points"',
#'          onlyDraw = selected.sites,
#'          zoom = TRUE)
#'
#' # type = "mixed"
#' pathwaysPCA(abundance = post_disturbance$abundance,
#'          sites = post_disturbance$sites,
#'          surveys = post_disturbance$surveys,
#'          cluster = color,
#'          type = "mixed", main = 'type = "mixed"',
#'          onlyDraw = selected.sites,
#'          zoom = TRUE)
#'
#' # If type = "arrows" a different color is used to represent
#' # post-disturbance dynamics of each site
#' pathwaysPCA(abundance = post_disturbance$abundance,
#'          sites = post_disturbance$sites,
#'          surveys = post_disturbance$surveys,
#'          type = "arrows", main = 'type = "arrows"',
#'          onlyDraw = selected.sites,
#'          zoom = TRUE)
#'
#'
#'
#'
#' @export

pathwaysPCA <- function(abundance, sites, surveys, type = "mixed", onlyDraw = NULL,
                        cluster = NULL, axes = c(1, 2), species = FALSE, sp.fac = 1,
                        zoom = FALSE, add.legend = FALSE, cluster.name = NULL, ...){

  if(!all(sapply(abundance, is.numeric))){
    stop("'abundance' must be numeric")
  }

  if(length(sites) != nrow(abundance)){
    stop("Length of 'sites' must be equal to the number of abundance rows")
  }

  if(length(surveys) != nrow(abundance)){
    stop("Length of 'surveys' must be equal to the number of abundance rows")
  }

  if(any(table(sites) != 3)){
    stop("Each site must be surveyed three times")
  }

  if(!is.null(cluster) & length(cluster) != nrow(abundance)){
    stop("Length of 'cluster' must be equal to the number of abundance rows")
  }

  # Sites that will be drawn
  if(!is.null(onlyDraw)){
    siteID <- onlyDraw
  } else {siteID <- unique(sites)}

  drawID <- which(sites %in% siteID)
  drawSites <- sites[drawID]
  drawSurveys <- surveys[drawID]

  # PCA of abundance matrix
  pca <- vegan::rda(abundance)

  # Coordinates of PCA
  scores.pca <- vegan::scores(pca, choices = axes, display = c("sites", "species"))
  coord <- as.data.frame(scores.pca$sites)
  coord.sp <- as.data.frame(scores.pca$species)

  # Set plot space if add.legend == TRUE
  if(add.legend == TRUE){
    h.leg = max(max(length(unique(cluster)), 1)*0.07, 0.2)
    graphics::layout(mat = matrix(c(1, 1, 1, 1, 0, 2), nrow = 2, ncol = 3),
                     widths = c(1, 1, 1), heights = c(1-h.leg, h.leg))
    graphics::par(mar = c(4, 4, 3, 0))
  }

  # Draw empty plot with zoom (if zoom == TRUE)
  if(zoom == TRUE){
    graphics::plot(pca, type = "n",
                   xlab = paste0("PC ", axes[1], " (", round(100*pca$CA$eig[axes[1]]/sum(pca$CA$eig)),"%)"),
                   ylab = paste0("PC ", axes[2], " (", round(100*pca$CA$eig[axes[2]]/sum(pca$CA$eig)),"%)"),
                   xlim = c(min(coord[drawID, 1]), max(coord[drawID, 1])),
                   ylim = c(min(coord[drawID, 2]), max(coord[drawID, 2])), ...)
  } else {
    graphics::plot(pca, type = "n",
                   xlab = paste0("PC ", axes[1], " (", round(100*pca$CA$eig[axes[1]]/sum(pca$CA$eig)),"%)"),
                   ylab = paste0("PC ", axes[2], " (", round(100*pca$CA$eig[axes[2]]/sum(pca$CA$eig)),"%)"), ...)
  }

  # Draw species arrows
  if(species == T){
    for(sp in 1:nrow(coord.sp)){
      shape::Arrows(x0 = 0, y0 = 0, x1 = sp.fac*(coord.sp[sp, 1]), y1 = sp.fac*(coord.sp[sp, 2]),
                    col = "darkgrey", arr.adj = 1)
      graphics::text(x = sp.fac*(coord.sp[sp, 1]), y = sp.fac*(coord.sp[sp, 2]), labels = row.names(coord.sp)[sp],
                     col = "darkgrey")
    }
  }

  for(i in 1:length(siteID)){

    # sites and surveys may not be in order
    ind <- drawID[which(drawSites == siteID[i])]
    ind <- ind[order(drawSurveys[which(drawSites == siteID[i])])]

    # Linking lines
    if(type == "mixed"){
      if(!is.null(cluster)) graphics::lines(x = coord[ind, 1], y = coord[ind, 2],
                                            col = cluster[ind][2])
      else graphics::lines(x = coord[ind, 1], y = coord[ind, 2])
    }

    # Disturbed states
    if(type == "points"){
      if(!is.null(cluster)) graphics::points(x = coord[ind[2], 1], y = coord[ind[2], 2],
                                             pch = 15, col = cluster[ind][2])

      else graphics::points(x = coord[ind[2], 1], y = coord[ind[2], 2], pch = 15)


    }

    if(type %in% c("mixed", "points")){
      # Pre-disturbance states
      if(!is.null(cluster)) graphics::points(x = coord[ind[1], 1], y = coord[ind[1], 2],
                                             pch = 16, col = cluster[ind][1])
      else graphics::points(x = coord[ind[1], 1], y = coord[ind[1], 2], pch = 16)

      # Final points
      if(!is.null(cluster)) graphics::points(x = coord[ind[3], 1], y = coord[ind[3], 2],
                                             pch = 17, col = cluster[ind][3])
      else graphics::points(x = coord[ind[3], 1], y = coord[ind[3], 2], pch = 17)
    }

    # Arrows
    if(type == "arrows"){
      graphics::lines(x = coord[ind, 1], y = coord[ind, 2], col = ind[1])
      shape::Arrows(x0 = coord[ind[2], 1], y0 = coord[ind[2], 2],
                    x1 = coord[ind[3], 1], y1 = coord[ind[3], 2],
                    col = ind[1])
    }

  }

  # Add the legend
  if(add.legend == TRUE){
    if(is.null(cluster.name)) cluster.name = rep("", length(drawSites))

    # Define legend in a dataframe
    if(!is.null(cluster)){
      drawCluster <- cluster[drawID]
      drawClustername <- cluster.name[drawID]
      legend.table <- reshape2::dcast(data.frame(drawSurveys, drawClustername, drawCluster),
                                      drawClustername + drawCluster ~ drawSurveys,
                                      value.var = "drawClustername", length)
    } else {legend.table <- data.frame(drawClustername = "", drawCluster = 1, X1 = 1, X2 = 1, X3 = 1)}

    # Convert colors into character if needed
    if(is.character(cluster)) legend.table$drawCluster <- as.character(legend.table$drawCluster)

    # Specify colors only when they are plotted
    for(i in 1:nrow(legend.table)){
      legend.table[i, 3] <- ifelse(legend.table[i, 3] > 0, legend.table$drawCluster[i], 0)
      legend.table[i, 4] <- ifelse(legend.table[i, 4] > 0, legend.table$drawCluster[i], 0)
      legend.table[i, 5] <- ifelse(legend.table[i, 5] > 0, legend.table$drawCluster[i], 0)
    }

    graphics::par(mar = c(0, 1, 0, 0))
    graphics::plot.new()

    # Draw state names
    graphics::text("Pre-disturbance\nstate", x=0.3, y=0.75, adj = c(0.5, 0.5), cex = 1)
    graphics::text("Disturbed\nstate", x=0.62, y=0.75, adj = c(0.5, 0.5), cex = 1)
    graphics::text("Final\nstate", x=0.87, y=0.75, adj = c(0.5, 0.5), cex = 1)

    # Draw symbols
    if(nrow(legend.table) <= 5){
      ypos = seq(0.5, 0.1, length.out = nrow(legend.table))
    } else {ypos = seq(0.65, 0.05, length.out = nrow(legend.table))}

    for(i in 1:nrow(legend.table)){
      graphics::text(legend.table[i, 1], x = 0.05, y = ypos[i], adj = c(0, 0), cex = 1)
      graphics::points(x = 0.3, y = ypos[i], pch = 16, col = legend.table[i, 3])
      graphics::points(x = 0.87, y = ypos[i], pch = 17, col = legend.table[i, 5])

      if(type == "points"){
        graphics::points(x = 0.62, y = ypos[i], pch = 15, col = legend.table[i, 4])
      } else if(type == "mixed"){
        graphics::points(x = 0.62, y = ypos[i], pch = "__", col = legend.table[i, 4])
      }
    }
  }

}

