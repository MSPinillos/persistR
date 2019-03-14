#' @title Average pathways of post-disturbance dynamics
#'
#' @description The function \code{avg.pathways} summarizes and represents average
#' post-disturbance dynamics of a set of disturbed communities as pathways displayed
#' in an ordination space resulting from applying principal component analyses (PCA)
#' on an abundance matrix.
#'
#' @param abundance A data frame containing species abundance of target communities
#' for three states relative to the occurrence of a disturbance (pre-disturbance,
#' disturbed, and post-disturbance states). The three surveys of each site are in
#' rows and species are in columns.
#' @param sites A vector indicating the site of each community state.
#' @param surveys A vector indicating the survey of each community state.
#' @param cluster A numeric or a vector of colors indicating the cluster to which
#' each community state belongs or the color used for each state.
#' @param type A character value indicating the type of pathways to be plotted.
#' \code{"mixed"} circular dots representing pre-disturbance states are linked by
#' lines to post-disturbance states represented by triangles. \code{"points"}
#' pre-disturbance, disturbed and post-disturbance states are represented by circular,
#' squared and triangular points, respectively. \code{"arrows"} post-disturbance
#' dynamics are represented by arrows.
#' @param onlyPaths A vector of length three containing the cluster identifier of
#' each state (i.e. pre-disturbance, disturbed, and post-disturbance states) to be
#' plotted. If all clusters are considered for one or two states, corresponding value/s
#' must be replaced by \code{NA}. If all pathways must be considered, set \code{onlyPaths = NULL}.
#' @param lwd.factor A positive number to use line widths proportional to the number
#' of plots following the same pathways. The default \code{NULL} uses the same \code{lwd = 1}
#' for all pathways.
#' @param axes The pair of the principal coordinates to be used in the plot.
#' @param species A logical value indicating if species arrows will be plotted.
#' @param zoom A logical value to zoom the plot to the community state coordinates.
#' @param add.legend A logical flag to add the plot legend.
#' @param cluster.name A character vector indicating the name of the cluster to which each
#' community state belongs. Values are specified in the legend if \code{add.legend = TRUE}.
#' @param ... Additional parameters to \code{plot}
#'
#' @details Average pathways are constructed by linking in chronological order the
#' centroids calculated for each group of pre-disturbance, disturbed, and post-disturbance
#' states of the set of communities following equivalent dynamics (i.e. exhibiting
#' the same combination of clusters along the three states). Clusters corresponding
#' to each state and community must be previously defined and provided through the
#' argument \code{clusters}. For \code{type = "mixed} or \code{type = "points}
#' pre-disturbance, disturbed, and post-disturbance average states are represented
#' by the color corresponding to \code{cluster}.
#'
#'
#' @return Average post-disturbance dynamics of a set of communities are displayed
#' as pathways represented in an ordination space (PCA) by linking the pre-disturbance,
#' disturbed, and post-disturbance states in chronological order.
#'
#' @author Martina Sánchez-Pinillos (Forest Sciences Centre of Catalonia)
#' \email{martina.sanchezpinillos@@gmail.com}
#'
#' @references Sánchez-Pinillos, M., Leduc, A., Ameztegui, A. Kneeshaw, D., Lloret, F.,
#' Coll, L. Resistance, resilience or change: post-disturbance dynamics of boreal forests
#' after insect outbreaks. (Submitted to Ecosystems)
#'
#' @seealso \code{\link{pathways}}
#'
#' @examples
#'
#' data("post_disturbance")
#'
#' # Display individual pathways using a different color for each
#' # cluster
#' head(post_disturbance$cluster)
#' pathways(abundance = post_disturbance$abundance,
#'          sites = post_disturbance$sites,
#'          surveys = post_disturbance$surveys,
#'          cluster = post_disturbance$cluster,
#'          type = "mixed",
#'          main = "Individual pathways",
#'          zoom = TRUE)
#'
#' # Summarize individual into average pathways
#' avg.pathways(abundance = post_disturbance$abundance,
#'              sites = post_disturbance$sites,
#'              surveys = post_disturbance$surveys,
#'              cluster = post_disturbance$cluster,
#'              type = "mixed",
#'              main = "Average pathways",
#'              lwd.factor = 20,
#'              zoom = TRUE)
#'
#' # Display average dynamics of communities classified in cluster 4
#' # in the pre-disturbance state.
#' avg.pathways(abundance = post_disturbance$abundance,
#'              sites = post_disturbance$sites,
#'              surveys = post_disturbance$surveys,
#'              cluster = post_disturbance$cluster,
#'              type = "mixed",
#'              onlyPaths = c(4, NA, NA),
#'              main = "Average pathways of 'cluster 4'",
#'              lwd.factor = 5,
#'              zoom = TRUE)
#'
#'
#' @export
#'


avg.pathways <- function(abundance, sites, surveys, cluster, type = "mixed",
                         onlyPaths = NULL, lwd.factor = NULL,
                         axes = c(1, 2), species = FALSE, zoom = FALSE,
                         add.legend = FALSE, cluster.name = NULL, ...){

  if(length(sites) != nrow(abundance)){
    stop("Length of 'sites' must be equal to the number of abundance rows")
  }

  if(length(surveys) != nrow(abundance)){
    stop("Length of 'surveys' must be equal to the number of abundance rows")
  }

  if(length(cluster) != nrow(abundance)){
    stop("Length of 'cluster' must be equal to the number of abundance rows")
  }

  if(any(table(sites) != 3)){
    stop("Each site must be surveyed three times")
  }



  # PCA of abundance matrix
  pca <- vegan::rda(abundance)

  # Coordinates of PCA
  scores.pca <- vegan::scores(pca, choices = axes, display = c("sites", "species"))
  coord <- as.data.frame(scores.pca$sites)
  coord.sp <- as.data.frame(scores.pca$species)



  # Define pathway codes
  sit.surv <- data.frame(sites = sites, surveys = surveys, cluster = cluster)
  cluster.df <- reshape2::dcast(
    sit.surv,
    sites ~ surveys, value.var = "cluster"
  )
  cluster.df$pathways <- apply(cluster.df[, -1], 1, paste, collapse = "")
  cluster.df <- merge(sit.surv[, 1:2], cluster.df, all.x = T)
  cluster.df <- cluster.df[order(match(
    paste(cluster.df$sites, cluster.df$surveys),
    paste(sit.surv$sites, sit.surv$surveys)
  )), ]

  # Pathways following each site and survey
  pathways <- cluster.df$pathways

  # Subset selected pathways
  if(!is.null(onlyPaths)){
    for(i in 1:length(onlyPaths)){
      if(!is.na(onlyPaths[i])){
        cluster.df <- cluster.df[cluster.df[, (i+2)] == onlyPaths[i], ]
      }
    }
  }

  # Proportion of plots for each pathway to set lwd
  if(!is.null(lwd.factor)){
    prop <- lwd.factor * table(unique(cluster.df)$pathways)/sum(table(unique(cluster.df)$pathways))
  } else {prop <- table(unique(cluster.df)$pathways)/table(unique(cluster.df)$pathways)}

  # Calculate average coordinates by state and cluster
  avg.coord <- reshape2::melt(unique(cluster.df[, 3:ncol(cluster.df)]),
                              id = "pathways",
                              variable.name = "surveys", value.name = "cluster")
  avg.coord <- avg.coord[order(avg.coord$pathways, avg.coord$surveys),]

  for(i in 1:nrow(avg.coord)){
    id <- which(pathways == avg.coord$pathways[i] & surveys == avg.coord$surveys[i])
    avg.coord$x[i] <- mean(coord[id, 1])
    avg.coord$y[i] <- mean(coord[id, 2])
  }
  avg.coord <- merge(avg.coord, data.frame(prop),
                     by.x = "pathways", by.y = "Var1", all.x = T)

  # Pathways that will be drawn
  pathID <- unique(cluster.df$pathways)
  drawID <- which(avg.coord$pathways %in% pathID)
  drawPath1 <- avg.coord$pathways[drawID]
  drawSurveys <- avg.coord$surveys[drawID]
  drawCluster <- avg.coord$cluster[drawID]

  # Set plot space if add.legend == TRUE
  if(add.legend == TRUE){
    graphics::layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2, byrow = T), widths = c(1, 1/2))
  }

  graphics::par(mar = c(4, 4, 3, 0))

  # Draw empty plot with zoom (if zoom = TRUE)
  if(zoom == TRUE){
    graphics::plot(pca, type="n",
                   xlab = paste0("PC ", axes[1], " (", round(100*pca$CA$eig[axes[1]]/sum(pca$CA$eig)),"%)"),
                   ylab = paste0("PC ", axes[2], " (", round(100*pca$CA$eig[axes[2]]/sum(pca$CA$eig)),"%)"),
                   xlim = c(min(avg.coord[drawID, "x"]), max(avg.coord[drawID, "x"])),
                   ylim = c(min(avg.coord[drawID, "y"]), max(avg.coord[drawID, "y"])), ...)
  } else {
    graphics::plot(pca, type = "n",
                   xlab = paste0("PC ", axes[1], " (", round(100*pca$CA$eig[axes[1]]/sum(pca$CA$eig)),"%)"),
                   ylab = paste0("PC ", axes[2], " (", round(100*pca$CA$eig[axes[2]]/sum(pca$CA$eig)),"%)"), ...)
  }

  # Draw species arrows
  if(species == T){
    for(sp in 1:nrow(coord.sp)){
      shape::Arrows(x0 = 0, y0 = 0, x1 = coord.sp[sp, 1], y1 = coord.sp[sp, 2],
                    col = "darkgrey", arr.adj = 1)
      graphics::text(x = coord.sp[sp, 1], y = coord.sp[sp, 2], labels = row.names(coord.sp)[sp],
                     col = "darkgrey")
    }
  }


  for(i in 1:length(pathID)){

    # sites and surveys may not be in order
    ind = drawID[which(drawPath1 == pathID[i])]
    ind = ind[order(drawSurveys[which(drawPath1 == pathID[i])])]

    # Linking lines
    if(type == "mixed"){
      graphics::lines(x = avg.coord[ind, "x"], y = avg.coord[ind, "y"],
                      col = avg.coord$cluster[ind][2], lwd = avg.coord$Freq[ind][1])
    }

    # Disturbed states
    if(type == "points"){
      graphics::points(x = avg.coord[ind[2], "x"], y = avg.coord[ind[2], "y"],
                       pch = 15, col = avg.coord$cluster[ind][2])
    }

    if(type %in% c("mixed", "points")){
      # Pre-disturbance states
      graphics::points(x = avg.coord[ind[1], "x"], y = avg.coord[ind[1], "y"],
                       pch = 16, col = avg.coord$cluster[ind][1])

      # Post disturbance points
      graphics::points(x = avg.coord[ind[3], "x"], y = avg.coord[ind[3], "y"],
                       pch = 17, col = avg.coord$cluster[ind][3])
    }

    # Arrows
    if(type == "arrows"){
      graphics::lines(x = avg.coord[ind, "x"], y = avg.coord[ind, "y"],
                      lwd = avg.coord$Freq[ind][1], col = ind[1])
      shape::Arrows(x0 = avg.coord[ind[2], "x"], y0 = avg.coord[ind[2], "y"],
                    x1 = avg.coord[ind[3], "x"], y1 = avg.coord[ind[3], "y"],
                    arr.lwd = avg.coord$Freq[ind][1], col = ind[1])
    }

  }

  # Add the legend
  if(add.legend == TRUE){
    drawSites <- sites[drawID]
    if(is.null(cluster.name)) cluster.name = rep("", length(drawSites))

    # Define legend in a dataframe
    drawClustername <- cluster.name[drawID]
    legend.table <- reshape2::dcast(data.frame(drawSurveys, drawClustername, drawCluster),
                                    drawClustername + drawCluster ~ drawSurveys,
                                    value.var = "drawClustername", length)

    # Convert colors into character if needed
    if(is.character(cluster)) legend.table$drawCluster <- as.character(legend.table$drawCluster)

    # Specify colors only when they are plotted
    for(i in 1:nrow(legend.table)){
      legend.table[i, 3] <- ifelse(legend.table[i, 3] > 0, legend.table$drawCluster[i], 0)
      legend.table[i, 4] <- ifelse(legend.table[i, 4] > 0, legend.table$drawCluster[i], 0)
      legend.table[i, 5] <- ifelse(legend.table[i, 5] > 0, legend.table$drawCluster[i], 0)
    }

    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()

    # Draw state names
    graphics::text("Pre-disturbance\nstate", x=0.3, y=0.75, adj = c(0.5, 0.5), cex = 0.7)
    graphics::text("Disturbed\nstate", x=0.62, y=0.75, adj = c(0.5, 0.5), cex = 0.7)
    graphics::text("Final\nstate", x=0.87, y=0.75, adj = c(0.5, 0.5), cex = 0.7)

    # Draw symbols
    ypos = seq(0.65, 0.05, length.out = nrow(legend.table))
    for(i in 1:nrow(legend.table)){
      graphics::text(legend.table[i, 1], x = 0.05, y = ypos[i], adj = c(0, 0), cex = 0.7)
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

