
#' Plot the first-class stability of a node
#' @param nid a integer denotes the node id.
#' @param firstLevelStability a dateframe includes the pmf of covariates.
#'
#' @return a barplot show the pmf of all covariates.
#' @export
#'
plotFirstLevelStability <- function(nid, firstLevelStability) {
  # par(mfrow = c(1,1))
  barplot(firstLevelStability, ylim = c(0, 1), main = paste("Split id: ", as.character(nid)))
}

#
#' Plot the second-class stability of a node
#' @param nid a integer denotes the node id.
#' @param secondLevelStability a list includes all the information related to the second level stability.
#'
#' @return a density/barplot show the pdf/pmf of all covariates.
#' @export
#'
plotSecondLevelStability <- function(nid, secondLevelStability) {

  strMain <- paste("Variable name: ", secondLevelStability$name, ",  Split id: ", as.character(nid))

  if (secondLevelStability$isNumeric) {
    plot(secondLevelStability$density,
      lwd = 2, col = "red",
      main = strMain,
      cex.lab = 1.5, cex.axis = 1.5
    )
    abline(v = secondLevelStability$max)
    legend("topleft",
      legend = parse(text = sprintf(
        "Mode == %s",
        round(secondLevelStability$max, 2)
      )), bty = "n", cex = 1.5
    )
  } else {
    barplot(secondLevelStability$prop.table,
      ylim = c(0, 1), main = strMain, cex = 1.5
    )
  }
  # par(mfrow = c(1,1))
}

#====================================================================================
# legacy code
#====================================================================================
plotSplitStability <- function(node_list, nid, nTopCovs = 0, numericVarsName = NULL) {
  n <- sum(node_list[[1]]$pmf)
  n.c <- length(node_list[[1]]$pmf)
  node_list[[as.character(nid)]] -> node
  if (is.null(node)) {
    return("nid is not valid")
  }

  if (nTopCovs < 0 | nTopCovs > n.c) {
    return("nTopCovs is out of bounds.")
  } else if (nTopCovs > 0 && !is.null(numericVarsName)) {
    cov_names <- names(node)[-1]
    top_names <- (colnames(node$pmf)[order(node$pmf, decreasing = T)])[1:nTopCovs]
    name_set <- top_names[match(intersect(top_names, cov_names), top_names)]
    par(mfrow = c(ceiling((length(name_set) + 1) / 2), 2))
    barplot(node$pmf / n, ylim = c(0, 1), main = paste("Node id: ", as.character(nid)))
    for (nm in name_set) {
      temp_mat <- do.call(rbind, node[[nm]])
      if (nm %in% numericVarsName) {
        x <- as.numeric(temp_mat[, 2])
        if (length(x) > 1) {
          dx <- density(x)
          plot(dx, lwd = 2, col = "red", main = nm)
          abline(v = dx$x[which.max(dx$y)])
          legend("topright", legend = parse(text = sprintf(
            "Mode == %s",
            round(dx$x[which.max(dx$y)], 2)
          )), bty = "n")
        }
      } else {
        barplot(prop.table(table(temp_mat[, 2])), main = nm, ylim = c(0, 1))
      }
    }
  } else {
    par(mfrow = c(1, 1))
    barplot(node$pmf / n, ylim = c(0, 1), main = paste("Node id: ", as.character(nid)))
  }
  par(mfrow = c(1, 1))
}

# Find the stable sample size
plotStableSize <- function(distance_list, name_vector, distance_type, title) {
  par(mfrow = c(1, 1))
  n <- length(distance_list)
  distance <- distance_list[[1]]
  plot(distance$sampleSize, distance$distance,
       xlab = "Sample Size",
       ylab = paste("Distance: ", distance_type),
       xlim = c(0, max(distance$sampleSize) * 1.1),
       ylim = c(0, max(distance$distance) * 1.1),
       main = title,
       type = "o"
  )
  if (n > 1) {
    for (i in 2:n) {
      lines(distance_list[[i]]$sampleSize, distance_list[[i]]$distance, col = i)
      points(distance_list[[i]]$sampleSize, distance_list[[i]]$distance, col = i)
    }
  }
  legend("topright", name_vector,
         lty = rep(1, n),
         col = c(1:n)
  )
}

