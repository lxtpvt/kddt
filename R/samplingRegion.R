# =========================================================================================================
# sampling region
# =========================================================================================================
#' Set the sampling region
#' @param X_range a list includes the supports of covariates.
#' @param treeInfo a dataframe include the information of the DDT, default value is NULL.
#' @param id a integer denotes the id of the sampling region.
#'
#' @return a list includes the supports of all covariate at this node which is identified by it id.
#' @export
#'
setSamplingRegion <- function(X_range, treeInfo = NULL, id = NA) {
  if (is.na(id)) {
    return(list(data_range = X_range, id_range = X_range))
  } else {
    return(samplingRegion(id, X_range, treeInfo))
  }
}

#' The supports of covariates
#' @param X a dataframe includes all the observations of covariates.
#'
#' @return a list includes the names of covariates, the supports of continuous variables, the supports of discrete variables.
#' @export
#'
dataRange <- function(X) {
  if (!is.data.frame(X)) {
    return("Please input a data frame!")
  }

  namesX.continuous <- vector()
  namesX.categorical <- vector()
  sapply(X, class) -> X.class
  for (i in 1:length(X.class)) {
    if ("numeric" %in% X.class[[i]] || "integer" %in% X.class[[i]]) {
      namesX.continuous <- c(namesX.continuous, names(X.class)[i])
    } else if ("factor" %in% X.class[[i]]) {
      namesX.categorical <- c(namesX.categorical, names(X.class)[i])
    } else {
      return("Wrong format in the dataframe. numeric, integer or factor are needed.")
    }
  }
  numeric <- NULL
  factor <- NULL
  # int_levels <- function(x){as.numeric(levels(x))}
  if (length(namesX.continuous) != 0) {
    X.continuous <- data.frame(X[, namesX.continuous])
    colnames(X.continuous) <- namesX.continuous
    numeric <- sapply(X.continuous, range)
  }
  if (length(namesX.categorical) != 0) {
    X.categorical <- data.frame(X[, namesX.categorical])
    colnames(X.categorical) <- namesX.categorical
    if (length(namesX.categorical) == 1) {
      factor <- list(levels(X.categorical[, 1]))
      names(factor) <- namesX.categorical
    } else {
      factor <- sapply(X.categorical, levels)
    }
  }
  return(list(names = colnames(X), numeric = numeric, factor = factor, n = nrow(X)))
}

#' Provide the sampling region of the node with its id.
#' @param id a integer denotes the node id.
#' @param X_range a list includes the covariates supports.
#' @param treeInfo a dataframe includes the tree information.
#'
#' @return a list includes the sampling region of the node with its id.
#' @export
#'
samplingRegion <- function(id, X_range, treeInfo) {
  id_X_range <- X_range
  df_conditions <- rangeConditions(id, treeInfo)
  df_numeric <- df_conditions$numeric
  df_factor <- df_conditions$factor
  x <- sapply(df_factor, is.factor)
  df_factor[x] <- lapply(df_factor[x], as.character)
  # categorical variables
  if (dim(df_factor)[1] > 0 && length(X_range$factor) > 0) {
    for (i in 1:dim(df_factor)[1]) {
      for (j in 1:length(X_range$factor)) {
        if (df_factor$var[i] == names(X_range$factor)[j]) {
          unlist(strsplit(df_factor$split[i], ",")) -> tp_i
          id_X_range$factor[[j]] <- intersect(id_X_range$factor[[j]], tp_i)
        }
      }
    }
  }
  # continuous variables
  if (dim(df_numeric)[1] > 0 && dim(X_range$numeric)[1] > 0) {
    for (i in 1:dim(df_numeric)[1]) {
      for (j in 1:dim(X_range$numeric)[2]) {
        if (df_numeric$var[i] == colnames(X_range$numeric)[j]) {
          if (df_numeric$sign[i] %in% c(">=", ">")) {
            if (X_range$numeric[1, j] < df_numeric$split[i]) {
              id_X_range$numeric[1, j] <- df_numeric$split[i]
            }
          } else if (df_numeric$sign[i] %in% c("<=", "<")) {
            if (X_range$numeric[2, j] > df_numeric$split[i]) {
              id_X_range$numeric[2, j] <- df_numeric$split[i]
            }
          }
        }
      }
    }
  }
  return(list(data_range = X_range, id_range = id_X_range))
}
