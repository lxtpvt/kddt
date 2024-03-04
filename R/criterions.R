#' Criterion
#'
#' @param y
#' @param pred_y
#' @param type
#'
#' @return MSE or mis-classification rate
#' @export
#'
criterion <- function(y, pred_y, type) {
  # assume y and pred_y are vectors
  n <- length(y)
  if(type){
    return(sum((y - pred_y)^2) / n)
  }else{
    return(sum(y == pred_y) / n)
  }
}

gini_index <- function(actual,predicted){

  # Calculate the confusion matrix
  cm <- table(actual, predicted)

  # Calculate the total number of observations
  n <- sum(cm)

  # Calculate the proportion of observations in each class
  prop_class <- apply(cm, 2, sum) / n

  # Calculate and return the Gini index
  return(1 - sum(prop_class^2))

}


#' @export
sse <- function(actual, predicted){
  return(sum((actual-predicted)^2)) # SSE
}


mse <- function(actual,predicted){
  return(sum((actual-predicted)^2)/length(actual))  # MSE
}

cIndex <- function(predicted, time, status){
  return(1-survival::survConcordance(survival::Surv(time, status) ~ predicted)$concordance)
}
