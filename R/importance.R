#
#' @param X a dateframe includes the obsercations of covariates.
#' @param y a vector for regression or classification, NULL for survival.
#' @param tree a rpart object is the fitted tree.
#' @param type a strings in the vector of c("classification","regression","survival").
#'
#' @return a list includes the variable importance of the interpretable part of kddt.
#' @export
#'
importance.interpretable <- function(X, y, tree, type)
{
  if(type == "survival"){
    if(length(system.file(package='survival'))==0){
      return("The package 'survival' is needed.")
    }
    y = tree$survival.y
  }
  n.y = length(y)
  treeInfo = treeInfo(tree)

  impSplitList <- tree$impList
  n <- length(impSplitList)
  split.nid = NULL
  split.var = NULL
  split.delta = NULL
  split.dataNumber = NULL
  imp.splits = list()
  for (i in 1:n) {
    # index of data points in parent node
    impSplitList[[i]]$sub_ids_nids[,'ids'] -> ids.parents
    n.parent <- length(ids.parents)
    if(n.parent>1){ # the split node must have real data points > 1

      temp.df <- as.data.frame(impSplitList[[i]]$sub_ids_nids)

      split.dataNumber = c(split.dataNumber, n.parent)
      # split nid
      snid = impSplitList[[i]]$nid
      split.nid = c(split.nid, snid)
      # split variable name
      impSplitList[[i]]$split_name -> splitVarName
      split.var = c(split.var, splitVarName)
      # real y values in parent node
      temp.df$actual.parent <- y[ids.parents]
      # predict values in parent node
      temp.df$predicted.parent <- tree$frame[which(treeInfo$nid==snid),"yval"]

      # initial children predicted
      which(temp.df$where==2*snid)->ids.lchild
      which(temp.df$where==(2*snid+1))->ids.rchild

      temp.df$predicted.children <- 0

      if(length(ids.lchild)>0){
        temp.df$predicted.children[ids.lchild] <- tree$frame[which(treeInfo$nid==2*snid),"yval"]
      }
      if(length(ids.rchild)>0){
        temp.df$predicted.children[ids.rchild] <- tree$frame[which(treeInfo$nid==(2*snid+1)),"yval"]
      }

      inpurity.delta = 0
      if(type %in% c("regression","survival")){ # if split variable is numeric
        inpurity.parent.delta <- mse(temp.df$actual.parent, temp.df$predicted.parent)
        inpurity.children.delta <- mse(temp.df$actual.parent, temp.df$predicted.children)
      }else if(type == "classification"){ # if split variable is categorical
        inpurity.parent.delta <- gini_index(temp.df$actual.parent, temp.df$predicted.parent)
        inpurity.children.delta <- gini_index(temp.df$actual.parent, temp.df$predicted.children)
      }
      inpurity.delta = abs((inpurity.parent.delta-inpurity.children.delta)*(n.parent/n.y))
      split.delta = c(split.delta, inpurity.delta)

      split.imp.df = impSplitList[[i]]$imp.df
      split.imp.df$imp.percent <- split.imp.df$imp.percent*inpurity.delta
      imp.splits = append(imp.splits, list(split.imp.df))
    }
  }
  # all the information for calculating variable importance
  inpurity.reduction <- data.frame(split.nid, split.var, split.delta, split.dataNumber)
  return(list(inpurity.reduction=inpurity.reduction, imp.splits=imp.splits))
}

# Fit a large number of stumps
#' Insert the child tree's rpart.object$frame into its parent tree's rpart.object$frame.
#' @param kddt a rpart object is the fitted tree.
#' @param type a strings in the vector of c("classification","regression","survival").
#'
#' @return a list includes the variable importance of the kddt.
#' @export
#'
importance.kddt <- function(kddt, type){

  if(!(type %in% c("classification","regression","survival"))){
    return("Please select 'type' in the set (classification, regression, survival)")
  }

  X = kddt$data$X
  y = kddt$data$y
  interpretableTree=kddt$interpretableTree
  treeInfo = treeInfo(interpretableTree)
  predictiveSubTrees=kddt$predictiveSubTrees

  if(type == "survival"){
    if(length(system.file(package='survival'))==0){
      return("The package 'survival' is needed.")
    }
    y = kddt$interpretableTree$survival.y
  }

  # the prediction of the entire kddt
  if(type %in% c("regression","survival")){
    y.kddt <- kddtPredict(kddt, X, "vector")
  }else{
    y.kddt <- kddtPredict(kddt, X, "class")
  }

  imp.interpretableTree = importance.interpretable(X, y, interpretableTree, type)

  interpretabilityIndex = imp.interpretableTree$inpurity.reduction

  n.y = length(y)

  # inpurity.subTrees.delta = NULL
  imp.subTrees = list()

  oids <- interpretableTree$where
  nids <- predictiveSubTrees$ids

  for (nid in nids) {
    temp.oids = 0
    temp.oids <- which(oids==nid)
    if(length(temp.oids)>0){
      actual <- y[temp.oids]
      interpretable.predict <- interpretableTree$frame[which(treeInfo$nid==nid),"yval"]
      predictive.predict <- y.kddt[temp.oids]

      inpurity.delta = 0
      if(type %in% c("regression","survival")){ # if split variable is numeric
        inpurity.interpretable.delta <- mse(actual, interpretable.predict)
        inpurity.predictive.delta <- mse(actual, predictive.predict)
      }else if(type == "classification"){ # if split variable is categorical
        inpurity.interpretable.delta <- gini_index(actual, interpretable.predict)
        inpurity.predictive.delta <- gini_index(actual, predictive.predict)
      }
      inpurity.delta = abs((inpurity.interpretable.delta-inpurity.predictive.delta)*(length(temp.oids)/n.y))
      # inpurity.subTrees.delta = c(inpurity.subTrees.delta, inpurity.delta)

      predictiveSubTrees$trees[[which(nids==nid)]]$variable.importance -> imp.subTree
      splitCovNames = names(imp.subTree)
      imp.percent = (imp.subTree/sum(imp.subTree))*inpurity.delta
      subTree.imp.df = data.frame(splitCovNames,imp.percent)
      imp.subTrees = append(imp.subTrees, list(subTree.imp.df))

      interpretabilityIndex[nrow(interpretabilityIndex) + 1,] <- list(nid,"<leaf>",inpurity.delta,length(temp.oids))
    }
  }

  interpretabilityIndex$split.delta = interpretabilityIndex$split.delta/sum(interpretabilityIndex$split.delta)
  names(interpretabilityIndex)<-c('id','variable','(P)XP','obs')

  imp.dfs = append(imp.interpretableTree$imp.splits, imp.subTrees)

  dataRange(X)$names -> variable.names
  variable.importance = numeric(length(variable.names))
  for (i in 1:length(variable.names)) {
    for (imp.df in imp.dfs) {
      if(variable.names[i] %in% imp.df$splitCovNames){
        variable.importance[i] = variable.importance[i] + imp.df$imp.percent[which(imp.df$splitCovNames==variable.names[i])]
      }
    }
  }
  variable.importance = abs(variable.importance/sum(variable.importance))
  df.all = data.frame(variable.names,variable.importance)
  df.all = df.all[order(df.all$variable.importance, decreasing = T),]
  res=list(importance = df.all, XIs=interpretabilityIndex)
  return(res)
}









