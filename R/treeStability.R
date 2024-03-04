
#===============================================================================
# Single stump stability for tree induction
#===============================================================================

# analysis split (stump) stability
# extract simulation results to a matrix

stumpsToMat <- function(stumps_list){
  if(length(stumps_list)==1){
    stumps_list = stumps_list[[1]]
  }
  simRes = list()
  for (stump in stumps_list) {
    #print(stump)
    treeInfo(stump) -> stumpInfo
    var = as.character(stumpInfo$var[1])
    condition = stumpInfo$conditions[2]
    parseCondition(condition)->aaa
    split = aaa[1]
    simRes = append(simRes, list(c(var = var, split = split)))
  }
  simResMatStump = do.call(rbind, simRes)
  return(simResMatStump)
}

# create the pmf of covariates
firstLevelStabilityStump <- function(namesCovariates, simResMatStump){

  length(namesCovariates) -> p
  pmf_mat = matrix(0,nrow = 1,ncol = p)
  colnames(pmf_mat)<-namesCovariates
  n = dim(simResMatStump)[1]
  if(is.null(simResMatStump) | n==0){
    return("simResMatStump can't be empty!")
  }
  for (i in 1:p) {

    id_flag = (simResMatStump[,"var"]==namesCovariates[i])
    pmf_mat[1,namesCovariates[i]] = sum(id_flag)/n

  }
  return(pmf_mat)

}

# Select the optimal split value according the split variable and the second stability.
# Need to return a stump (rpart.object).
secondLevelStabilityStump <- function(nameCov, isNumeric, simResMatStump){
  if(is.null(simResMatStump)){
    return("simResMat can't be empty!")
  }
  n = dim(simResMatStump)[1]
  ids = (simResMatStump[,"var"]==nameCov)
  if(isNumeric){
    if(n>1){
      x = as.numeric(simResMatStump[ids,"split"])
      dx <- density(x)
      max = dx$x[which.max(dx$y)]
      # find the nearest one
      abs_diff <- abs(x - max)
      # if there are more than one nearest values, select the first one.
      id = which(simResMatStump[,"split"]==as.character(x[which.min(abs_diff)]))[1]
    }else{
      # if n=1
      dx = NULL
      max=simResMatStump[1, "split"]
      id = 1
    }
    return(list(name=nameCov, isNumeric=isNumeric, density=dx, max=max, stumpId=id))
  }else{
    if(n>1){
      prop.table(table(simResMatStump[ids,"split"])) -> a
      maxCat = names(which.max(a))[1]
      id = which(simResMatStump[,"split"]==maxCat)[1]
    }else{
      # if n=1
      a = NULL
      max=simResMatStump[1, "var"]
      id = 1
    }
    return(list(name=nameCov, isNumeric=isNumeric, prop.table=a, max=maxCat, stumpId=id))
  }
}




#===============================================================================
# The old functions, please ignore them!
#===============================================================================
stableSampleSize <- function(fitedModel, samplingStrategy, rpartParas, sampleSize_v, criterion){

  m = length(sampleSize_v)
  distance = numeric(m)

  for (i in 1:m) {
    # training data
    X_tr = do.call(samplingStrategy$samplingMethod,
                   list(samplingStrategy$samplingRegion,
                        samplingStrategy$samplingParameters, sampleSize_v[i]))
    y_tr = do.call("predict",list(fitedModel,X_tr))
    df_sm = data.frame(X_tr,y_tr)
    # fit tree
    tree_m = rpart(y_tr~., data = df_sm, method = rpartParas$method,
                   control = rpartParas$control)
    # testing data
    X_te = do.call(samplingStrategy$samplingMethod,
                   list(samplingStrategy$samplingRegion,
                        samplingStrategy$samplingParameters, floor(sampleSize_v[i]*0.3)))
    y_te = do.call("predict",list(fitedModel,X_te))
    # predit
    pred = predict(object = tree_m, newdata = X_te, type = rpartParas$predict_type)
    # calculate distance
    distance[i] = do.call(criterion, list(y_te, pred, rpartParas$predict_type))
  }
  return(data.frame(sampleSize = sampleSize_v, distance = distance))
}


# do simulations: fit distillation tree
simStability <- function(nSim, fitedModel, samplingStrategies, sampleSize, rpartParas){

  nStrategies = length(samplingStrategies)
  length(samplingStrategies[[1]]$samplingRegion$data_range$names) -> p
  samplingStrategies[[1]]$samplingRegion$data_range$names -> namesCov

  spT = list()

  for (i in 1:nSim) {
    X = NULL
    for (j in 1:nStrategies) {
      temp_X = do.call(samplingStrategies[[j]]$samplingMethod,
                       list(samplingStrategies[[j]]$samplingRegion,
                            samplingStrategies[[j]]$samplingParameters, sampleSize[[j]]))
      X = rbind(X, temp_X)
    }
    y = do.call("predict",list(fitedModel,X))
    df = data.frame(X,y)
    tree = rpart(y~., data = df, method = rpartParas$method, control = rpartParas$control)
    spT = append(spT, list(tree))
    print(i)
  }

  return(simRes = list(treeTables = spT, namesCov=namesCov))

}


# analysis tree stability
stabilityAnalyze <- function(simRes){

  length(simRes$namesCov) -> p
  pmf_mat = matrix(0,nrow = 1,ncol = p)
  colnames(pmf_mat)<-simRes$namesCov
  nid_v = c()
  sp_list = list()
  for (tree in simRes$treeTables) {
    sp = splitTable(treeInfo(tree))
    nid_v = union(nid_v,sp$nid)
    sp_list = append(sp_list, list(sp))
  }
  nodes_list = split(nid_v,nid_v)

  for (i in 1:length(nid_v)) {
    nodes_list[[i]] = list(pmf = pmf_mat)
  }

  tree_in_sim = 1
  for (sp in sp_list) {
    n = nrow(sp)
    for (j in 1:n) {
      ind = which(colnames(nodes_list[[as.character(sp$nid[j])]]$pmf)==sp$var[j])
      nodes_list[[as.character(sp$nid[j])]]$pmf[1,ind] =
        nodes_list[[as.character(sp$nid[j])]]$pmf[1,ind] + 1

      if(sp$var[j] %in% names(nodes_list[[as.character(sp$nid[j])]])[-1]){
        temp_id = which(names(nodes_list[[as.character(sp$nid[j])]])==sp$var[j])
        nodes_list[[as.character(sp$nid[j])]][[temp_id]] =
          append(nodes_list[[as.character(sp$nid[j])]][[temp_id]],list(c(tree_in_sim,sp$split[j])))
      }else{
        nodes_list[[as.character(sp$nid[j])]] =
          append(nodes_list[[as.character(sp$nid[j])]],list(l = list(c(tree_in_sim,sp$split[j]))))
        names(nodes_list[[as.character(sp$nid[j])]])[length(names(nodes_list[[as.character(sp$nid[j])]]))]=sp$var[j]
      }
    }
    tree_in_sim = tree_in_sim+1
  }

  return(list(node_list = nodes_list, split_list = sp_list))

}

# Find the stable splits by a criterion (cutoff of the pmf of covariates)
#      and index of trees in the simulation

stableSplits <- function(node_list, criteria){

  n_sim = sum(node_list[[1]][[1]])
  n_node = length(names(node_list))
  stable_node_list = list()
  stable_tree_list = list()
  for (i in 1:n_node) {
    node_list[[i]][[1]]/n_sim->pmf
    max_id = which.max(pmf[1,])
    if(pmf[1,max_id]>=criteria){
      stable_node_list = append(stable_node_list, names(node_list[i]))

      max_name = colnames(pmf)[max_id]
      node_list[[i]][[max_name]]->temp_stb
      stable_tree_list = append(stable_tree_list, list(as.integer(as.data.frame(temp_stb)[1,])))
    }
  }
  return(list(stableNodes = stable_node_list, stableTrees = stable_tree_list))
}

# Find the list of stable trees satisfied the criteria.
stableTrees <- function(stableSplits_stableTrees, simRes=NULL){

  stableSplits_stableTrees[[1]] -> temp_ids
  for (i in 2:length(stableSplits_stableTrees)) {
    temp_ids = intersect(temp_ids,stableSplits_stableTrees[[i]])
  }
  if(is.null(simRes)){
    return(list(ids = temp_ids))
  }else{
    return(list(ids = temp_ids, trees = simRes[[1]][temp_ids]))
  }
}

minStableTrees <- function(stableSplits, simRes=NULL){
  n_list = length(stableSplits$stableTrees)
  for (i in n_list:1) {
    res = stableTrees(stableSplits$stableTrees[1:i],simRes)
    if(length(res$ids)>0){
      return(list(minStableNodes= unlist(stableSplits$stableNodes[1:i]), minStableTrees = res))
    }
  }
}

# Find top n unstable splits
nTopUnstableSplits <- function(node_list, criteria=0.8, n=1){

  n_sim = sum(node_list[[1]][[1]])
  n_node = length(names(node_list))
  top_n_v = character(n)
  j = 1
  for (i in 1:n_node) {
    node_list[[i]][[1]]/n_sim->pmf
    max_id = which.max(pmf[1,])
    if(pmf[1,max_id]<criteria){
      if(j<=n){
        top_n_v[j] = names(node_list[i])
      }
      j=j+1
    }
  }
  return(top_n_v)
}

# stb_list <- stabilityAnalyze(simRes)
nextSamplingRegion <- function(nid, X_range, stb_list, simRes, first_class_criterion){

  stableSplits(stb_list$node_list,first_class_criterion)->ss
  minStableTrees(ss,simRes)->min_stbts
  snipNodes(min_stbts$minStableNodes)->snip_nodes_ids
  minStableTree <- snip.rpart(min_stbts$minStableTrees$trees[[1]], toss = snip_nodes_ids)

  nsr = samplingRegion(nid, X_range, treeInfo(minStableTree))
  return(nsr)

}
