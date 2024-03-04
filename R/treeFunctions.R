
#===============================================================================
# tree functions
#===============================================================================

#' The induction of DDT
#' @param method one of "anova", "poisson", "class" or "exp". If method is missing then the routine tries to make an intelligent guess. If y is a survival object, then method = "exp" is assumed, if y has 2 columns then method = "poisson" is assumed, if y is a factor then method = "class" is assumed, otherwise method = "anova" is assumed. It is wisest to specify the method directly, especially as more criteria may added to the function in future. Alternatively, method can be a list of functions named init, split and eval. Examples are given in the file ‘tests/usersplits.R’ in the sources, and in the vignettes ‘User Written Split Functions’.
#' @param control a list of options that control details of the rpart algorithm. See rpart.control.
#' @param predict_type character string, c("vector", "prob", "class", "matrix"), denoting the type of predicted value returned. If the rpart object is a classification tree, then the default is to return prob predictions, a matrix whose columns are the probability of the first, second, etc. class. (This agrees with the default behavior of tree). Otherwise, a vector result is returned.
#'
#' @return a list includes the parameters of rpart function.
#' @export
#'
setRpartPara <- function(method,control,predict_type){
  return(list(method = method, control = control, predict_type = predict_type))
}

# Get all ancestors of a node by its id
getRids <- function(id){
  rids=id
  temp = id
  while (temp!=1) {
    rids=c(rids,floor(temp/2))
    temp=floor(temp/2)
  }
  rids
}

# Get all node's ids of a complete tree by its depth
getAllNidByLevel <- function(level){
  return(1:(2^(level+1)-1))
}


rangeConditions <- function(id, treeInfo) {
  df_numeric <- data.frame(
    nid = integer(), var = character(), sign = character(), split = double(),
    stringsAsFactors = FALSE
  )
  df_factor <- data.frame(
    nid = integer(), var = character(), sign = character(), split = character(),
    stringsAsFactors = FALSE
  )
  treeInfo[treeInfo$nid %in% getRids(id), c("nid", "conditions")] -> conditions
  signs1 <- c(">=", "<=", "!=")
  signs2 <- c("<", "=", ">")
  for (i in 1:dim(conditions)[1]) {
    flag_signs1 <- FALSE
    for (j_1 in 1:length(signs1)) {
      if (grepl(signs1[j_1], conditions[i, "conditions"], fixed = T)) {
        strsplit(conditions[i, "conditions"], signs1[j_1])[[1]] -> temp_v
        if (signs1[j_1] != "!=") {
          df_numeric <- rbind(
            df_numeric,
            data.frame(
              nid = as.integer(conditions[i, "nid"]),
              var = trimws(temp_v[1]), sign = signs1[j_1],
              split = as.numeric(trimws(temp_v[2]))
            )
          )
        } else {
          df_factor <- rbind(
            df_factor,
            data.frame(
              nid = as.integer(conditions[i, "nid"]),
              var = trimws(temp_v[1]), sign = signs1[j_1],
              split = trimws(temp_v[2])
            )
          )
        }
        flag_signs1 <- TRUE
        break
      }
    }
    if (!flag_signs1) {
      for (j_2 in 1:length(signs2)) {
        if (grepl(signs2[j_2], conditions[i, "conditions"], fixed = T)) {
          strsplit(conditions[i, "conditions"], signs2[j_2])[[1]] -> temp_v
          if (signs2[j_2] != "=") {
            df_numeric <- rbind(
              df_numeric,
              data.frame(
                nid = as.integer(conditions[i, "nid"]),
                var = trimws(temp_v[1]), sign = signs2[j_2],
                split = as.numeric(trimws(temp_v[2]))
              )
            )
          } else {
            df_factor <- rbind(
              df_factor,
              data.frame(
                nid = as.integer(conditions[i, "nid"]),
                var = trimws(temp_v[1]), sign = signs2[j_2],
                split = trimws(temp_v[2])
              )
            )
          }
          break
        }
      }
    }
  }
  return(list(numeric = df_numeric, factor = df_factor))
}


treeInfo <- function(tree, digits = 5, minlength = 0L){
  frame=tree$frame
  frame$nid = as.numeric(row.names(tree$frame))
  frame$conditions<-labels(tree, digits = digits, minlength = minlength)
  frame
}


# Get the split information.
splitTable <- function(treeInfo){
  treeInfo[treeInfo$var!="<leaf>",c("nid","var")] -> table
  table$split = NA
  dim(table)[1]->n
  for (i in 1:n) {
    tempId = 2*table$nid[i]
    rangeConditions(tempId,treeInfo) -> conditions
    a = conditions$numeric[which(conditions$numeric$nid==tempId),"split"]
    b = conditions$factor[which(conditions$factor$nid==tempId),"split"]
    if(length(a>0)){
      table$split[i]=a
    }else if(length(b>0)){
      table$split[i]=b
    }
  }
  return(table[order(table$nid),])
}


# given ddt pseudo tree information and one observation of X return ddt leaf node id
getDdtLeafNid <- function(x, pseudoTreeInfo, dataRange){

  id = 1 # start from the root

  while (id %in% pseudoTreeInfo$nid[which(pseudoTreeInfo$var!="<leaf>")]) { # if the id is the internal node

    splitVarName = pseudoTreeInfo$var[which(pseudoTreeInfo$nid==id)]
    ps = which(colnames(x)==splitVarName) # find the position of split variable
    parseCondition(pseudoTreeInfo$conditions[which(pseudoTreeInfo$nid==2*id)])->aaa # find the split value
    split_value = aaa[1]
    sg = aaa[2]

    if(isNumeric(splitVarName,dataRange)){ # if split variable is continuous
      if(x[,ps] < as.numeric(split_value)){ # put into left child
        if(sg %in% c("<=", "<")){
          id = 2*id
        }else{
          id = 2*id+1
        }
      }else{ # put into right child
        if(sg %in% c("<=", "<")){
          id = 2*id+1
        }else{
          id = 2*id
        }
      }
    }else{
      # if split variable is categorical
      if(as.character(x[,ps]) %in% unlist(strsplit(split_value, ","))){
        id = 2*id
      }else{ # put into right child
        id = 2*id+1
      }
    }
  }
  return(id)
}


snipNodes <- function(stableNodes){
  n_list = length(stableNodes)
  stableNodes_v = as.integer(unlist(stableNodes))
  for (i in n_list:1) {
    id_v = getRids(as.integer(stableNodes[[i]]))
    for (id in id_v) {
      if(!(id %in% stableNodes_v)){
        stableNodes_v[i]=NA
      }
    }
  }
  as.vector(na.omit(stableNodes_v))->nodes_keep
  l_children = nodes_keep*2
  r_children = nodes_keep*2 + 1
  return(setdiff(union(l_children, r_children),nodes_keep))

}

#=====================================================
# version 0.1.0
#=====================================================


# Insert the child tree's rpart.object$frame into its parent tree's rpart.object$frame.
insertChildFrame <- function(tree, nid, stump){
  pFrame = tree$frame
  cFrame = stump$frame
  np = dim(pFrame)[1]
  # the node ids of child tree
  old_cids = as.integer(row.names(cFrame))
  # convert the child node ids to the parent node ids
  new_cids = old_cids+(nid-1)*2^(floor(log2(old_cids)))
  row.names(cFrame)<-as.character(new_cids)
  # find the row number with nid, insert cFrame after this row.
  row_num <- which(row.names(pFrame) == as.character(nid))
  if(row_num<np){
    tree$frame <- rbind(pFrame[1:(row_num-1),], cFrame, pFrame[(row_num+1):np,])
  }else{
    tree$frame <- rbind(pFrame[1:(row_num-1),], cFrame)
  }
  tree
}


insertChildSplits <- function(tree, nid, stump){
  pSplits = tree$splits
  max.pos = max(pSplits[pSplits[,2]>1,4])
  cSplits = stump$splits
  cSplits[cSplits[,2]>1,4]<-cSplits[cSplits[,2]>1,4]+max.pos
  treeFrame = tree$frame
  nonLeafTreeFrame = treeFrame[treeFrame$var!='<leaf>',]
  insertId = as.integer(row.names(nonLeafTreeFrame)[which(row.names(nonLeafTreeFrame)==as.character(nid))-1])
  # find the position of last insertId in pSplits
  row_num = max(which(pSplits[,"nids"]==insertId))
  np = nrow(pSplits)

  if(row_num<np){
    tree$splits <- rbind(pSplits[1:row_num,], cSplits, pSplits[(row_num+1):np,])
  }else{
    tree$splits <- rbind(pSplits[1:row_num,], cSplits)
  }
  tree

}


insertChildCsplits <- function(tree, nid, stump){
  pCsplits = tree$csplit
  if(is.null(pCsplits)){
    tree$csplit = stump$csplit
  }else{
    cCsplits = stump$csplit
    treeFrame = tree$frame
    nonLeafTreeFrame = treeFrame[treeFrame$var!='<leaf>',]
    insertId = as.integer(row.names(nonLeafTreeFrame)[which(row.names(nonLeafTreeFrame)==as.character(nid))-1])
    # find the position of last insertId in pCsplits
    row_num = max(which(as.integer(rownames(pCsplits))==insertId))
    np = nrow(pCsplits)

    if(row_num<np){
      tree$csplit <- rbind(pCsplits[1:row_num,], cCsplits, pCsplits[(row_num+1):np,])
    }else{
      tree$csplit <- rbind(pCsplits[1:row_num,], cCsplits)
    }
  }
  tree
}


# create a stump (node) from the results of Algorithm 1.
createStump <- function(stumpsRes, X_range){

  # a summary of the simulation results, a string matrix includes the split variable and split values.
  simResMatStump = stumpsToMat(stumpsRes$stump_list)
  # the number of simulation must greater than 2.
  if(ncol(simResMatStump)<2){
    return(NULL)
  }
  # first-level stability
  pmf_mat = firstLevelStabilityStump(X_range$names, simResMatStump)
  # select the splitting variable according first-level stability
  bestCovName = colnames(pmf_mat)[which.max(pmf_mat)]
  #print(pmf_mat)
  # second-level stability
  secStb = secondLevelStabilityStump(nameCov=bestCovName,
                                     isNumeric=isNumeric(bestCovName, X_range),
                                     simResMatStump)
  # update the stump with the stable split information
  stump = stumpsRes$stump_list[[secStb$stumpId]]
  return(list(stump=stump, firstLevelStability=pmf_mat, secondLevelStability=secStb))
}

# update the information in the rpart.object$frame, n, wt
getUpdateInfo <- function(nid, stump, tree, X){
  # the stump
  stump$variable.importance -> imp.stump
  splitCovNames = names(imp.stump)
  imp.percent = imp.stump/sum(imp.stump)
  stump.imp.df = data.frame(splitCovNames,imp.percent)
  treeInfo(stump)->tInfo
  splitTable(tInfo)->spInfo
  covname = spInfo[1,"var"]
  X_range = dataRange(X)
  isnumeric = covname %in% colnames(X_range$numeric)
  isfactor = covname %in% names(X_range$factor)
  if(!(isnumeric | isfactor)){
    return("no such split variable!")
  }
  sv = spInfo[1,"split"]
  rangeConditions(2,tInfo)->conditions
  # the data
  n = nrow(X)
  if(nid==1){
    where = rep(1,n)
  }else{
    where = tree$where
  }
  #print(tree$where)
  ids = c(1:n)
  whereMat = cbind(ids, where)
  # find the subset of data according the parent nid.
  subWhereMat = whereMat[(whereMat[,2]==nid),]
  length(subWhereMat)/2 -> snr # the number of rows of subset
  n_left = 0
  n_right = 0
  imp_list = NULL  # for calculating variable importance
  if(snr>1){
    # assign data to corresponding node.
    for (i in 1:snr) {
      if(isnumeric){
        conditions$numeric$sign[1]->sg
        if(sg %in% c("<=", "<")){
          if(X[subWhereMat[i,1],covname] < sv){
            whereMat[subWhereMat[i,1],2] = 2*nid
            subWhereMat[i,2] = 2*nid
            n_left=n_left+1
          }else{
            whereMat[subWhereMat[i,1],2] =  2*nid+1
            subWhereMat[i,2] = 2*nid+1
            n_right=n_right+1
          }
        }else{
          if(X[subWhereMat[i,1],covname] >= sv){
            whereMat[subWhereMat[i,1],2] = 2*nid
            subWhereMat[i,2] = 2*nid
            n_left=n_left+1
          }else{
            whereMat[subWhereMat[i,1],2] =  2*nid+1
            subWhereMat[i,2] = 2*nid+1
            n_right=n_right+1
          }
        }
      }else{
        #
        if(X[subWhereMat[i,1],covname] %in% unlist(strsplit(sv, ","))){
          whereMat[subWhereMat[i,1],2] = 2*nid
          subWhereMat[i,2] = 2*nid
          n_left=n_left+1
        }else{
          whereMat[subWhereMat[i,1],2] =  2*nid+1
          subWhereMat[i,2] = 2*nid+1
          n_right=n_right+1
        }
      }
    }
    imp_list = list(nid=nid, split_name = covname, sub_ids_nids = subWhereMat, imp.df=stump.imp.df)
  }
  list(sNmae = covname, n=snr, nl=n_left, nr=n_right, where = whereMat[,2], impList = imp_list)
}

# insert the new stump into the tree
insertStump <- function(nid, createStumpRes, tree, X){
  stump <- createStumpRes$stump
  # get the information for updating
  updateInfo <- getUpdateInfo(nid, stump, tree, X)
  stump$frame[1,c("n","wt")] = updateInfo$n
  stump$frame[2,c("n","wt")] = updateInfo$nl
  stump$frame[3,c("n","wt")] = updateInfo$nr
  stump$splits[,"count"] = updateInfo$n
  nids = rep(nid,nrow(stump$splits))
  stump$splits = cbind(stump$splits,nids)
  if(!is.null(stump$csplit)){
    rownames(stump$csplit)=rep(nid,nrow(stump$csplit))
  }
  if(nid==1){ # the root split
    stump$where = updateInfo$where
    stump$dsplits = c(nid,updateInfo$where)
    stump$splitNames = updateInfo$sNmae
    stump$impList = list(updateInfo$impList)
    stump$stability = list(list(nid=1, first=createStumpRes$firstLevelStability,
                           second=createStumpRes$secondLevelStability))
    return(stump)
  }
  tree$where = updateInfo$where
  tree$dsplits = rbind(tree$dsplits, c(nid,updateInfo$where))
  tree$splitNames = c(tree$splitNames, updateInfo$sNmae)
  if(!is.null(updateInfo$impList)){
    tree$impList = append(tree$impList,list(updateInfo$impList))
  }
  tree = insertChildFrame(tree,nid,stump)
  tree = insertChildSplits(tree, nid, stump)
  if(!is.null(stump$csplit)){
    tree = insertChildCsplits(tree, nid, stump)
  }
  tree$stability = append(tree$stability, list(list(nid=nid,first=createStumpRes$firstLevelStability,
                         second=createStumpRes$secondLevelStability)))
  return(tree)
}
























