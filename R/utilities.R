#=========================================================================================================
# utilities functions
#=========================================================================================================

parseCondition <- function(condition){

  signs1 = c(">=","<=","!=")
  signs2 = c("<","=",">")

  for (j_1 in 1:length(signs1)) {
    if(grepl(signs1[j_1],condition,fixed = T)){
      strsplit(condition,signs1[j_1])[[1]] -> temp_v
      split=trimws(temp_v[2])
      return(c(split,signs1[j_1]))
    }
  }
  for (j_2 in 1:length(signs2)) {
    if(grepl(signs2[j_2],condition,fixed = T)){
      strsplit(condition,signs2[j_2])[[1]] -> temp_v
      split=trimws(temp_v[2])
      return(c(split,signs2[j_2]))
    }
  }
}

isNumeric <- function(strName, dR){

  if(strName %in% dR$names){
    if(strName %in% names(dR$factor)){
      return(F)
    }else{
      return(T)
    }
  }else{
    return(NULL)
  }
}

#===============================================================================
get_ancestors_and_siblings <- function(node_id) {
  # initialize empty list to store ancestor and sibling ids
  ids <- list()

  # add current node id as first element
  ids[[1]] <- node_id

  # loop through ancestors and siblings until root node is reached
  while (node_id != 1) {
    # get parent node id
    parent_id <- floor(node_id/2)

    # get sibling node ids
    if (node_id %% 2 == 0) {
      sibling_id <- node_id + 1
    } else {
      sibling_id <- node_id - 1
    }

    # add parent and sibling node ids to list
    ids <- c(list(parent_id), list(sibling_id), ids)

    # update node id to parent node id for next iteration
    node_id <- parent_id
  }

  # return list of ancestor and sibling ids
  return(unlist(ids))
}

#===============================================================================
depth_first <- function(node_ids) {
  # Initialize an empty result vector and a stack for DFS
  res <- c()
  stack <- c()

  # Start with the root node
  stack <- c(stack, node_ids[1])

  # Traverse the tree in depth-first order
  while (length(stack) > 0) {
    # Pop the next node from the stack
    curr_node <- stack[length(stack)]
    stack <- stack[-length(stack)]

    # Add the current node to the result vector
    res <- c(res, curr_node)

    # Find the left and right child nodes (if they exist)
    left_child <- 2 * curr_node
    right_child <- 2 * curr_node + 1

    # Push the right child node onto the stack (if it exists)
    if (right_child %in% node_ids) {
      stack <- c(stack, right_child)
    }

    # Push the left child node onto the stack (if it exists)
    if (left_child %in% node_ids) {
      stack <- c(stack, left_child)
    }
  }

  return(res)
}

# node_ids <- get_ancestors_and_siblings(38)
# depth_first(node_ids)

#===============================================================================
width_first <- function(node_ids) {
  # Initialize queue with root node
  queue <- list(node_ids[1])
  visited <- list()
  i <- 2

  # Perform level-order traversal
  while(length(queue) > 0) {
    # Dequeue next node
    current_node <- queue[[1]]
    queue <- queue[-1]

    # Add current node to visited list
    visited <- c(visited, current_node)

    # Enqueue children of current node
    left_child <- node_ids[i]
    if(!is.na(left_child)) {
      queue <- c(queue, left_child)
    }
    i <- i + 1

    right_child <- node_ids[i]
    if(!is.na(right_child)) {
      queue <- c(queue, right_child)
    }
    i <- i + 1
  }

  return(unlist(visited))
}

# treeIds(c(5, 11, 38))->a
# width_first(sort(a))
#===============================================================================

treeIds <- function(r_ids,res=NULL){
  r_ids <- sort(r_ids,decreasing=T)
  get_ancestors_and_siblings(r_ids[1])->e_ids
  r_ids[!r_ids %in% e_ids]->r_ids
  # print(e_ids)
  # print(r_ids)
  res=c(res,e_ids)
  if(length(r_ids)>0){
    treeIds(r_ids,res=res)
  }else{
    return(unique(res))
  }
}

# treeIds(c(5, 11, 38))


#===============================================================================
remove_bottom_nodes <- function(node_ids) {
  # Define function to check if a node is a leaf
  is_leaf <- function(node_id) {
    # If the node has no children, it is a leaf
    left_child <- 2 * node_id
    right_child <- 2 * node_id + 1
    !(left_child %in% node_ids || right_child %in% node_ids)
  }

  # Identify leaf nodes
  leaf_nodes <- node_ids[sapply(node_ids, is_leaf)]

  # Remove leaf nodes
  pruned_nodes <- setdiff(node_ids, leaf_nodes)

  return(pruned_nodes)
}

# remove_bottom_nodes(depth_first(node_ids))

#===============================================================================
# sample() function  can instead of this function.
# C is a vector
runifCategorical <- function(n, C){

  r_v = runif(n,min = 0,max = length(C))
  res = character(n)
  for (j in 1:n) {
    res[j] = C[ceiling(r_v[j])]
  }
  return(res)
}

