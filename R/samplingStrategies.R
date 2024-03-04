#' Set sampling parameters
#' @param n.interpretable the number of pseudo samples for fit the interpretable tree.
#' @param n.predictive the number of pseudo samples for fit the predictive subtrees.
#' @param samplingMethod a string denoting the sampling method, c("uniform","empirical").
#' @param intCovariates a vector of string includes the name of covariates whose type are integer.
#' @param exploreNodes a list includes the nodes need to explore more, each element of the list includes a node id and the sample size.
#' @param seed a integer denotes the number of random sampling
#'
#' @return a list of sampling parameters
#' @export
#'
setSamplingParameters <- function(n.interpretable, n.predictive, samplingMethod="uniform",
                                  intCovariates=NULL, exploreNodes=NULL, seed=123){
  if(!(samplingMethod %in% c("uniform","empirical"))){
    return("Please use the option 'uniform' !")
  }
  if(samplingMethod=='uniform'){
    res = list(samplingMethod='randomSampling', sampleSize=c(n.interpretable,n.predictive),
               exploreNodes=exploreNodes, seed = seed)
  }else{
    res = list(samplingMethod='empiricalSampling', sampleSize=c(n.interpretable,n.predictive),
               exploreNodes=exploreNodes, seed = seed)
  }
  return(res)
}

#
setSamplingStrategy <- function(samplingRegion, samplingParameters, is.interpretable, X=NULL, PCApercentageCutoff=0.8){
  if(is.null(samplingParameters)){
    return("Please set sampling parameters!")
  }
  sampleSize = 0
  if(is.interpretable){
    sampleSize = samplingParameters$sampleSize[1]
  }else{
    sampleSize = samplingParameters$sampleSize[2]
  }

  samplingMethod = samplingParameters$samplingMethod

  seed = samplingParameters$seed

  intCovariates = samplingParameters$intCovariates

  if(samplingMethod %in% c("randomSampling","empiricalSampling")){
    res = list(samplingMethod = samplingMethod, samplingRegion = samplingRegion,
               sampleSize = sampleSize, intCovariates=intCovariates, seed = seed)
  }else{
    res = list(samplingMethod = samplingMethod, samplingRegion = samplingRegion,
               sampleSize = sampleSize, percentageVariance = PCApercentageCutoff, X = X,
               intCovariates=intCovariates, seed = seed)
  }
  return(res)
}

#===============================================================================
# uniform random sampling
randomSampling <- function(samplingStrategy, dX){

  n = samplingStrategy$sampleSize
  samplingRegion = samplingStrategy$samplingRegion

  X_sampled <- data.frame(matrix(NA,    # Create empty data frame
                                 nrow = n,
                                 ncol = length(samplingRegion$data_range$names)))
  colnames(X_sampled) <- samplingRegion$data_range$names


  n.factor = length(samplingRegion$id_range$factor)
  n.numeric = dim(samplingRegion$id_range$numeric)[2]
  set.seed(samplingStrategy$seed)
  if(n.factor>0){
    for (i in 1:n.factor) {
      sample(samplingRegion$id_range$factor[[i]], n, replace = TRUE) ->
        X_sampled[,which(colnames(X_sampled) == names(samplingRegion$id_range$factor)[i])]
    }
  }
  if(!is.null(n.numeric)){
    for (i in 1:n.numeric) {
      if(colnames(samplingRegion$id_range$numeric)[i] %in% samplingStrategy$intCovariates){
        round(runif(n, min = floor(samplingRegion$id_range$numeric[1,i]),
                    max = ceiling(samplingRegion$id_range$numeric[2,i]))) ->
          X_sampled[,which(colnames(X_sampled) == colnames(samplingRegion$id_range$numeric)[i])]
      }else{
        runif(n, min = samplingRegion$id_range$numeric[1,i],max = samplingRegion$id_range$numeric[2,i]) ->
          X_sampled[,which(colnames(X_sampled) == colnames(samplingRegion$id_range$numeric)[i])]
      }
    }
  }
  # factor()
  for (i in 1:length(samplingRegion$data_range$names)) {
    if(colnames(X_sampled)[i] %in% names(samplingRegion$data_range$factor)){
      X_sampled[,i] = factor(X_sampled[,i],levels = samplingRegion$data_range$factor[[colnames(X_sampled)[i]]])
    }
  }
  return(X_sampled)
}
#===============================================================================
# empirical sampling
empiricalSampling <- function(samplingStrategy, dX){

  n = samplingStrategy$sampleSize
  samplingRegion = samplingStrategy$samplingRegion

  X_sampled <- data.frame(matrix(NA,    # Create empty data frame
                                 nrow = n,
                                 ncol = length(samplingRegion$data_range$names)))
  colnames(X_sampled) <- samplingRegion$data_range$names


  n.factor = length(samplingRegion$id_range$factor)
  n.numeric = dim(samplingRegion$id_range$numeric)[2]
  set.seed(samplingStrategy$seed)
  if(n.factor>0){
    for (i in 1:n.factor) {
      if(length(samplingRegion$id_range$factor[[i]])>1){
        names(samplingRegion$id_range$factor)[i] -> name_i
        dX[,c(name_i)] -> x_i
        rows_include = (x_i %in% samplingRegion$id_range$factor[[i]])
        sprob = as.vector(prop.table(table(dX[rows_include,c(name_i)])))
        sample(samplingRegion$id_range$factor[[i]], n, replace = TRUE, prob = sprob) ->
          X_sampled[,which(colnames(X_sampled) == names(samplingRegion$id_range$factor)[i])]
      }else{
        sample(samplingRegion$id_range$factor[[i]], n, replace = TRUE) ->
          X_sampled[,which(colnames(X_sampled) == names(samplingRegion$id_range$factor)[i])]
      }
    }
  }
  if(!is.null(n.numeric)){
    for (i in 1:n.numeric) {
      colnames(samplingRegion$id_range$numeric)[i] -> name_i
      dX[,c(name_i)] -> x_i
      new_x_i = x_i[x_i>=samplingRegion$id_range$numeric[1,i] & x_i<=samplingRegion$id_range$numeric[2,i]]
      pct = length(new_x_i)/length(x_i)
      if(pct<=0.3){
        runif(n, min = samplingRegion$id_range$numeric[1,i],max = samplingRegion$id_range$numeric[2,i]) ->
          X_sampled[,which(colnames(X_sampled) == colnames(samplingRegion$id_range$numeric)[i])]
      }else{
        my_r = new_r(new_x_i,type = "discrete")
        my_r(n) -> X_sampled[,which(colnames(X_sampled) == colnames(samplingRegion$id_range$numeric)[i])]
      }
    }
  }
  # factor()
  for (i in 1:length(samplingRegion$data_range$names)) {
    if(colnames(X_sampled)[i] %in% names(samplingRegion$data_range$factor)){
      X_sampled[,i] = factor(X_sampled[,i],levels = samplingRegion$data_range$factor[[colnames(X_sampled)[i]]])
    }
  }
  return(X_sampled)
}



#===============================================================================
# PCA sampling, please ignore
#===============================================================================
pcaCutOff <- function(res.pca, percentageVariance, p){
  # ratio of sample size
  for (i in 1:p) {
    if((percentageVariance-sum((res.pca$sdev^2/sum(res.pca$sdev^2))[1:i]))<=0){
      cfInd = i
      break
    }
  }
  ceiling(res.pca$sdev^2/(res.pca$sdev^2)[cfInd]) -> size_v
  if(0 %in% size_v){
    return("In X, observations' number n is less than covariates number p.")
  }else{
    return(list(cfInd = cfInd, size_v = size_v))
  }
}

# samplingRegionNumeric: a matrix for the hypercube support of continuous covariates. id_range
pcaSamplingContinuous <- function(res.pca, samplingRegionNumeric, pcaRegion, pcaCf, p, seed) {

  size_v = pcaCf$size_v
  cfInd = pcaCf$cfInd

  set.seed(seed)
  # sampling in the pca region
  temp_s = list()
  for (i in 1:p) {
    temp_s[[i]] <- runif(size_v[i],min = pcaRegion[1,i], pcaRegion[2,i])
  }
  #
  if(cfInd>2){
    for (i in 2:(cfInd-1)) {
      temp_s[[i]] <- sample(temp_s[[i]],size_v[1],replace = TRUE)
    }
  }
  for (i in cfInd:p) {
    temp_s[[i]] <- rep(temp_s[[i]],size_v[1])
  }

  # construct the matrix
  pcaSample = matrix(nrow = size_v[1],ncol = 0)
  for (i in 1:p) {
    pcaSample = cbind(pcaSample,temp_s[[i]])
  }

  # pca filtering
  X_temp = pcaSample%*%t(res.pca$rotation)
  X_sampled = t(t(X_temp) * res.pca$scale + res.pca$center)

  good_flag = rep(TRUE,size_v[1])
  for (j in 1:p) {
    good_flag = good_flag & (X_sampled[, j] >= samplingRegionNumeric[1,j] &
                               X_sampled[, j] <= samplingRegionNumeric[2,j])
  }
  X_sampled = X_sampled[good_flag,]

  if(sum(good_flag)==1){
    X_sampled = matrix(X_sampled,nrow = 1, ncol = length(X_sampled))
  }
  return(as.data.frame(X_sampled))

}

# pca sampling include both categorical and numerical covariates
pcaSamplingSizeOne <- function(res.pca, samplingRegion, pcaRegion, pcaCf, p, seed){

  # do pca sampling
  pca_sampled = pcaSamplingContinuous(res.pca, samplingRegion$id_range$numeric, pcaRegion, pcaCf, p, seed)
  n_good_numeric_sampled = dim(pca_sampled)[1]

  set.seed(seed)
  # if pcaSamplingContinuous' return is not null
  if(n_good_numeric_sampled>0){
    n.factor = length(samplingRegion$data_range$factor)
    # if there are some categorical covariates, draw them and combined them with numeric covariates
    if(n.factor>0){
      for (i in 1:n.factor) {
        pca_sampled <- cbind(pca_sampled,
                             rep(sample(samplingRegion$id_range$factor[[i]], 1, replace = TRUE),
                                 n_good_numeric_sampled))
      }
    }
    colnames(pca_sampled)<-c(colnames(samplingRegion$data_range$numeric),
                             names(samplingRegion$data_range$factor))
    return(pca_sampled)
  }else{
    # if pcaSamplingContinuous' return is null
    return(NULL)
  }
}

# marginal random sampling based on pca random sampling
marginalPCArandomSampling <- function(samplingStrategy){

  samplingRegion = samplingStrategy$samplingRegion

  # first random sampling
  if(is.null(samplingStrategy$X)){
    X = randomSampling(samplingStrategy)
  }else{
    X = samplingParameters$X
  }

  n.numeric = dim(samplingRegion$data_range$numeric)[2]
  # if numerical covariates don't exist
  if(is.null(n.numeric)){
    # if there are only categorical covariates, do random sampling only in categorical ones.
    return(randomSampling(samplingStrategy))
  }else{
    # if there are continuous covariates, first, do pca analysis.
    # find continuous covariates
    X.continuous = X[,colnames(samplingRegion$data_range$numeric)]
    res.pca <- prcomp(X.continuous,scale=TRUE)
    pcaRegion = (dataRange(as.data.frame(res.pca$x)))$numeric
    p.con = dim(X.continuous)[2]
    pcaCf = pcaCutOff(res.pca, samplingStrategy$percentageVariance, p.con)
    #print(p.con)
    # draw PCA samples
    pca_sampled = NULL
    sz = 0
    while (sz<n) {

      pca_sampled = rbind(pca_sampled, pcaSamplingSizeOne(res.pca, samplingRegion, pcaRegion,
                                                          pcaCf, p.con, samplingStrategy$seed))
      sz = dim(pca_sampled)[1]
      if(is.null(sz)){
        sz = 0
      }
    }
    # re-factorize all factors
    for (i in 1:length(samplingRegion$data_range$names)) {
      if(colnames(pca_sampled)[i] %in% names(samplingRegion$data_range$factor)){
        pca_sampled[,i] = factor(pca_sampled[,i],
                                 levels = samplingRegion$data_range$factor[[colnames(pca_sampled)[i]]])
      }
    }
    return(pca_sampled)
  }
}

