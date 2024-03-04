library(rpart)
library(kddt)
library(pdp)     # for partial dependence plots
library(rpart.plot)
library(randomForest)
#===============================================================================
# Helper function
#===============================================================================
# change max and min values for unifying plot color scale
newMaxMin2data <- function(max, min, data){
  if(dim(data)[2]!=3){
    return("Please input right data!")
  }
  if(data[which.max(data[,3]),3] < max){
    data[which.max(data[,3]),3] <- max
  }
  if(data[which.min(data[,3]),3] > min){
    data[which.min(data[,3]),3] <- min
  }
  return(data)
}

#===============================================================================
# Load simulated data
#===============================================================================
# setwd(getwd())  set the correct directory
load(file="./others/data_sim.RData")

head(data_sim)
dim(data_sim)
data.max = max(data_sim$yhat)
data.min = min(data_sim$yhat)
# show the 2D function
rwb <- colorRampPalette(c("blue", "white", "red"))
plotPartial(data_sim, contour = TRUE, col.regions = rwb, contour.color='darkgray')
# decision tree rpart general settings
my.ctrl <- rpart.control(xval=10, minbucket=2, minsplit=4,  cp=0.001)
my.toss <- c(8,9,10,11,24,50,51,13,14,15)

#===============================================================================
# The optimal partition fitted using decision tree by using the true function
# (entire dataset)
#===============================================================================
# A decision tree fitted using the entire dataset
boston_tree <- rpart(yhat~.,data = data_sim, control = my.ctrl)
# rpart.plot(boston_tree)
boston_tree <- snip.rpart(boston_tree, toss=my.toss)
rpart.plot(boston_tree)
true.pd <- partial(boston_tree, pred.var = c("x1", "x2"))
scaled.true.pd = newMaxMin2data(data.max, data.min, true.pd)
plotPartial(scaled.true.pd, contour = T, col.regions = rwb, contour.color='darkgray')
#===============================================================================
# simulating data from the true function (entire dataset)
#===============================================================================
set.seed(6995)
n = length(data_sim$yhat)
sample(1:n,50)->ids
df.sample = data_sim[ids,]  # sample 50 observations
scaled.df.sample = newMaxMin2data(data.max, data.min, df.sample)
# show the samples
plotPartial(scaled.df.sample, contour = TRUE, col.regions = rwb)
#===============================================================================
# fit a ordinary decision tree with the 50 samples
#===============================================================================
set.seed(1)
odt.sim = rpart(yhat~.,data = df.sample, control = my.ctrl)
odt.sim <- snip.rpart(odt.sim, toss=my.toss)
rpart.plot(odt.sim,roundint=FALSE)
odt.sim.pd <- partial(odt.sim, pred.var = c("x1", "x2"))
scaled.odt.sim.pd = newMaxMin2data(data.max, data.min, odt.sim.pd)
plotPartial(scaled.odt.sim.pd, contour = T, col.regions = rwb, contour.color='darkgray')

#===============================================================================
# fit a random forest with the 50 samples
#===============================================================================
set.seed(1)
rf.sim = randomForest(yhat~.,data = df.sample)
pd.rf.sim = partial(rf.sim,pred.var = c("x1","x2"))
# show the prediction of RF
scaled.pd.rf.sim = newMaxMin2data(data.max, data.min, pd.rf.sim)
plotPartial(scaled.pd.rf.sim, contour = TRUE, col.regions = rwb)

#===============================================================================
# fit a distillation decision tree from the RF
#===============================================================================
fitedModel = rf.sim
predMethod='predict'

X <- df.sample[,-3] # the covariates
y <- df.sample[,3]

dim(df.sample)[1] -> nobs
nSim=30
samplingParameters <- setSamplingParameters(n.interpretable=500*nobs, n.predictive=10000*nobs,
                                            exploreNodes=list(c(12,1000*nobs)))

rf.kddt.interpretive <- inductionInterpretableTree(c(8,10,50,14), "regression", X, fitedModel, predMethod,
                                                     samplingParameters, nSim)
rpart.plot(rf.kddt.interpretive,roundint=FALSE)

# fit the full KDDT
sim.rf.kddt <- kddtInductionFromInterpretable(X, y, rf.kddt.interpretive, fitedModel, predMethod, samplingParameters)

kddtPredict(sim.rf.kddt, newX = data_sim[,c(1,2)], predict_type='vector')->sim.rf.kddt.pred
df.kddt.pred<-data.frame(cbind(data_sim[,c(1,2)],sim.rf.kddt.pred))
colnames(df.kddt.pred)<-c('x1','x2','yhat')
class(df.kddt.pred)<-c("data.frame","partial")
scaled.df.kddt.pred = newMaxMin2data(data.max, data.min, df.kddt.pred)
plotPartial(scaled.df.kddt.pred, contour = T, col.regions = rwb, contour.color='darkgray')

# split stability
for (i in 1:9) {
  print(rf.kddt.interpretive$stability[[i]]$first)
  plotSecondLevelStability(nid=rf.kddt.interpretive$stability[[i]]$nid, rf.kddt.interpretive$stability[[i]]$second)
}


#===============================================================================
# show the interpretable partition
#===============================================================================
my.toss <- c(8,9,10,11,12,13,14,15)
# The best partition
boston_tree <- snip.rpart(boston_tree, toss=my.toss)
predict(boston_tree,newdata = data_sim[,c(1,2)])->best.pred
df.best.pred<-data.frame(cbind(data_sim[,c(1,2)],best.pred))
colnames(df.best.pred)<-c('x1','x2','yhat')
class(df.best.pred)<-c("data.frame","partial")
scaled.df.best.pred = newMaxMin2data(data.max, data.min, df.best.pred)
plotPartial(scaled.df.best.pred, contour = F, col.regions = rwb, contour.color='darkgray')

# The odt partition
odt.sim <- snip.rpart(odt.sim, toss=my.toss)
predict(odt.sim,newdata = data_sim[,c(1,2)])->odt.pred
df.odt.pred<-data.frame(cbind(data_sim[,c(1,2)],odt.pred))
colnames(df.odt.pred)<-c('x1','x2','yhat')
class(df.odt.pred)<-c("data.frame","partial")
scaled.df.odt.pred = newMaxMin2data(data.max, data.min, df.odt.pred)
plotPartial(scaled.df.odt.pred, contour = F, col.regions = rwb, contour.color='darkgray')

# the kddt partition
rf.kddt.interpretive <- inductionInterpretableTree(c(8,10,12,14), "regression", X, fitedModel, predMethod,
                                                  samplingParameters, nSim)
kddtPredict(rf.kddt.interpretive,newX = data_sim[,c(1,2)], predict_type='vector')->kddt.predictive.pred
kddt.predictive.pred<-data.frame(cbind(data_sim[,c(1,2)],kddt.predictive.pred))
colnames(kddt.predictive.pred)<-c('x1','x2','yhat')
class(kddt.predictive.pred)<-c("data.frame","partial")
scaled.kddt.predictive.pred = newMaxMin2data(data.max, data.min, kddt.predictive.pred)
plotPartial(scaled.kddt.predictive.pred, contour = F, col.regions = rwb, contour.color='darkgray')

df.comp.kddt<-df.best.pred
df.comp.kddt$yhat<-abs(df.comp.kddt$yhat-kddt.predictive.pred$yhat)

df.comp.odt<-df.best.pred
df.comp.odt$yhat<-abs(df.comp.odt$yhat-df.odt.pred$yhat)
# make same color scale to compare
df.comp.odt$yhat[which.max(df.comp.odt$yhat)] <- df.comp.kddt$yhat[which.max(df.comp.kddt$yhat)]

rwb.new <- colorRampPalette(c("white", "darkgreen"))
plotPartial(df.comp.odt, contour = F, col.regions = rwb.new, contour.color='darkgray')
plotPartial(df.comp.kddt, contour = F, col.regions = rwb.new, contour.color='darkgray')

sum(df.comp.odt$yhat^2)/length(df.comp.odt$yhat)
sum(df.comp.kddt$yhat^2)/length(df.comp.kddt$yhat)

