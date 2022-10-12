library(evaluomeR)
library(RSKC)
library(sparcl)


evaluomeRSupportedCBI()
# TODO, to add this in stabilityRange, and quality/qualityRange
dataFrame <- stability(data=ontMetrics, cbi="rskc", k=3, all_metrics=TRUE, bs=100)
assay(dataFrame)

dataFrame <- stability(data=ontMetrics, cbi="kmeans", k=5, bs=100)


sim <-
  function(mu,f){
    D<-matrix(rnorm(60*f),60,f)
    D[1:20,1:50]<-D[1:20,1:50]+mu
    D[21:40,1:50]<-D[21:40,1:50]-mu
    return(D)
  }
sim
d0<-sim(1,500)# generate a dataset
true<-rep(1:3,each=20) # vector of true cluster labels
d<-d0
ncl<-3
for ( i in 1 : 10){
  d[sample(1:60,1),sample(1:500,1)]<-rnorm(1,mean=0,sd=15)
}

# The generated dataset looks like this...
pairs(
  d[,c(1,2,3,200)],col=true,
  labels=c("clustering feature 1",
           "clustering feature 2","clustering feature 3",
           "noise feature1"),
  main="The sampling distribution of 60 cases colored by true cluster labels",
  lower.panel=NULL)

d

# RSKC works when more than 2 columns are provided

r3<-RSKC(d[,1:5],ncl,alpha=10/60,L1=6,nstart=200)
