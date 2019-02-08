setwd("xxxxxxxxxxxx")

library("ZIBR", lib.loc="~/R/win-library/3.5")
library("dplyr", lib.loc="~/R/win-library/3.5")
library("nlme", lib.loc="C:/Program Files/R/R-3.5.1/library")



############ Control vs LOS vs NEC #############
#### Filter OTUs with low presence and low abundance
taxa.data1=read.delim("taxa.data1.csv", sep=',', row.names=1)
taxa.data1 = t(taxa.data1/27994)

sparse0.05_taxa.data1=which(rowMeans(taxa.data1[,1:68])<0.0005 & rowMeans(taxa.data1[,69:80])<0.0005 & rowMeans(taxa.data1[,81:96])<0.0005)
taxa.data1=taxa.data1[-sparse0.05_taxa.data1,]
sparse_p50_taxa.data1=which(rowSums(taxa.data1[,1:68]>0)<34 & rowSums(taxa.data1[,69:80]>0)<6 & rowSums(taxa.data1[,81:96]>0)<8)
taxa.data1=taxa.data1[-sparse_p50_taxa.data1,]

taxa.data1 <- t(taxa.data1)
write.table(taxa.data1,file="filtered_matrix1.csv",row.names = TRUE,col.names = TRUE,sep= ",")


#### create covariates
cov1=read.delim("cov1.csv", sep=',',stringsAsFactors = FALSE)
X <- data.frame(cov1[,c('Time','Treat')])
rownames(X) <- cov1$Sample
Z <- X
subject.ind <- cov1$Subject
time.ind    <- as.numeric(cov1$Time)

#### Fit ZIBR model for NEC vs LOS vs Control
spe.all <- colnames(taxa.data1)
p.species.list.zibr <- list()

for (spe in spe.all){
  ### Print process on screen
  cat(spe,'\n')
  Y <- taxa.data1[, spe]
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  
  ### run ZIBR
    est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
                subject.ind=subject.ind,
                time.ind=time.ind,
                quad.n=30,verbose=TRUE)
    p.species.list.zibr[[spe]] <- est$joint.p
  #break
}


#### Fit LME model for NEC vs LOS vs Control
p.species.list.lme <- list()

for (spe in spe.all){
  Y <- taxa.data1[, spe]
  tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=subject.ind)
  lme.fit <- lme(Y.tran ~ Time + Treat,random=~1| SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
  p.species.list.lme[[spe]] <- coef.mat[,2]
  #break
}


#### export p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
write.csv(p.species.zibr,file=paste('Results_ZIBR_NECvsLOSvsControl.csv',sep=''))

p.species.lme <- t(as.data.frame(p.species.list.lme))
write.csv(p.species.lme,file=paste('Results_LME_NECvsLOSvsControl.csv',sep=''))






############ Control vs LOS #############
#### Filter OTUs with low presence and low abundance
taxa.data1=read.delim("taxa.data1_LOS.csv", sep=',', row.names=1)
taxa.data1 = t(taxa.data1/27994)

sparse0.05_taxa.data1=which(rowMeans(taxa.data1[,1:68])<0.0005 & rowMeans(taxa.data1[,69:80])<0.0005)
taxa.data1=taxa.data1[-sparse0.05_taxa.data1,]
sparse_p50_taxa.data1=which(rowSums(taxa.data1[,1:68]>0)<34 & rowSums(taxa.data1[,69:80]>0)<6)
taxa.data1=taxa.data1[-sparse_p50_taxa.data1,]

taxa.data1 <- t(taxa.data1)
write.table(taxa.data1,file="filtered_matrix1_LOS.csv",row.names = TRUE,col.names = TRUE,sep= ",")


#### create covariates
cov1=read.delim("cov1_LOS.csv", sep=',',stringsAsFactors = FALSE)
X <- data.frame(cov1[,c('Time','Treat')])
rownames(X) <- cov1$Sample
Z <- X
subject.ind <- cov1$Subject
time.ind    <- as.numeric(cov1$Time)

#### Fit ZIBR model for LOS vs Control
spe.all <- colnames(taxa.data1)
p.species.list.zibr <- list()

for (spe in spe.all){
  ### Print process on screen
  cat(spe,'\n')
  Y <- taxa.data1[, spe]
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  
  ### run ZIBR
  est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
              subject.ind=subject.ind,
              time.ind=time.ind,
              quad.n=30,verbose=TRUE)
  p.species.list.zibr[[spe]] <- est$joint.p
  #break
}


#### Fit LME model for LOS vs Control
p.species.list.lme <- list()

for (spe in spe.all){
  Y <- taxa.data1[, spe]
  tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=subject.ind)
  lme.fit <- lme(Y.tran ~ Time + Treat,random=~1| SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
  p.species.list.lme[[spe]] <- coef.mat[,2]
  #break
}


#### export p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
write.csv(p.species.zibr,file=paste('Results_ZIBR_ControlvsLOS.csv',sep=''))

p.species.lme <- t(as.data.frame(p.species.list.lme))
write.csv(p.species.lme,file=paste('Results_LME_ControlvsLOS.csv',sep=''))




############ Control vs NEC #############
#### Filter OTUs with low presence and low abundance
taxa.data1=read.delim("taxa.data1_NEC.csv", sep=',', row.names=1)
taxa.data1 = t(taxa.data1/27994)

sparse0.05_taxa.data1=which(rowMeans(taxa.data1[,1:68])<0.0005 & rowMeans(taxa.data1[,69:84])<0.0005)
taxa.data1=taxa.data1[-sparse0.05_taxa.data1,]
sparse_p50_taxa.data1=which(rowSums(taxa.data1[,1:68]>0)<34 & rowSums(taxa.data1[,69:84]>0)<8)
taxa.data1=taxa.data1[-sparse_p50_taxa.data1,]

taxa.data1 <- t(taxa.data1)
write.table(taxa.data1,file="filtered_matrix1_NEC.csv",row.names = TRUE,col.names = TRUE,sep= ",")


#### create covariates
cov1=read.delim("cov1_NEC.csv", sep=',',stringsAsFactors = FALSE)
X <- data.frame(cov1[,c('Time','Treat')])
rownames(X) <- cov1$Sample
Z <- X
subject.ind <- cov1$Subject
time.ind    <- as.numeric(cov1$Time)

#### Fit ZIBR model for NEC vs Control
spe.all <- colnames(taxa.data1)
p.species.list.zibr <- list()

for (spe in spe.all){
  ### Print process on screen
  cat(spe,'\n')
  Y <- taxa.data1[, spe]
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  
  ### run ZIBR
  est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
              subject.ind=subject.ind,
              time.ind=time.ind,
              quad.n=30,verbose=TRUE)
  p.species.list.zibr[[spe]] <- est$joint.p
  #break
}


#### Fit LME model for NEC vs Control
p.species.list.lme <- list()

for (spe in spe.all){
  Y <- taxa.data1[, spe]
  tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=subject.ind)
  lme.fit <- lme(Y.tran ~ Time + Treat,random=~1| SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
  p.species.list.lme[[spe]] <- coef.mat[,2]
  #break
}


#### export p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
write.csv(p.species.zibr,file=paste('Results_ZIBR_ControlvsNEC.csv',sep=''))

p.species.lme <- t(as.data.frame(p.species.list.lme))
write.csv(p.species.lme,file=paste('Results_LME_ControlvsNEC.csv',sep=''))






############ LOS vs NEC #############
#### Filter OTUs with low presence and low abundance
taxa.data1=read.delim("taxa.data2.csv", sep=',', row.names=1)
taxa.data1 = t(taxa.data1/27994)

sparse0.05_taxa.data1=which(rowMeans(taxa.data1[,1:18])<0.0005 & rowMeans(taxa.data1[,19:36])<0.0005)
taxa.data1=taxa.data1[-sparse0.05_taxa.data1,]
sparse_p50_taxa.data1=which(rowSums(taxa.data1[,1:18]>0)<9 & rowSums(taxa.data1[,19:36]>0)<9)
taxa.data1=taxa.data1[-sparse_p50_taxa.data1,]

taxa.data1 <- t(taxa.data1)
write.table(taxa.data1,file="filtered_matrix2.csv",row.names = TRUE,col.names = TRUE,sep= ",")


#### create covariates
cov1=read.delim("cov2.csv", sep=',',stringsAsFactors = FALSE)
X <- data.frame(cov1[,c('Time','Treat')])
rownames(X) <- cov1$Sample
Z <- X
subject.ind <- cov1$Subject
time.ind    <- as.numeric(cov1$Time)

#### Fit ZIBR model for NEC vs LOS
spe.all <- colnames(taxa.data1)
p.species.list.zibr <- list()

for (spe in spe.all){
  ### Print process on screen
  cat(spe,'\n')
  Y <- taxa.data1[, spe]
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  
  ### run ZIBR
  est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
              subject.ind=subject.ind,
              time.ind=time.ind,
              quad.n=30,verbose=TRUE)
  p.species.list.zibr[[spe]] <- est$joint.p
  #break
}


#### Fit LME model for NEC vs LOS
p.species.list.lme <- list()

for (spe in spe.all){
  Y <- taxa.data1[, spe]
  tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=subject.ind)
  lme.fit <- lme(Y.tran ~ Time + Treat,random=~1| SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
  p.species.list.lme[[spe]] <- coef.mat[,2]
  #break
}


#### export p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
write.csv(p.species.zibr,file=paste('Results_ZIBR_LOSvsNEC.csv',sep=''))

p.species.lme <- t(as.data.frame(p.species.list.lme))
write.csv(p.species.lme,file=paste('Results_LME_LOSvsNEC.csv',sep=''))




