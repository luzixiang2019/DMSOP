#model
#\vz_t=(\mI_n \otimes \mGamma)\vz_(t-1)+\mX_t\vbeta+(\mDelta \otimes \mI_S)\vdelta+\vxi
#\vz_0=\mX_0\vbeta_0+(\mDelta \otimes \mI_S)\vdelta+\vxi
#\vdelta=\rho*(W \otimes \mI_S)\vdelta+\epsilon
#\epsilon \sim \rN (0, \mI_{m} \otimes \diag(\vsigma^2))
#\vxi \sim \rN (0, \mI_{n} \otimes \mV)
#y_{ikst}=l_s if \gamma_{s,l_s-1} < z_{ikst} \le \gamma_{s,l_s}


#packages
####################################
#packages
library(msm)
library(mvtnorm) 
library(coda)
library(Matrix)
library(spatialprobit)

#simulated data 
######################################################################################

# set weight matrix
####################################

#set.seed(1)
library(readxl)
X4_empirical_data <- read_excel("6 empirical data_use.xlsx", 
                                sheet = "Weight", col_names = FALSE)

W <-as.matrix(X4_empirical_data)
lam<-eigen(W)
lam<-lam$values
ilammin<-1/min(lam)
ilammax<-1/max(lam)

W==t(W)

save(W, file = "Weight-original.Rdata")

W=W/rowSums(W)
lam<-eigen(W)
lam<-lam$values
ilammin<-1/min(lam)
ilammax<-1/max(lam)
save(W, file = "Weight.Rdata")
rm(list = ls())
####################################



#empirical data
########################################################################
#packages
library(msm)
library(mvtnorm) 
library(coda)
library(Matrix)
library(spatialprobit)
library(readxl)
####################################

# load Weight matrix
  load("Weight.Rdata")
  library(readxl)
  data<- read_excel("6 empirical data_use.xlsx", sheet = "data-three-use")
  data<-as.matrix(data)
  S<-2
  TM<-2
  qs<-6
  m<-length(unique(data[,2]))
  n<-dim(data)[1]
  prov<-unique(data[,2])
  
  nob<-c()
  for (i in 1:m) 
  {
    nob[i] <- length(which(data[,2]==prov[i]))
  }
  
  n=sum(nob)
  nS=n*S
  mS=m*S
  data<-as.matrix(data,n,)
  yt01<-as.matrix(data[,4])
  yt02<-as.matrix(data[,5])
  yt0<-matrix(NA,nS,1)
  yt0[c(T,F)]<-yt01
  yt0[c(F,T)]<-yt02
  
  plus<-qs+S+3
  
  yt11<-as.matrix(data[,(4+plus)])
  yt12<-as.matrix(data[,(5+plus)])
  yt1<-matrix(NA,nS,1)
  yt1[c(T,F)]<-yt11
  yt1[c(F,T)]<-yt12
  
  yt21<-as.matrix(data[,(4+plus*2)])
  yt22<-as.matrix(data[,(5+plus*2)])
  yt2<-matrix(NA,nS,1)
  yt2[c(T,F)]<-yt21
  yt2[c(F,T)]<-yt22
  
  y<-as.matrix(c(yt1,yt2))
  
  q<-qs*S
  
  xt0<-matrix(0,nS,q)
  xt01<-data[,6:(6+qs-1)]
  xt02<-data[,6:(6+qs-1)]
  xt0[c(T,F), 1:qs]<-xt01
  xt0[c(F,T),-(1:qs)]<-xt02
  
  
  xt1<-matrix(0,nS,q)
  xt11<-data[,(6+plus):(6+plus+qs-1)]
  xt12<-data[,(6+plus):(6+plus+qs-1)]
  xt1[c(T,F), 1:qs]<-xt11
  xt1[c(F,T),-(1:qs)]<-xt12
  
  xt2<-matrix(0,nS,q)
  xt21<-data[,(6+plus*2):(6+plus*2+qs-1)]
  xt22<-data[,(6+plus*2):(6+plus*2+qs-1)]
  xt2[c(T,F), 1:qs]<-xt21
  xt2[c(F,T),-(1:qs)]<-xt22
  
  
  x<-rbind(xt1,xt2)
  
  tim=matrix(c(rep(1,nS),rep(2,nS)),nS*TM,1)
  
  tst<-matrix(1,nob[1]*S,1)
  for (i in 2:m) {
    tst<-rbind(tst,matrix(i,nob[i]*S,1))
  }
  region<-matrix(rep(tst,TM),nS*TM,1)
  
  
  tst1=matrix(NA,nS,1)
  indi<-matrix(seq(1,nob[1]),nob[1],1)
  for (i in 2:m) {
    indi<-rbind(indi,matrix(seq(1,nob[i]),nob[i],1))
  }
  tst1[c(T,F)]=indi
  tst1[c(F,T)]=indi
  individual<-matrix(rep(tst1,TM),nS*TM,1)
  
  #dat=data.frame(cbind(tst,tst1,y,mX))
  dat=cbind(tim,region,individual,y,x)
  
  save(W,dat,yt0,xt0,file = "All data_use.Rdata")
  
  province<-table(data[,2])
  summaryresults=function(x){
    c(mean=mean(x),sd=sd(x),quantile(x,c(0,0.5,1)),effSize=effectiveSize(x))
  }
  
  Age_12<-summaryresults(xt01[,1])
  Income_12<-summaryresults(xt01[,6])
  gender_12<-table(xt01[,2])
  urban_12<-table(xt01[,3])
  education_12<-summaryresults(xt01[,4])
  employ_12<-table(xt01[,5])
  Age_12
  Income_12
  gender_12
  urban_12
  education_12
  employ_12
  
  Age_14<-summaryresults(xt11[,1])
  Income_14<-summaryresults(xt11[,6])
  gender_14<-table(xt11[,2])
  urban_14<-table(xt11[,3])
  education_14<-summaryresults(xt11[,4])
  employ_14<-table(xt11[,5])
  
  
  
  Age_16<-summaryresults(xt21[,1])
  Income_16<-summaryresults(xt21[,6])
  gender_16<-table(xt21[,2])
  urban_16<-table(xt21[,3])
  education_16<-summaryresults(xt21[,4])
  employ_16<-table(xt21[,5])
  
  LF_12<-table(yt01[,1])
  Health_12<-table(yt02[,1])
  LF_14<-table(yt11[,1])
  Health_14<-table(yt12[,1])
  LF_16<-table(yt21[,1])
  Health_16<-table(yt22[,1])
  
  LF_12
  LF_14
  LF_16
  
  Health_12
  Health_14
  Health_16

  
  Age_12
  Age_14
  Age_16
  
  Income_12
  Income_14
  Income_16
  
  gender_12
  gender_14
  gender_16
  
  urban_12
  urban_14
  urban_16
  
  education_12
  education_14
  education_16
  
  employ_12
  employ_14
  employ_16


