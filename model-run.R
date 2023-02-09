#model
#\vz_t=(\mI_n \otimes \mGamma)\vz_(t-1)+\mX_t\vbeta+(\mDelta \otimes \mI_S)\vdelta+\vxi
#\vz_0=\mX_0\vbeta_0+(\mDelta \otimes \mI_S)\vdelta+\vxi
#\vdelta=\rho*(W \otimes \mI_S)\vdelta+\epsilon
#\epsilon \sim \rN (0, \mI_{m} \otimes \diag(\vsigma^2))
#\vxi \sim \rN (0, \mI_{n} \otimes \mV)
#y_{ikst}=l_s if \gamma_{s,l_s-1} < z_{ikst} \le \gamma_{s,l_s}

#packages
####################################
library(msm)
library(mvtnorm) 
library(coda)
library(Matrix)
library(spatialprobit)
library(spatialreg)
#library(stats)
####################################

####################################
#load data
#dat<- read.csv(file=file.path(("."), "test.csv"), header=TRUE, sep=",")
#W<- read.csv(file=file.path(("."), "testFW.csv"), header=FALSE, sep=",")

load("All data_use.Rdata")


####################################

####################################
# Variable assignment
data=na.omit(dat)
S<-2

time<-as.matrix(data[,1])
TM<-length(unique(time))

y<-as.matrix(data[,4])
y1<-as.matrix(y[c(T,F)])
y2<-as.matrix(y[c(F,T)])
yt01<-as.matrix(yt0[c(T,F)])
yt02<-as.matrix(yt0[c(F,T)])
n<-dim(y)[1]/S/TM
nS<-n*S
for (tm in 1:TM) {
  assign(paste("yt",tm,sep = ""),as.matrix(y[((tm-1)*nS+1):(tm*nS),]))
}


x<-as.matrix(data[,-c(1,2,3,4)])
for (tm in 1:TM) {
  assign(paste("xt",tm,sep = ""),as.matrix(x[((tm-1)*nS+1):(tm*nS),]))
}
q<-dim(x)[2]
qs<-q/S
x1<-x[c(T,F), 1:qs]
x2<-x[c(F,T),-(1:qs)]
xt01<-as.matrix(xt0[c(T,F),])
xt02<-as.matrix(xt0[c(F,T),])
W<-as.matrix(W)
#W=matrix(c(0,1,1,1,0,0,1,0,0),3,3,byrow=TRUE)
#W=W/rowSums(W)


n1<-dim(x)[1]/S/TM
n2<-dim(W)[1]
n3<-dim(W)[2]

#obtain the number of region and observations in region i, and level of y
index<-data[,2:3]
index<-as.matrix(index)
m<-index[which.max(index[,1])]
L1<-length(unique(y1))
L2<-length(unique(y2))
nob<-c()
for (i in 1:m) 
{
  nob[i] <- length(which(index[,1]==i))/S/TM
}

#obtain the Delta
one<-matrix()
Delta<-matrix()
for (i in 1:m){
  one<-matrix(1,nob[i],1)
  Delta<-as.matrix(bdiag(Delta,one))
}
Delta<-Delta[-1,-1]
####################################

####################################
# error check

#check size of observation
# check for more than one region, else no check
if (m>1) {
  # check for region 1:m-1
  for (ck in 1:(m-1)) #index[sm,2]) is not an number, all ways a space + number
  {
    sm<-sum(nob[1:ck])
    if(!isTRUE(all.equal(as.integer(index[sm*S,2]),nob[ck]))||!isTRUE(all.equal(as.integer(index[(sm+1)*S,2]),1))){
      stop('wrong size of observation')}
  }
  #check for region m
  sm<-sum(nob[1:m])
  if(!isTRUE(all.equal(as.integer(index[sm*S,2]),nob[m]))){
    stop('wrong size of observation')}
}
# check the whole number of obeservation in different regions
if(!isTRUE(all.equal(sum(nob),n))){
  stop('wrong size of observation')
}

#check matrix X
if(!isTRUE(all.equal(n1,n))){
  stop('wrong size of X matrix')
}

if(!isTRUE(all.equal(sum(x1[,1]),n*TM))) 
{
  tst<-apply(x1, 2, sum)
  ind<-which(tst==n*TM)
  if(length(ind)>0) 
  {
    stop('intercept term must be in first column')
  } 
  else
  { if (length(ind)==0) {
    cflag1<-0 # we do not have intercept term
  }
  } 
} 

if(isTRUE(all.equal(sum(x1[,1]),n*TM))) {
  cflag1<-1 # we have intercept term
} 
if(!isTRUE(all.equal(sum(x2[,1]),n*TM))) 
{
  tst<-apply(x2, 2, sum)
  ind<-which(tst==n)
  if(length(ind)>0) 
  {
    stop('intercept term must be in first column')
  } 
  else
  { if (length(ind)==0) {
    cflag2<-0 # we do not have intercept term
  }
  } 
} 
if(isTRUE(all.equal(sum(x2[,1]),n*TM))) {
  cflag2<-1 # we have intercept term
} 

#check matrix W
if(!isTRUE(all.equal(n2,m))){
  stop('wrong size of W matrix')
}
if(!isTRUE(all.equal(n2,n3))){
  stop('wrong size of W matrix')
}
####################################

####################################
#priors
vc<-matrix(0,q,1)
mT<-diag(q)*(1e+12) 
e<-0
g<-0
lam<-eigen(W)
lam<-lam$values
ilammin<-1/min(lam)
ilammax<-1/max(lam)
ifelse(ilammin<(-1), ilammin<-(-1), ilammin<-ilammin) 
ifelse(ilammax>1, ilammax<-1, ilammax<-ilammax)
kappa10<-rep(0,L1-3)
K10<-100*diag(L1-3)
kappa20<-rep(0,L2-3)
K20<-100*diag(L2-3)
n0<-1
Q0<-matrix(c(0.01,0.01,0.01,0.01),S,S)
d0<-matrix(0,S,1)
D0<-diag(S)*(1e+12) 
#prior<-list(mbeta=vc,vbeta=mT,rval=iota,sigmsp=e,sigmscl=g,rmin=ilammin,rmax=ilammax)
####################################


##############################
#starting values
mS=m*S
nST=nS*TM
oneT<-matrix(1,TM,1)
DeltaIS<-kronecker(Delta,diag(S))
#DeltaIS<-Matrix(kronecker(Delta,diag(S)),sparse=TRUE)
OneTDeltaIS<-kronecker(oneT,DeltaIS)

#OneTDeltaIS<-Matrix(kronecker(oneT,DeltaIS),sparse=TRUE)
WIS<-kronecker(W,diag(S))

beta0<-matrix(1,q,1)
beta<-matrix(1,q,1)
delta<-matrix(0.1,mS,1)
vsigm2<-c(0.01,0.01)
rho<-0.5
a<-0.2 #for MH step of rho
#starting value for V
mV=matrix(c(1,0.1,0.1,1),S,S)
kappa1<-c(runif(L1-3))
kappa2<-c(runif(L2-3))
gamm1 <- c(-Inf, seq(0,1,length.out = L1-1), +Inf)
gamm2 <- c(-Inf, seq(0,1,length.out = L2-1), +Inf)
z<-matrix(0,nST,1)
zt0<-matrix(0,nS,1)
Lam<-matrix(c(0,0,0,0),S,S)
##############################


#Funcions
##############################
#log FCD of rho# 
logFCDrho=function(rho,delta,vsigm2,m,S,W){
  Brho<-diag(m*S)-rho*kronecker(W,diag(S))
  tst<-t(delta)%*%t(Brho)%*%solve(kronecker(diag(m),diag(vsigm2)))%*%Brho%*%delta
  tst1<-as.double(log(det(Brho))-tst/2)
  return(tst1)
}
##############################

#---------------------------------------------------

# PURPOSE: compute the log determinant |I_n - rho*W|
# using the user-selected (or default) method
# ---------------------------------------------------
# USAGE: detval = sar_lndet(lflag,W,rmin,rmax)
# where ldetflag,rmin,rmax,W contains input flags 
# ---------------------------------------------------
# Code is now based on method spatialreg::do_ldet written by Roger Bivand
# ldetflag = 0 : Pace and Barry (1997) grid
# Pace, R. and Barry, R. (1997). Quick computation of spatial autoregressive estimators.
# Geographics Analysis, 29(3):232-247.
#
# ldetflag = 1 : Pace and LeSage (2004) Chebyshev approximation
# Pace, R. and LeSage, J. (2004). Chebyshev approximation of log-determinants of
# spatial weight matrices. Computational Statistics & Data Analysis, 45(2):179-
# 196.
#
# ldetflag = 2 : Barry and Pace (1999) MC approximation
# Barry, R. and Pace, R. (1999). Monte Carlo estimates of the log determinant of
# large sparse matrices. Linear Algebra and its Applications, 289(1-3):41-54.
sar_lndet <- function(ldetflag,W,rmin,rmax){
  results <- NULL
  env <- new.env(parent=globalenv())
  listw <- mat2listw(W)
  assign("n", nrow(W), envir=env)
  assign("listw", listw, envir=env)
  assign("family", "SAR", envir=env)
  assign("similar", FALSE, envir=env)
  
  # do lndet approximation calculations if needed
  if( ldetflag == 0 ){ # no approximation
    t0            <- Sys.time()
    SE_classic_setup(env)
    # fine grid is already computed
    results$detval  <- get("detval1", env)
  }else if( ldetflag == 1 ) { 
    t0   <- Sys.time()
    cheb_setup(env, q=4)  # Chebychev approximation, q = 4
    # For Chebyshev-Approximation we must compute the fine grid with do_ldet
    # to be later used for numerical integration
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else if( ldetflag == 2 ) {
    t0   <- Sys.time()
    mcdet_setup(env, p=16, m=30, which=1)
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else{
    # TODO: Pace and Barry, 1998 spline interpolation
    stop('sar_lndet: method not implemented')
  }
  results$time   <- Sys.time() - t0
  return( results )
}
# ---------------------------------------------------


# ---------------------------------------------------
# detval = an ngrid x 2 matrix with rho-values and lndet values
#detval<-sar_lndet(0,W,rmin,rmax)$detval# no approximation
#detval<-sar_lndet(0,W,-0.9,0.9)$detval# no approximation
#detval<-sar_lndet(1,W,-1,1)$detval#Chebyshev approximation
#detval<-sar_lndet(1,W,-0.9,0.9)$detval#Chebyshev approximation
#detval<-sar_lndet(2,W,-1,1)$detval #MC approximation
# ---------------------------------------------------

# ---------------------------------------------------
draw_rho_uni=function(detval,epe0,eped,epe0d,n,k,rho,a1,a2){
  #      epe0 <- t(e0) %*% e0
  #      eped <- t(ed) %*% ed
  #      epe0d<- t(ed) %*% e0 
  #PURPOSE: draws rho-values using griddy Gibbs and inversion
  # ---------------------------------------------------
  #  USAGE: rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)
  # where info contains the structure variable with inputs 
  # and the outputs are either user-inputs or default values
  # ---------------------------------------------------
  # REFERENCES: LeSage and Pace (2009) Introduction to Spatial Econometrics
  # Chapter 5, pp 136-141 on Bayesian spatial regression models.
  
  # where:       detval = an ngrid x 2 matrix with rho-values and lndet values
  #                  e0 = y - x*b0;
  #                 ed = Wy - x*bd;
  #               epe0 = e0'*e0
  #               eped = ed'*ed
  #              epe0d = ed'*e0
  #               nobs = # of observations
  #               nvar = # of explanatory variables
  #            logdetx = log(det(x'*x))
  #                 a1 = parameter for beta prior on rho
  #                 a2 = parameter for beta prior on rho
  # detval1 umbenannt in => rho_grid = grid/vector of rhos
  # detval2 umbenannt in => lndet = vector of log-determinants
  # detval1sq = rvec^2
  # lnbprior = vector of log prior densities for rho (default: Beta(1,1))
  #
  lnbprior <- log(dbeta(detval[,1],a1,a2)) #less then 0, the value is -Inf
  #u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval)  # do_ldet() liefert nur 2000 statt 2001 Gridpoints
  nmk      <- (n-k)/2
  rho_grid   <- detval[,1]  # rho grid values
  lndet    <- detval[,2]  # log-determinant grid values
  rho_gridsq <- rho_grid * rho_grid
  #rho_gridsq=rho_grid.*rho_grid
  
  iota <- matrix(1,nrho,1)
  z <- epe0*iota - 2*rho_grid*epe0d + rho_gridsq*eped
  z <- -nmk*log(z)
  den <- lndet + z + lnbprior         # vector of log posterior densities for rho vector
  # posterior density post(rho | data) \propto likelihood(data|rho,beta,z) * prior(rho)
  # log posterior density post(rho | data) \propto loglik(data|rho,beta,z) + log prior(rho)
  n <- nrho
  adj <- max(den)
  den <- den - adj                     # adjustieren/normieren der log density auf maximum 0; 
  x <- exp(den)                        # density von p(rho) --> pdf
  
  # trapezoid rule
  yy <- (rho_grid[2:nrho] + rho_grid[1:(nrho-1)])
  isum <- sum(yy * (x[2:nrho] - x[1:(nrho - 1)])/2)
  z <- abs(x/isum)
  den <- cumsum(z)
  u<-runif(1)
  rnd <- u * sum(z)
  ind <- which(den <= rnd)
  idraw <- max(ind)
  if (idraw > 0 && idraw < nrho) {
    results <- rho_grid[idraw]
  }else {
    results <- rho
  }
  return(results)
}

# ---------------------------------------------------
#muti
draw_rho_multi=function(detval,epe01,eped1,epe0d1,epe02,eped2,epe0d2,n,k,rho,a1,a2){
  #      epe0 <- t(e0) %*% e0
  #      eped <- t(ed) %*% ed
  #      epe0d<- t(ed) %*% e0 
  #PURPOSE: draws rho-values using griddy Gibbs and inversion
  # ---------------------------------------------------
  #  USAGE: rho = draw_rho(detval,epe0,eped,epe0d,n,k,rho,a1,a2)
  # where info contains the structure variable with inputs 
  # and the outputs are either user-inputs or default values
  # ---------------------------------------------------
  # REFERENCES: LeSage and Pace (2009) Introduction to Spatial Econometrics
  # Chapter 5, pp 136-141 on Bayesian spatial regression models.
  
  
  # where:       detval = an ngrid x 2 matrix with rho-values and lndet values
  #                  e0 = y - x*b0;
  #                 ed = Wy - x*bd;
  #               epe0 = e0'*e0
  #               eped = ed'*ed
  #              epe0d = ed'*e0
  #               nobs = # of observations
  #               nvar = # of explanatory variables
  #            logdetx = log(det(x'*x))
  #                 a1 = parameter for beta prior on rho
  #                 a2 = parameter for beta prior on rho
  # detval1 umbenannt in => rho_grid = grid/vector of rhos
  # detval2 umbenannt in => lndet = vector of log-determinants
  # detval1sq = rvec^2
  # lnbprior = vector of log prior densities for rho (default: Beta(1,1))
  #
  lnbprior <- log(dbeta(detval[,1],a1,a2)) #less then 0, the value is -Inf
  #u        <- runif(thinning * ndraw + burn.in)   # u ~ U(0, 1)
  nrho     <- nrow(detval)  # do_ldet() liefert nur 2000 statt 2001 Gridpoints
  nmk      <- (n-k)/2
  rho_grid   <- detval[,1]  # rho grid values
  lndet    <- detval[,2]  # log-determinant grid values
  rho_gridsq <- rho_grid * rho_grid
  #rho_gridsq=rho_grid.*rho_grid
  
  iota <- matrix(1,nrho,1)
  z1 <- epe01*iota - 2*rho_grid*epe0d1 + rho_gridsq*eped1
  z1 <- -nmk*log(z1)
  z2 <- epe02*iota - 2*rho_grid*epe0d2 + rho_gridsq*eped2
  z2 <- -nmk*log(z2)
  z=z1+z2
  den <- 2*lndet + z + lnbprior         # vector of log posterior densities for rho vector
  # posterior density post(rho | data) \propto likelihood(data|rho,beta,z) * prior(rho)
  # log posterior density post(rho | data) \propto loglik(data|rho,beta,z) + log prior(rho)
  n <- nrho
  adj <- max(den)
  den <- den - adj                     # adjustieren/normieren der log density auf maximum 0; 
  x <- exp(den)                        # density von p(rho) --> pdf
  
  # trapezoid rule
  yy <- (rho_grid[2:nrho] + rho_grid[1:(nrho-1)])
  isum <- sum(yy * (x[2:nrho] - x[1:(nrho - 1)])/2)
  z <- abs(x/isum)
  den <- cumsum(z)
  u<-runif(1)
  rnd <- u * sum(z)
  ind <- which(den <= rnd)
  idraw <- max(ind)
  if (idraw > 0 && idraw < nrho) {
    results <- rho_grid[idraw]
  }else {
    results <- rho
  }
  return(results)
}


##############################
#conditional mean of zt0#

CMZt0=function(zt0,xt0){
  ht0<-xt0%*%beta0+kronecker(Delta,diag(S))%*%delta
  
  hatht0=matrix(NA,nS,1)
  
  for (i in 1:nS) {
    if (i%%2==1) {
      s=1
      hatht0[i]<-ht0[i]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(zt0[i+1]-ht0[i+1])
    }
    if (i%%2==0) {
      s=2
      hatht0[i]<-ht0[i]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(zt0[i-1]-ht0[i-1])
    }
  }
  return(hatht0)
}
##############################

##############################
#conditional mean of z#

CMZ=function(z,zj1,Lam,x){
  
  InLam<-kronecker(diag(n),Lam)
  ITInLam<-kronecker(diag(TM),InLam)
  
  h<-ITInLam%*%zj1+x%*%beta+OneTDeltaIS%*%delta
  
  hath=matrix(NA,nST,1)
  
  for (i in 1:nST) {
    if (i%%2==1) {
      s=1
      hath[i]<-h[i]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[i+1]-h[i+1])
    }
    if (i%%2==0) {
      s=2
      hath[i]<-h[i]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[i-1]-h[i-1])
    }
  }
  return(hath)
}
##############################


##############################
#-log FCD of kappa1#
FCDkappa1=function(kappa1){
  gamm1[1]<--Inf
  gamm1[2]<-0
  for (l1 in 3:(L1-1)) {
    gamm1[l1]<-(gamm1[l1-1]+exp(kappa1[l1-2]))/(1+exp(kappa1[l1-2]))
  }
  
  SUM<-0
  hatht01<-as.matrix(hatht0[c(T,F),])
  hath1<-as.matrix(hath[c(T,F),])
  hatv1<-mV[1,1]-mV[1,-1]%*%solve(mV[-1,-1])%*%mV[-1,1]
  for (l1 in 2:(L1-1)) {
    logphi<-log(pnorm((gamm1[l1+1]-hath1[y1==l1])/c(sqrt(hatv1)))-
                  pnorm((gamm1[l1]-hath1[y1==l1])/c(sqrt(hatv1))))
    SUM<-SUM+as.double(sum(logphi))
  }
  
  for (l1 in 2:(L1-1)) {
    logphit0<-log(pnorm((gamm1[l1+1]-hatht01[yt01==l1])/c(sqrt(hatv1)))-
                    pnorm((gamm1[l1]-hatht01[yt01==l1])/c(sqrt(hatv1))))
    SUM<-SUM+as.double(sum(logphit0))
  }
  
  SUMk<-0
  for (l1 in 2:(L1-2)) {
    logJc<-log(1-gamm1[l1])+kappa1[l1-1]-2*log(1+exp(kappa1[l1-1]))
    SUMk<-SUMk+logJc
  }
  
  logP1=SUM+SUMk+(-0.5*t(kappa1-kappa10)%*%solve(K10)%*%(kappa1-kappa10))
  -logP1
}
##############################



##############################
#-log FCD of kappa2#
FCDkappa2=function(kappa2){
  gamm2[1]<--Inf
  gamm2[2]<-0
  for (l2 in 3:(L2-1)) {
    gamm2[l2]<-(gamm2[l2-1]+exp(kappa2[l2-2]))/(1+exp(kappa2[l2-2]))
  }
  
  SUM<-0
  hatht02<-as.matrix(hatht0[c(F,T),])
  hath2<-as.matrix(hath[c(F,T),])
  hatv2<-mV[2,2]-mV[2,-2]%*%solve(mV[-2,-2])%*%mV[-2,2]
  for (l2 in 2:(L2-1)) {
    logphi<-log(pnorm((gamm2[l2+1]-hath2[y2==l2])/c(sqrt(hatv2)))-
                  pnorm((gamm2[l2]-hath2[y2==l2])/c(sqrt(hatv2))))
    SUM<-SUM+as.double(sum(logphi))
  }
  
  for (l2 in 2:(L2-1)) {
    logphit0<-log(pnorm((gamm2[l2+1]-hatht02[yt02==l2])/c(sqrt(hatv2)))-
                    pnorm((gamm2[l2]-hatht02[yt02==l2])/c(sqrt(hatv2))))
    SUM<-SUM+as.double(sum(logphit0))
  }
  
  SUMk<-0
  for (l2 in 2:(L2-2)) {
    logJc<-log(1-gamm2[l2])+kappa2[l2-1]-2*log(1+exp(kappa2[l2-1]))
    SUMk<-SUMk+logJc
  }
  
  logP2=SUM+SUMk+(-0.5*t(kappa2-kappa20)%*%solve(K20)%*%(kappa2-kappa20))
  -logP2
}
##############################



##############################
#MCMC
#draw number and omit number
ndraw<-11000
#storage
sbeta0<-matrix(NA,ndraw,q) 
sbeta<-matrix(NA,ndraw,q) 
sdelta<-matrix(NA,ndraw,mS)
srho<-matrix(NA,ndraw,1)
#srho1<-matrix(NA,ndraw,1)
svsigm2<-matrix(NA,ndraw,S)
sV<-matrix(NA,ndraw,S*S)
skappa1<-matrix(NA,ndraw,L1-3)
skappa2<-matrix(NA,ndraw,L2-3)
sgamm1<-matrix(NA,ndraw,L1+1)
sgamm2<-matrix(NA,ndraw,L2+1)
sz<-matrix(NA,ndraw,nST)
szt0<-matrix(NA,ndraw,nS)
sLam<-matrix(NA,ndraw,S*S)

#necessary values
mTi<-solve(mT)
mTivc<-mTi%*%vc
D0i<-solve(D0)
D0id0<-solve(D0)%*%d0
nu1=L1-3
ac1=rep(0,L1-3)
nu2=L2-3
ac2=rep(0,L2-3)
st<-proc.time()
iter<-1
acc<-0
adj<-1901
zomit<-10
#detval<-sar_lndet(0,W,rmin,rmax)$detval#no approximation
detval<-sar_lndet(1,W,ilammin,ilammax)$detval#Chebyshev approximation
detvall<-detval[1:adj,] #Adjust according to the situation

InmV<-kronecker(diag(n),mV)
#InmV<-Matrix(kronecker(diag(n),mV),sparse = TRUE)
InmVi<-kronecker(diag(n),solve(mV))
#InmVi<-Matrix(kronecker(diag(n),solve(mV)),sparse=TRUE)
ITnmVi<-kronecker(diag(TM*n),solve(mV))
#ITnmVi<-Matrix(kronecker(diag(TM*n),solve(mV)),sparse=TRUE)

Imvsigm2<-kronecker(diag(m),diag(vsigm2))

Brho<-diag(mS)-rho*WIS
Arho<-diag(m)-rho*W

InLam<-kronecker(diag(n),Lam)
#InLam<-Matrix(kronecker(diag(n),Lam),sparse=TRUE)
ITInLam<-kronecker(diag(TM),InLam)
#ITInLam<-kronecker(diag(TM),InLam)

zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
zLam<-z-ITInLam%*%zj1


while(iter <= ndraw)
{
  #update z# #z should update first, otherwise there will be error
  lower<-rep(NA,nST)
  upper<-rep(NA,nST)
  lower1<-gamm1[y1]
  upper1<-gamm1[y1+1]
  lower2<-gamm2[y2]
  upper2<-gamm2[y2+1]
  lower[c(T,F)]<-lower1
  lower[c(F,T)]<-lower2
  upper[c(T,F)]<-upper1
  upper[c(F,T)]<-upper2
  
  #  InmVi<-kronecker(diag(n),solve(mV))
  #InmVi<-Matrix(kronecker(diag(n),solve(mV)),sparse=TRUE)
  InLammVi<-kronecker(diag(n),t(Lam)%*%solve(mV))
  #InLammVi<-Matrix(kronecker(diag(n),t(Lam)%*%solve(mV)),sparse=TRUE)
  #  InLam<-kronecker(diag(n),Lam)
  #InLam<-Matrix(kronecker(diag(n),Lam),sparse=TRUE)  
  Onedelta<-DeltaIS%*%delta
  
  Bzt<-solve(mV)+t(Lam)%*%solve(mV)%*%Lam
  InBzt<-kronecker(diag(n),Bzt)
  InBzti<-kronecker(diag(n),solve(Bzt))
  Ht<-Matrix(InBzt,sparse = TRUE)
  
  #t=1#
  ztj1tmp<-as.matrix(zt0)
  ztja1tmp<-as.matrix(c(z[(nS+1):(2*nS)]))
  
  xttmp<-x[1:nS,]
  xtja1tmp<-x[(nS+1):(2*nS),]
  
  Inbzt<-InmVi%*%(InLam%*%ztj1tmp+xttmp%*%beta+Onedelta)+InLammVi%*%(ztja1tmp-xtja1tmp%*%beta-Onedelta)
  ht<-InBzti%*%Inbzt
  zt1 <- as.double(rtmvnorm.sparseMatrix(n=1, mean=ht, H=Ht, 
                                         lower=lower[1:nS], upper=upper[1:nS], burn.in=zomit))
  z[1:nS,]<-zt1
  
  # for (t in 2:(TM-1)) {
  #   ztj1tmp<-as.matrix(z[(nS*(t-2)+1):(nS*(t-1))])
  #   ztja1tmp<-as.matrix(z[(nS*t+1):(nS*(t+1))])
  #   
  #   xttmp<-as.matrix(x[(nS*(t-1)+1):(nS*t),])
  #   xtja1tmp<-as.matrix(x[(nS*t+1):(nS*(t+1)),])
  #   
  #   
  #   Inbzt<-InmVi%*%(InLam%*%ztj1tmp+xttmp%*%beta+Onedelta)+InLammVi%*%(ztja1tmp-xtja1tmp%*%beta-Onedelta)
  #   ht<-InBzti%*%Inbzt
  #   zttmp <- as.double(rtmvnorm.sparseMatrix(n=1, mean=ht, H=Ht, 
  #                                            lower=lower[(nS*(t-1)+1):(nS*t)], upper=upper[(nS*(t-1)+1):(nS*t)], burn.in=zomit))
  #   z[(nS*(t-1)+1):(nS*t),]<-zttmp
  # }  
  
  #t=T#
  zTj1=z[(nS*(TM-2)+1):(nS*(TM-1)),]
  xT=x[(nS*(TM-1)+1):nST,]
  hT<-InLam%*%zTj1+xT%*%beta+DeltaIS%*%delta
  HT<-Matrix(InmVi,sparse = TRUE)
  zT <- as.double(rtmvnorm.sparseMatrix(n=1, mean=hT, H=HT, 
                                        lower=lower[(nS*(TM-1)+1):nST], upper=upper[(nS*(TM-1)+1):nST], burn.in=zomit))
  
  z[(nS*(TM-1)+1):nST,]<-zT
  
  zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
  zLam<-z-ITInLam%*%zj1
  
  
  #update beta 
  tst<-t(x)%*%ITnmVi
  Bi<-solve(tst%*%x+mTi)
  b<-tst%*%(zLam-OneTDeltaIS%*%delta)+mTivc
  beta<-t(rmvnorm(1,Bi%*%b,Bi))
  
  #update delta# ##transform sigm2 from matrix to numeric
  tst<-t(OneTDeltaIS)%*%ITnmVi
  Omega<-tst%*%OneTDeltaIS+t(Brho)%*%solve(Imvsigm2)%*%Brho+t(DeltaIS)%*%InmVi%*%DeltaIS
  Omegai<-solve(Omega)
  eta<- tst%*%(zLam-x%*%beta)+t(DeltaIS)%*%InmVi%*%(zt0-xt0%*%beta0)
  delta<-t(rmvnorm(1,Omegai%*%eta,Omegai))
  #  mdelta<-Omegai%*%eta
  #  vdelta<-Omegai
  #  mdelta1<-mdelta[c(T,F)]
  #  mdelta2<-mdelta[c(F,T)]
  #  vdelta1<-vdelta[c(T,F)]
  #  vdelta2<-vdelta[c(F,T)]
  #  delta1<-t(rmvnorm(1,mdelta1,vdelta1))
  #  delta2<-t(rmvnorm(1,mdelta2,vdelta2))
  #  delta<-matrix(NA,mS,1)
  #  delta[c(T,F)]<-delta1
  #  delta[c(F,T)]<-delta2
  
  #comment multirow£¬press Ctrl+shift+C 
  #metropolis step to get rho update and Brho, Arho update
  # accept<-0
  # rho2<-rho+a*rnorm(1)
  # while(accept==0) {
  #   if (rho2 >ilammin && rho2< ilammax) accept<-1
  #   else rho2<-rho+a*rnorm(1)
  # }
  # logprho<-logFCDrho(rho,delta,vsigm2,m,S,W)
  # logprho2<-logFCDrho(rho2,delta,vsigm2,m,S,W)
  # ratio<-exp(logprho2-logprho)
  # alpha<-min(ratio,1)
  # u<-runif(1)
  # if (u < alpha) {
  #   rho<-rho2
  #   acc<-acc+1
  # }
  # acc_rate<-acc/iter
  # 
  # if (acc_rate <0.4) a<-a/1.1 #it is necessary
  # if (acc_rate >0.6) a<-a*1.1 #it is necessary
  # Brho<-diag(mS)-rho*WIS
  # Arho<-diag(m)-rho*W
 
  
  # update rho using griddy Gibbs
  a1=a2=1
  e01=delta[c(T,F)]
  ed1=W%*%e01
  e02=delta[c(F,T)]
  ed2=W%*%e02
  
  # epe01 <- as.numeric(t(e01) %*% e01)
  # eped1 <- as.numeric(t(ed1) %*% ed1)
  # epe0d1<- as.numeric(t(ed1) %*% e01)
  # 
  # epe02 <- as.numeric(t(e02) %*% e02)
  # eped2 <- as.numeric(t(ed2) %*% ed2)
  # epe0d2<- as.numeric(t(ed2) %*% e02)
  
  epe01 <- as.double(crossprod(e01))  # slightly faster than t(e0) %*% e0
  eped1 <- as.double(crossprod(ed1))
  epe0d1<- as.double(crossprod(ed1, e01))
  
  epe02 <- as.double(crossprod(e02))  # slightly faster than t(e0) %*% e0
  eped2 <- as.double(crossprod(ed2))
  epe0d2<- as.double(crossprod(ed2, e02))
  
  nt<-m
  kt<-0
  rho <- draw_rho_multi(detvall,epe01,eped1,epe0d1,epe02,eped2,epe0d2,nt,kt,rho,a1,a2)
  Brho<-diag(mS)-rho*WIS
  Arho<-diag(m)-rho*W
  
  #update sigm2#  # sigma is different
  delta1<-delta[c(T,F)]
  delta2<-delta[c(F,T)]
  tst<-rchisq(S,m+2*e)
  tst1<-as.numeric(t(delta1)%*%t(Arho)%*%Arho%*%delta1)+2*g
  tst2<-as.numeric(t(delta2)%*%t(Arho)%*%Arho%*%delta2)+2*g
  vsigm2[1]<-tst1/tst[1]
  vsigm2[2]<-tst2/tst[2]
  Imvsigm2<-kronecker(diag(m),diag(vsigm2))
  #update sigm2# sigma is the same
  #   tst<-rchisq(1,mS+2*e)
  #   tst1<-as.numeric(t(delta)%*%t(Brho)%*%Brho%*%delta)+2*g
  #   vsigm2<-rep(tst1/tst,S)
  #   Imvsigm2<-kronecker(diag(m),diag(vsigm2))
  
  #update mV
  errt0<-zt0-xt0%*%beta0-DeltaIS%*%delta
  tmpt0<-errt0[1:2]%*%t(errt0[1:2])
  for (i in 2:n) {
    etmpt0<-errt0[(2*i-1):(2*i)]%*%t(errt0[(2*i-1):(2*i)])
    tmpt0<-tmpt0+etmpt0
  }
  
  err<-zLam-x%*%beta-OneTDeltaIS%*%delta
  tmp<-err[1:2]%*%t(err[1:2])
  for (i in 2:(n*TM)) {
    etmp<-err[(2*i-1):(2*i)]%*%t(err[(2*i-1):(2*i)])
    tmp<-tmp+etmp
  }
  hatn<-n0+n*(TM+1)
  hatQ<-Q0+tmp+tmpt0
  mVi<-rWishart(1, hatn, solve(hatQ))
  mVi<-mVi[,,1]
  mV<-as.matrix(solve(mVi))
  InmV<-kronecker(diag(n),mV)
  InmVi<-kronecker(diag(n),solve(mV))
  ITnmVi<-kronecker(diag(TM*n),solve(mV))
  
  
  #update kappa1 and gamma1#
  zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
  hatht0<-CMZt0(zt0,xt0)
  hath<-CMZ(z,zj1,Lam,x)
  kappa_o1<-kappa1
  opt_kappa1<-optim(kappa_o1,FCDkappa1,method="BFGS",hessian=TRUE)
  kappa_h1<-opt_kappa1$par
  Gh1<-solve(opt_kappa1$hessian)
  kappa_n1<-c(rmvt(1,kappa_h1,sigma=Gh1,df=nu1))
  while (kappa_n1>4 || kappa_n1<(-6)) {
    kappa_n1<-c(rmvt(1,kappa_h1,sigma=Gh1,df=nu1))
  }
  alpha1<--FCDkappa1(kappa_n1)+FCDkappa1(kappa_o1)+
    dmvt(t(kappa_o1),delta=kappa_h1,sigma=Gh1,df=nu1,log=TRUE)-
    dmvt(t(kappa_n1),delta=kappa_h1,sigma=Gh1,df=nu1,log=TRUE)
  
  if(log(runif(1))<alpha1){
    kappa1<-kappa_n1
    ac1<-ac1+1
  }
  
  for (l1 in 3:(L1-1)) {
    gamm1[l1]<-(gamm1[l1-1]+exp(kappa1[l1-2]))/(1+exp(kappa1[l1-2]))
  }
  
  #update kappa2 and gamma2#
  kappa_o2<-kappa2
  opt_kappa2<-optim(kappa_o2,FCDkappa2,method="BFGS",hessian=TRUE)
  kappa_h2<-opt_kappa2$par
  Gh2<-solve(opt_kappa2$hessian)
  kappa_n2<-c(rmvt(1,kappa_h2,sigma=Gh2,df=nu2))
  while (kappa_n2>4 || kappa_n2<(-6)) {
    kappa_n2<-c(rmvt(1,kappa_h2,sigma=Gh2,df=nu2))
  }
  alpha2<--FCDkappa2(kappa_n2)+FCDkappa2(kappa_o2)+
    dmvt(t(kappa_o2),delta=kappa_h2,sigma=Gh2,df=nu2,log=TRUE)-
    dmvt(t(kappa_n2),delta=kappa_h2,sigma=Gh2,df=nu2,log=TRUE)
  
  if(log(runif(1))<alpha2){
    kappa2<-kappa_n2
    ac2<-ac2+1
  }
  
  for (l2 in 3:(L2-1)) {
    gamm2[l2]<-(gamm2[l2-1]+exp(kappa2[l2-2]))/(1+exp(kappa2[l2-2]))
  }
  
  
  #update gamma
  #z1<-z[c(T,F),]
  #for (l1 in 2:(L1-2)){
  #  left<-max(max(z1[y1==l1]),gamm1[l1-1+1])
  #  right<-min(min(z1[y1==l1+1]),gamm1[l1+1+1])
  # gamm1[l1+1]<-runif(1,left,right)
  #}
  
  #z2<-z[c(F,T),]
  #for (l2 in 2:(L2-2)){
  #  left<-max(max(z2[y2==l2]),gamm2[l2-1+1])
  #  right<-min(min(z2[y2==l2+1]),gamm2[l2+1+1])
  # gamm2[l2+1]<-runif(1,left,right)
  #}
  #update zt0#
  lower<-rep(NA,nS)
  upper<-rep(NA,nS)
  lower1<-gamm1[yt01]
  upper1<-gamm1[yt01+1]
  lower2<-gamm2[yt02]
  upper2<-gamm2[yt02+1]
  lower[c(T,F)]<-lower1
  lower[c(F,T)]<-lower2
  upper[c(T,F)]<-upper1
  upper[c(F,T)]<-upper2
  
  zt1<-as.matrix(z[1:nS,])
  xt1<-as.matrix(x[1:nS,])
  Bz0<-t(Lam)%*%solve(mV)%*%Lam+D0i+solve(mV)
  Bz0i<-solve(Bz0)
  InBz0<-kronecker(diag(n),Bz0)
  InBz0i<-kronecker(diag(n),solve(Bz0))
  Ht0<-Matrix(InBz0,sparse = TRUE)
  
  InLammVi<-kronecker(diag(n),t(Lam)%*%solve(mV))
  #  InLammVi<-Matrix(kronecker(diag(n),t(Lam)%*%solve(mV)),sparse=TRUE)
  Inbz0<-InLammVi%*%(zt1-xt1%*%beta-DeltaIS%*%delta)+as.matrix(rep(D0id0,n))+InmVi%*%(xt0%*%beta0+DeltaIS%*%delta)
  ht0<-InBz0i%*%Inbz0
  zt0 <- as.double(rtmvnorm.sparseMatrix(n=1, mean=ht0, H=Ht0, 
                                         lower=lower, upper=upper, burn.in=zomit))
  
  zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
  zLam<-z-ITInLam%*%zj1
  
  #update beta0 
  tst<-t(xt0)%*%InmVi
  Bi<-solve(tst%*%xt0+mTi)
  b<-tst%*%(zt0-DeltaIS%*%delta)+mTivc
  beta0<-t(rmvnorm(1,Bi%*%b,Bi))
  
  #update Lam#
  hatv1<-mV[1,1]-mV[1,-1]%*%solve(mV[-1,-1])%*%mV[-1,1]
  hatv2<-mV[2,2]-mV[2,-2]%*%solve(mV[-2,-2])%*%mV[-2,2]
  
  A<-matrix(NA,nST,S)
  zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
  deltatmp<-OneTDeltaIS%*%delta
  
  #Lam[1,1] s=1,g0=1#
  s=1
  A[c(T,F),1]<-z[c(T,F)]-(Lam[1,2]*zj1[c(F,T)] +x[c(T,F),]%*%beta+deltatmp[c(T,F)]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[c(F,T)]-Lam[2,1]*zj1[c(T,F)]-Lam[2,2]*zj1[c(F,T)]-x[c(F,T),]%*%beta-deltatmp[c(F,T)]))
  AZ<-A[c(T,F),1]*zj1[c(T,F)]
  Z2<-zj1[c(T,F)]^2
  mnl<-sum(AZ)/sum(Z2)
  vl<-hatv1/sum(Z2)
  Lam[1,1]<-rnorm(1,mnl,vl)
  accept<-0
  while(accept==0) {
    if (abs(Lam[1,1]) < 1  && abs(Lam[1,1])+abs(Lam[1,2])< 1) accept<-1
    #    if (abs(Lam[1,1]) < 1) accept<-1
    else Lam[1,1]<-rnorm(1,mnl,vl)
  }
  
  #Lam[2,2] s=2,g0=2#
  s=2
  A[c(F,T),2]<-z[c(F,T)]-(Lam[2,1]*zj1[c(T,F)] +x[c(F,T),]%*%beta+deltatmp[c(F,T)]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[c(T,F)]-Lam[1,1]*zj1[c(T,F)]-Lam[1,2]*zj1[c(F,T)]-x[c(T,F),]%*%beta-deltatmp[c(T,F)]))
  AZ<-A[c(F,T),2]*zj1[c(F,T),]
  Z2<-zj1[c(F,T)]^2
  mnl<-sum(AZ)/sum(Z2)
  vl<-hatv2/sum(Z2)
  Lam[2,2]<-rnorm(1,mnl,vl)
  accept<-0
  while(accept==0) {
    if (abs(Lam[2,2]) < 1  && abs(Lam[2,2])+abs(Lam[2,1])< 1) accept<-1
    else Lam[2,2]<-rnorm(1,mnl,vl)
  }
  
  #Lam[1,2] s=1,g0=2#
  s=1
  A[c(T,F),2]<-z[c(T,F)]-(Lam[1,1]*zj1[c(T,F)] +x[c(T,F),]%*%beta+deltatmp[c(T,F)]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[c(F,T)]-Lam[2,1]*zj1[c(T,F)]-Lam[2,2]*zj1[c(F,T)]-x[c(F,T),]%*%beta-deltatmp[c(F,T)]))
  AZ<-A[c(T,F),2]*zj1[c(F,T),]
  Z2<-zj1[c(F,T)]^2
  mnl<-sum(AZ)/sum(Z2)
  vl<-hatv1/sum(Z2)
  Lam[1,2]<-rnorm(1,mnl,vl)
  accept<-0
  while(accept==0) {
    if (abs(Lam[1,2]) < 1  && abs(Lam[1,2])+abs(Lam[1,1])< 1) accept<-1
    else Lam[1,2]<-rnorm(1,mnl,vl)
  }
  
  #Lam[2,1] s=2,g0=1#
  s=2
  A[c(F,T),1]<-z[c(F,T)]-(Lam[2,2]*zj1[c(F,T)] +x[c(F,T),]%*%beta+deltatmp[c(F,T)]+as.double(mV[s,-s]%*%solve(mV[-s,-s]))*(z[c(T,F)]-Lam[1,1]*zj1[c(T,F)]-Lam[1,2]*zj1[c(F,T)]-x[c(T,F),]%*%beta-deltatmp[c(T,F)]))
  
  AZ<-A[c(F,T),1]*zj1[c(T,F),]
  Z2<-zj1[c(T,F)]^2
  mnl<-sum(AZ)/sum(Z2)
  vl<-hatv2/sum(Z2)
  Lam[2,1]<-rnorm(1,mnl,vl)
  accept<-0
  while(accept==0) {
    if (abs(Lam[2,1]) < 1  && abs(Lam[2,1])+abs(Lam[2,2])< 1) accept<-1
    #    if (abs(Lam[2,1]) < 1) accept<-1
    else Lam[2,1]<-rnorm(1,mnl,vl)
  }
  
  InLam<-kronecker(diag(n),Lam)
  ITInLam<-kronecker(diag(TM),InLam)
  
  zj1<-as.matrix(c(zt0,z[1:(nS*(TM-1))]))
  zLam<-z-ITInLam%*%zj1
  
  
  
  sbeta0[iter,1:q]<-t(beta0)
  sbeta[iter,1:q]<-t(beta)
  sdelta[iter,1:mS]<-t(delta)
  srho[iter,1]<-rho
  svsigm2[iter,1:S]<-vsigm2
  sV[iter,1:(S*S)]<-c(mV)
  skappa1[iter,1:(L1-3)]<-t(kappa1)
  skappa2[iter,1:(L2-3)]<-t(kappa2)
  sgamm1[iter,1:(L1+1)]<-t(gamm1)
  sgamm2[iter,1:(L2+1)]<-t(gamm2)
  sz[iter,1:nST]<-t(z)
  szt0[iter,1:nS]<-t(zt0)
  sLam[iter,1:(S*S)]<-c(Lam)
  
  cat(100*iter/ndraw,"% done \n",sep="")
  iter <- iter + 1
}

cat(" End of MCMC \n",date(),"\n")
en<-proc.time()
gtime<-en-st
nomit<-1000
mnbeta0<-apply(sbeta0[(nomit+1):ndraw,],2,mean)
mnbeta<-apply(sbeta[(nomit+1):ndraw,],2,mean)
mndelta<-apply(sdelta[(nomit+1):ndraw,],2,mean)
mnrho<-mean(srho[(nomit+1):ndraw,])
mnvsigm2<-apply(svsigm2[(nomit+1):ndraw,],2,mean)
mnmV<-apply(sV[(nomit+1):ndraw,],2,mean)
mngamm1<-apply(sgamm1[(nomit+1):ndraw,],2,mean)
mngamm2<-apply(sgamm2[(nomit+1):ndraw,],2,mean)
mnz<-apply(sz[(nomit+1):ndraw,],2,mean)
mnzt0<-apply(szt0[(nomit+1):ndraw,],2,mean)
mnLam<-apply(sLam[(nomit+1):ndraw,],2,mean)

summaryresults=function(x){
  c(mean=mean(x),sd=sd(x),quantile(x,c(0.025,0.05,0.5,0.95,0.975)),
    effSize=effectiveSize(x))
}
sumbeta0<-round(apply(sbeta0[(nomit+1):ndraw,],2,summaryresults),4)
sumbeta<-round(apply(sbeta[(nomit+1):ndraw,],2,summaryresults),4)
sumdelta<-round(apply(sdelta[(nomit+1):ndraw,],2,summaryresults),4)
sumrho<-t(round(summaryresults(srho[(nomit+1):ndraw]),4))
sumvsigm2<-round(apply(svsigm2[(nomit+1):ndraw,],2,summaryresults),4)
sumV<-round(apply(sV[(nomit+1):ndraw,],2,summaryresults),4)
if (L1>3) {sumgamm1<-round(apply(sgamm1[(nomit+1):ndraw,3:(L1-1)],2,summaryresults),4)}
if (L2>3) {sumgamm2<-round(apply(sgamm2[(nomit+1):ndraw,3:(L2-1)],2,summaryresults),4)}
sumLam<-round(apply(sLam[(nomit+1):ndraw,],2,summaryresults),4)

#MSEbeta0<-mean(apply((sbeta0[(nomit+1):ndraw,]-matrix(rep(tbeta0,ndraw-nomit),ndraw-nomit, byrow = TRUE))^2,1,sum))
#MSEbeta<-mean(apply((sbeta[(nomit+1):ndraw,]-matrix(rep(tbeta,ndraw-nomit),ndraw-nomit, byrow = TRUE))^2,1,sum))

#save.image("results.Rdata")
#save.image("results-gibby rho.Rdata")
save(sumbeta,sumdelta,sumLam,sumrho,sumvsigm2,sumV,sumgamm1,sumgamm2,file = "results.Rdata")


#sumbeta0
sumbeta
sumdelta
sumLam
sumrho
sumvsigm2
sumV
sumgamm1
sumgamm2

# rm(list = ls())

