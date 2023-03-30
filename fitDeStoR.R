####################################################################################################
### 	 R program for multivariate density deconvolution  by stochastic rotation of replicates	   ### 
###	       (last modified on Dec 15, 2020, at The University of Florida, Gainesville)            ###
####################################################################################################


setwd("") # Set working directory

###################
### Set Seed No ###	### Seed must be set inside the loops
###################

seedno <- 1
set.seed(seedno)  


######################################
### Add Libraries and Source Files ###
######################################


# These libraries must be available to R
library(corpcor)
library(corrplot)
library(foreach)
library(mvtnorm) 
library(movMF)
library(MCMCpack) 
library(msm)
library(RColorBrewer)

source("Decon_Functs.R")




#################################
### Set Simulation Parameters ###
#################################

seedno = 1
set.seed(seedno)


#################
### Load Data ###
#################

d <- 3
N <- 3*1000
n <- 1000

ws <- read.csv("Synthetic_Data.csv", header = T)
ws <- as.matrix(ws[, - 1])
ws <- t(ws)
###########################
### Set Size Parameters ###
###########################

reps <- 3
mis = rep(reps,times=n)
inds <- rep(1:n,times=mis)
ids <- inds#inds
N = sum(mis)

inds1 = inds2 = numeric(n)
inds1[1] = 1
inds2[1] = inds1[1]+mis[1]-1
for(ii in 2:n)
{
  inds1[ii] = inds1[ii-1]+mis[ii-1]
  inds2[ii] = inds1[ii] + mis[ii]-1
}

#########################################################################
### Adjust Range to Bring the Replicates for All Components in [0,20] ###
#########################################################################

range.xs = 20
x.lwr = rep(0,d)
x.upr = rep(10,d)
# ws <- ws.true
for(jj in 1:d)
{
  ws[jj,] <- range.xs*(ws[jj,]-min(ws[jj,]))/(max(ws[jj,])-min(ws[jj,]))
}
ws.obs <- ws
ws.true <- ws


W <- array(0, dim = c(3, d, N/d))
for(i in 1:3){
  W[i,,] <- ws.true[, (1:n - 1)*3 + i]
}

################################################################
### Run Univariate Samplers First to Generate Initial Values ###
################################################################

#################################
### Priors and Initial Values ###
#################################

### Initialization and prior of xs and us
inds = rep(1:n,times=mis)
xs = current.xs = proposed.xs = start.xs = s2is = wbars = matrix(0,d,n)
xs.trans = matrix(0,d,n)
us = matrix(0,d,N)
x.grid.length = 500
x.grid = matrix(0,d,x.grid.length)
for(jj in 1:d)
{
  xs[jj,] = wbars[jj,] = tapply(ws[jj,],inds,"mean")
  xs[jj,xs[jj,]>x.upr[jj]] = max(xs[jj,xs[jj,]<x.upr[jj]])
  s2is[jj,] = as.vector(tapply(ws[jj,],inds,var))
  x.grid[jj,] = seq(x.lwr[jj],x.upr[jj],length=x.grid.length)
}
range.start.xs = apply(t(apply(xs,1,range)),1,diff)

### Prior and initialization for mixture
alpha.x = 0.1

# Independent Normal Inverse-Gamma
if(d==1)
{
  mu0.x = mean(apply(as.matrix(t(xs[1:d,])),1,mean))
  sigmasq0.x = mean(apply(as.matrix(t(xs[1:d,])),1,var))
}else
{
  mu0.x = mean(apply(xs[1:d,],1,mean))
  sigmasq0.x = mean(apply(xs[1:d,],1,var))
}
a.sigmasq.x = 1
b.sigmasq.x = 1

z.x.max = 10
z.x = matrix(0,d,n)
for(jj in 1:d)
  z.x[jj,] = sample(1:z.x.max,size=n,replace=TRUE,prob=rep(1/z.x.max,z.x.max))
mu.x = sigmasq.x = matrix(0,d,z.x.max)
pi.x = matrix(0,d,z.x.max)
n.kk.x = matrix(0,d,z.x.max)
prob.x = pi.x
params.x = array(0,dim=c(d,z.x.max,2))

### Prior and initialization of s2t and thetas
K.t = 10 			# num of knots
alpha.t = 100
beta.t = 1
s2t = rep(0.01,d)
P.t = P.mat(K.t+1) 	# penalty matrix
knots.t = matrix(0,d,K.t)
delta.t = numeric(d)
var.es = numeric(d)
thetas = current.thetas = proposed.thetas = matrix(0,d,K.t+1)
vars = current.vars = proposed.vars = matrix(0,d,n)
prop.sig.thetas = array(0,dim=c(d,K.t+1,K.t+1))
var.grid = matrix(0,d,100)

### Prior and initialization for mixture
z.u.max = 10
z.u = matrix(1,d,N)
alpha.u = 0.1
a.u = 3
b.u = 1
sigma0.u = 4
params.u = array(0,dim=c(d,z.u.max,4))
for(jj in 1:d)
  params.u[jj,,] = matrix(c(0.5,0,1,1),nrow=z.u.max,ncol=4,byrow=T)		# unique values
pi.u = matrix(0,d,z.u.max) #matrix(1/z.u.max,d,z.u.max)
n.kk.u = matrix(0,d,z.u.max)
prob.u = pi.u
e.grid = seq(-4,4,length=100)
density.e.est = matrix(0,d,length(e.grid))
simsize.mh.u = 10

simsize_univ = 1000
burnin_univ = 500
z.x.max.univ = 5      # must be <= z.x.max
z.u.max.univ = 5      # must be <= z.u.max

print("Running univariate models! Thanks for your patience!")
for(jj in 1:d)
{
  filename = paste("Simulation",jj,"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
  filename = paste(filename,".RData",sep="") 
  univ_results = UNIV_DECON_REGULAR(ws[jj,],x.lwr[jj],x.upr[jj],mis,z.x.max.univ,z.u.max.univ,K.t,simsize_univ,burnin_univ,FALSE,FALSE) 
  save(file=filename,univ_results)
}
print("Done running univariate models for the components!")
for(jj in 1:d)
{
  filename = paste("Simulation",jj,"simsize",simsize_univ,"burnin",burnin_univ,"seedno",seedno,sep="_")
  filename = paste(filename,".RData",sep="") 
  print(filename)
  load(filename)
  
  xs[jj,] = univ_results$xs
  us[jj,] = univ_results$us
  knots.t[jj,] = univ_results$knots
  delta.t[jj] = knots.t[jj,2]-knots.t[jj,1]  
  thetas[jj,] = current.thetas[jj,] = univ_results$thetas
  
  z.x[jj,] = univ_results$z.xs
  pi.x[jj,1:z.x.max.univ] = univ_results$pi.xs
  params.x[jj,1:z.x.max.univ,] = univ_results$params.xs
  
  mu.x[jj,] = params.x[jj,,1]
  sigmasq.x[jj,] = params.x[jj,,2]
  
  z.u[jj,] = univ_results$z.us
  pi.u[jj,1:z.u.max.univ] = univ_results$pi.us
  params.u[jj,1:z.u.max.univ,] = univ_results$params.us 
}

for(jj in 1:d)
{
  current.xs[jj,] = start.xs[jj,] = xs[jj,]
  temp.vars = B.basis(xs[jj,],knots.t[jj,])%*%exp(thetas[jj,])
  temp.vars[temp.vars==0] = max(temp.vars)
  vars[jj,] = current.vars[jj,] = temp.vars
  prop.sig.thetas[jj,,] = make.positive.definite(prop.sig.thetas.fn(thetas[jj,],xs[jj,],mis,us[jj,],s2t[jj],K.t,P.t,knots.t[jj,],n))
  var.grid[jj,] = seq(min(knots.t[jj,]),max(knots.t[jj,]),length=length(var.grid[jj,]))
}  

### Redefine shared component parameters mu.x and sigma.sq.x running a small sampler keeping the xs fixed

mu.x.new = numeric(z.x.max)
mu.x.min = Inf
mu.x.max = -Inf
for(jj in 1:d)
{
  mu.x.min = min(mu.x.min,mu.x[jj,z.x[jj,]])
  mu.x.max = max(mu.x.max,mu.x[jj,z.x[jj,]])
}
mu.x.new = seq(mu.x.min,mu.x.max,len=z.x.max)
for(jj in 1:d)
{
  TempMat = abs(matrix(mu.x[jj,z.x[jj,]],n,z.x.max)-matrix(rep(mu.x.new,n),n,z.x.max,byrow=T))
  z.x[jj,] = apply(TempMat,1,which.min)
}
mu.x = mu.x.new
sigmasq.x = rep(2*var(as.vector(xs))/z.x.max,z.x.max)
for(iii in 1:100)
{
  for(jj in 1:d)
  {
    xs.trans.temp.num = pnorm((xs[jj,]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
    xs.trans.temp.den = pnorm((x.upr[jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))-pnorm((x.lwr[jj]-mu.x[z.x[jj,]])/sqrt(sigmasq.x[z.x[jj,]]))
    xs.trans.temp.den[xs.trans.temp.den<=0] = 0.001  # Numerical Stability
    xs.trans.temp = xs.trans.temp.num/xs.trans.temp.den
    xs.trans.temp[xs.trans.temp>=1] = 0.999  # Numerical Stability
    xs.trans.temp[xs.trans.temp<=0] = 0.001  # Numerical Stability
    xs.trans.temp = mu.x[z.x[jj,]]+sqrt(sigmasq.x[z.x[jj,]])*qnorm(xs.trans.temp)
    xs.trans.temp[xs.trans.temp<(x.lwr[jj]-10)] = x.lwr[jj] - 10  # Numerical Stability
    xs.trans.temp[xs.trans.temp>(x.upr[jj]+10)] = x.upr[jj] + 10  # Numerical Stability
    xs.trans[jj,] = xs.trans.temp
  }
  
  for(kk in 1:z.x.max)
  {
    xspool = NULL
    for(jj in 1:d)
    {
      temp = which(z.x[jj,]==kk)
      xspool = c(xspool,xs.trans[jj,temp])
    }
    
    sigmasq.temp = 1/(sum(n.kk.x[,kk])/sigmasq.x[kk] + 1/sigmasq0.x)
    mu.temp = (sum(xspool)/sigmasq.x[kk] + mu0.x/sigmasq0.x) * sigmasq.temp
    mu.x[kk] = rnorm(1,mu.temp,sqrt(sigmasq.temp))	
    
    post.a.sigmasq.x = a.sigmasq.x + length(xspool)/2
    post.b.sigmasq.x = b.sigmasq.x + sum((xspool-mu.x[kk])^2)/2
    sigmasq.x[kk] = 1/rgamma(1,shape=post.a.sigmasq.x,rate=post.b.sigmasq.x)
    
    density.numeric.check = dtnorm(x.grid[jj,],mu.x[kk],sqrt(sigmasq.x[kk]),x.lwr[jj],x.upr[jj])
    if((sum(is.nan(density.numeric.check))>0) || (sum(is.infinite(density.numeric.check))>0))
    {
      mu.x[kk] = mean(wbars)
      sigmasq.x[kk] = 10^10
    }
  }
  
  for(jj in 1:d)
  {
    for(ii in 1:n)
    {
      prob.x = pi.x[jj,] * dtnorm(rep(xs[jj,ii],z.x.max),mu.x,sqrt(sigmasq.x),lower=x.lwr[jj],upper=x.upr[jj])
      prob.x[is.nan(prob.x)]=0; 	prob.x[is.infinite(prob.x)]=max(prob.x[is.finite(prob.x)])   # Numerical Stability
      z.x[jj,ii] = sample(z.x.max,1,TRUE,prob.x)   # New z.x[ii] drawn
    }
    n.kk.x[jj,] = tabulate(z.x[jj,],nbins=z.x.max)
    pi.x[jj,] = rdirichlet(1,alpha.x/z.x.max+n.kk.x[jj,])
  }
}





#########################################################################
### Run Multivariate Sampler for Deconvolution by Stochastic Rotation ###
#########################################################################

#Initialize parameteres of the truncated normal mixture
muK   <- mu.x
sigK  <- sigmasq.x
IndK  <- z.x
probK <- pi.x

#Initialize unobserved X from univariate estimates
Xini  <- xs
kappaM <- matrix(0, d, n)

#Initialize kappa using movMF 
for(i in 1:n){
  mat <- W[,,i]
  mat <- mat / sqrt(rowSums(mat^2))
  temp <- movMF(mat, k=1)
  kappaM[, i] <- sqrt(sum(Xini[, i]^2))*array(temp$theta)/Xini[, i]
}

uniout <- list(mu = muK, sig = sigK, Ind = IndK, prob = probK)

fit  <- DeStoR(W, Xini =  Xini, Thini = kappaM, uniout = uniout, A = 0, B = 10, Kx = 10, knot = 10, norder = 4, Total_itr = 5000, burn = 3000, ncore = 1)
