#' @title The function to estimate deconvolution density based on the model
#' @description Takes the multivariate data and outputs from a univariate models and returns MCMC samples of the parameters
#' @references Roy and Sarkar (2021)
#'
#' @param W is the multivariate data replicates X dimension X individuals
#' @param Xini is the matrix of initial estimates of X from univariate sampler
#' @param Thini is the matrix of initial estimates of \kappa
#' @param uniout is the estimates of mixture components from the univariate sampler
#' @param Total_sample is the total number of iterations of MCMC
#' @param burn is the number of burn-in MCMC samples
#' @param Kx is the number of mixture components
#' @param A is the lower limit of truncated normal
#' @param B is the upper limit of truncated normal
#' @param knot is the number of knots in B-spline
#' @param norder is the order of B-spline, for example, norder=4 fits cubic splines
#' @param ncore is the number of machine cores to be used, set ncore=1 for Windows machines

#' @return DeStoR returns a list of the posterior samples of unobserved X, B-splines coefficients and mixture normal parameters



DeStoR <- function(W, Xini = NULL, Thini = NULL, uniout = NULL, Total_itr = 5000, burn = 2500, Kx = 5, A = 0, B = 10, knot = 6, norder = 4, ncore = 1)
{
  set.seed(1)
  library(fda)
  library(gtools)
  library(parallel)
  library(SphericalCubature)
  library(nnls)
  
  noupX = 500 
  noupT = 200
  nbl = 1
  sdX=1e-4
  
  
  #HMC sampler for B-spline coef for s^2()
  HMC = function (U, grad_U, epsilon, L = 30, current_q, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = q*(q>0)
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCth = function (U, grad_U, epsilon, L = 30, current_q, l, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      #if(TM>-Inf){
      q = q*(q>TM[l,]) + TM[l,]*(q<TM[l,])
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, l)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, l)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, l)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCX = function (U, grad_U, epsilon, L = 30, current_q, l, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = q * (q>=A) * (q <= B) + A*(q<A) + B*(q>B)
      
      q[which(q<1)] <- q[which(q<1)]+ 1e-10
      #q[which(q>5)] <- q[which(q>5)] - 1e-10
      
      if(length(which(q==Inf))>0)
        q[which(q==Inf)] = rep(B, length(which(q==Inf)))
      
      if(length(which(q==-Inf))>0)
        q[which(q==-Inf)] = rep(A, length(which(q==-Inf)))
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, l)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, l)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, l)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCthmu = function (U, grad_U, epsilon, L = 30, current_q, l, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, l)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, l)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, l)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCsig = function (U, grad_U, epsilon, L = 30, current_q, l, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = q*(q>0)
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, l)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, l) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, l)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, l)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  MvMFconst <- function(k, kappa, p=n2){
    kappaf <- kappa[, k]
    part1  <- 2^{-p*(p+5)/4+p^2/2}/{pi^(p/2)}
    part2  <- exp(sum(kappaf))
    part3  <- prod(gamma((p-1:p+1)/2))
    part4i <- sqrt(kappaf[index[1, ]] + kappaf[index[2, ]])
    part4  <- 1/prod(part4i)
    
    ret <- part1 * part2 * part3 * part4
    
    return(ret)
  } 
  
  MvMFconstder <- function(k, kappa, p = n2){
    kappaf <- kappa[, k]
    part1  <- 2^{-p*(p+5)/4+p^2/2}/{pi^(p/2)}
    part2  <- exp(sum(kappaf))
    part3  <- prod(gamma((p-1:p+1)/2))
    part4i <- sqrt(kappaf[index[1, ]] + kappaf[index[2, ]])
    part4  <- 1/prod(part4i)
    
    ret1 <- part1 * part2 * part3 * part4
    
    kappaf <- kappa[, k]
    part1 <- {-p*(p+5)/4+p^2/2}/{pi^(p/2)} * log(2)
    part2 <- sum(kappaf)
    part3 <- sum(digamma((p-1:p+1)/2))
    part4 <- -(0.5)*sum(log(kappaf[index[1, ]] + kappaf[index[2, ]]))
    
    ret2 <- part1 + part2 + part3 + part4
    
    return(ret1*ret2)
  }
  
  logMvMFconst <- function(k, kappa, p=n2){
    kappaf <- kappa[, k]
    #part1 <- {-p*(p+5)/4+p^2/2}/{pi^(p/2)} * log(2)
    part2 <- sum(kappaf)
    #part3 <- sum(digamma((p-1:p+1)/2))
    part4 <- -(0.5)*sum(log(kappaf[index[1, ]] + kappaf[index[2, ]]))
    
    ret <- part2 + part4
    
    return(ret)
  } 
  
  meanMvMF <- function(k, kappa){
    kappaf <- kappa[, k]
    part2 <- rep(1, n2)
    #part2[j] <- 1
    part4 <- rep(0, n2)
    
    for(i in 1:n2){
      part4[i] <- -(0.5)*sum(1/(kappaf[i] + kappaf[-i])) 
    }
    
    ret <- part2 + part4
    
    return(ret)
  } 
  
  meanMvMFder <- function(k, j, kappa){
    kappaf <- kappa[, k]
    part2 <- rep(0, n2)
    #part2[j] <- 1
    part4 <- rep(0, n2)
    
    for(i in 1:n2){
      part4[i] <- +(0.5)*sum(1/(kappaf[i] + kappaf[j])^2) 
    }
    
    ret <- part2 + part4
    
    return(ret)
  } 
  
  mathcal <- function(x){
    e1    <- rep(0, length(x))
    e1[1] <- 1
    
    Ac <- x + e1*sqrt(sum(x^2))
    A1 <- diag(length(x)) - 2*tcrossprod(Ac)/sum(Ac^2)
    
    return(A1)
  }
  
  mathcalder <- function(x){
    e1    <- rep(0, length(x))
    e1[1] <- 1
    ret   <- matrix(0, n2, n2)
    A1deri <- 4*tcrossprod(Ac)%*%Ac/(sum(Ac^2))^2 - 4*Ac/(sum(Ac^2))
    for(i in 1:n2)
      ret[,i] <- A1deri + rep(1, n2) + 0.5*e1*x[i]*(sum(x^2))^(-1/2)
    
    return(ret)
  }
  
  orthtr <- function(l, k, mat1, mat2){
    x  <- mat1[1:n2+n2*(l-1), k]
    y  <- mat2[1:n2+n2*(l-1)]
    
    Ac <- x + e1*sqrt(sum(x^2))
    A1 <- diag(length(x)) - 2*tcrossprod(Ac)/sum(Ac^2)
    Bc <- y + e1*sqrt(sum(y^2))
    B1 <- diag(length(y)) - 2*tcrossprod(Bc)/sum(Bc^2)
    
    # xn <- x / sqrt(sum(x^2))
    # yn <- y / sqrt(sum(y^2))
    # 
    # v <- xn-yn
    # v <- v/sqrt(sum(v^2))
    # mat <- diag(n2) - 2*v%*%t(v)
    s    <- sqrt(sum(x^2))/sqrt(sum(y^2))
    #return(diag(mat))
    return(s*colSums(A1*B1))
  }
  
  orthtrder <- function(l, k, i, mat1, mat2){
    x  <- mat1[1:n2+n2*(l-1), k]
    e1    <- rep(0, length(x))
    e1[1] <- 1
    ret   <- rep(0, n2)
    Ac <- x + e1*sqrt(sum(x^2))
    A1deri <- array(4*tcrossprod(Ac)%*%Ac/(sum(Ac^2))^2) - array(4*Ac/(sum(Ac^2)))
    A1 <- diag(length(x)) - 2*tcrossprod(Ac)/sum(Ac^2)
    
    y  <- mat2[1:n2+n2*(l-1)]
    Bc <- y + e1*sqrt(sum(y^2))
    B1 <- diag(length(y)) - 2*tcrossprod(Bc)/sum(Bc^2)
    
    
    s     <- sqrt(sum(x^2))/sqrt(sum(y^2))
    sderi <- x[i]/sqrt(sum(x^2))
    sderi <- sderi / sqrt(sum(y^2))
    
    
    temp     <- e1*x[i]*(sum(x^2))^(-1/2)
    temp[i]  <- temp[i] + 1
    tempm    <- array(temp) %*% t(array(A1deri))
    ret      <- colSums(tempm*B1)
    
    return(s*ret + sderi*colSums(A1*B1))
  }
  
  orthtrder1 <- function(l, k, i, mat1, mat2){
    x  <- mat2[1:n2+n2*(l-1)]
    e1    <- rep(0, length(x))
    e1[1] <- 1
    #ret   <- matrix(0, n2, n2)
    Ac <- x + e1*sqrt(sum(x^2))
    A1deri <- array(4*tcrossprod(Ac)%*%Ac/(sum(Ac^2))^2) - array(4*Ac/(sum(Ac^2)))
    
    y  <- mat1[1:n2+n2*(l-1), k]
    Bc <- y + e1*sqrt(sum(y^2))
    B1 <- diag(length(y)) - 2*tcrossprod(Bc)/sum(Bc^2)
    
    temp     <- e1*x[i]*(sum(x^2))^(-1/2)
    temp[i]  <- temp[i] + 1
    tempm    <- array(temp) %*% t(array(A1deri))
    ret       <- colSums(tempm*B1)
    return(ret)
  }
  
  Uth <- function(x, i){
    thetac      <- theta
    thetac[i, ] <- x#exp(x)
    BSXi        <- bsplineS(X[i, ],  breaks=seq(A,B,(B-A)/knot), norder=norder)
    
    kappac      <- kappa 
    kappac[i, ] <- BSXi %*% thetac[i, ]
    
    temp <- mcmapply(1:n3, FUN = meanMvMF, MoreArgs = list(kappa=kappac), mc.cores = ncore)
    
    Wmat <- t(matrix(array(W), n1))
    Wmat <- Wmat * rep(array(temp), n1)
    Wc <- array(t(Wmat), dim=c(n1, n2, n3))
    
    for(k in 1:n1){
      Wnorm[, k] <- array(sqrt(colSums(Wc[k,,]^2)))
    }
    
    M1 <- matrix(array(Wmat/Wnorm[inds, ]), n2*n3, n1) #* rep(array(temp), n1)
    #M1 <- M1/Umat[inds, ]
    
    tempmat <- array(0, dim = c(n2, n3, n1))
    
    Xnorm <- sqrt(colSums(X^2))/n2
    
    Xn <- X / matrix(n2*Xnorm, n2, n3, byrow = T)
    
    for(l in 1:n1){
      tempmat[,,l] <- mcmapply(1:n3, FUN = orthtr, MoreArgs = list(k = l, mat1=M1, mat2 = array(Xn)))
    }
    
    tempcst <- mcmapply(1:n3, FUN = logMvMFconst, MoreArgs = list(kappa=kappac), mc.cores = ncore)
    
    part1 <- n1 * sum(tempcst)
    part2 <- - sum(array(tempmat) * rep(array(kappac), n1)) #* rep(array(temp), n1))
    
    ret <- part1 + part2 + sum((x-mth[i,])^2/(2*dth[i]))#/(2*100)
    
    #Xnorm <- sqrt(colSums(X^2))/n2
    #BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
    #Uvar <- BSxn %*% beta#exp(x)
    Umat  <- (Wnorm/n2) / matrix(Xnorm, n3, n1)
    ret <- ret + sum(rowSums((log(Umat)+matrix(Uvar/2, n3, n1))^2)/Uvar)/2
    
    return(ret)
  }
  
  Ubt <- function(x){
    Xnorm <- sqrt(colSums(X^2))/n2
    BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
    Uvarc <- BSxn %*% x#exp(x)
    Umat  <- (Wnorm/n2) / matrix(Xnorm, n3, n1)
    ret <- sum(rowSums((log(Umat)+matrix(Uvarc/2, n3, n1))^2)/Uvarc) + n1*sum(log(Uvarc)) + sum((x-mbt)^2/dbt)
    
    return(ret/2)
  }
  
  grad_Ubt <- function(x){
    Xnorm <- sqrt(colSums(X^2))/n2
    BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
    Uvarc <- BSxn %*% x#exp(x)
    
    Umat  <- (Wnorm/n2) / matrix(Xnorm, n3, n1)
    
    part1 <- rowSums((log(Umat)+matrix(Uvarc/2, n3, n1)))/Uvarc - rowSums((Umat+matrix(Uvarc/2, n3, n1))^2)/(Uvarc^2) + n1/Uvarc
    #part2 <- 2*x/100
    part1 <- array(t(BSxn) %*% part1)#*exp(x)
    ret   <- part1#+part2
    
    return(ret/2 + (x-mbt)/dbt)
  }
  
  
  
  UXX <- function(x){
    Xc     <- x
    kappac <- kappa 
    
    for(i in 1:n2){
      BSXi  <- bsplineS(Xc[i, ],  breaks=seq(A,B,(B-A)/knot), norder=norder) 
      kappac[i, ] <- BSXi %*% theta[i, ] 
    }
    
    Xnorm <- sqrt(colSums(Xc^2))/n2
    BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
    Uvarc <- BSxn %*% beta
    
    temp <- mcmapply(1:n3, FUN = meanMvMF, MoreArgs = list(kappa=kappac), mc.cores = ncore)
    
    logtempcst <- mcmapply(1:n3, FUN = logMvMFconst, MoreArgs = list(kappa=kappac), mc.cores = ncore)
    
    Wmat <- t(matrix(array(W), n1))
    
    Wmat <- Wmat * rep(array(temp), n1)
    Wc <- array(t(Wmat), dim=c(n1, n2, n3))
    
    for(i in 1:n1){
      Wnorm[, i] <- array(sqrt(colSums(Wc[i,,]^2)))
    }
    
    M1 <- matrix(array(Wmat/Wnorm[inds, ]), n2*n3, n1) #* rep(array(temp), n1)
    #M1 <- M1/Umat[inds, ]
    
    tempmat <- array(0, dim = c(n2, n3, n1))
    
    Xcn <- Xc / matrix(n2*Xnorm, n2, n3, byrow = T)
    
    for(l in 1:n1){
      tempmat[,,l] <- mcmapply(1:n3, FUN = orthtr, MoreArgs = list(k = l, mat1=M1, mat2 = array(Xcn)))
    }
    
    part1 <- n1 * sum(logtempcst)
    part2 <- - sum(array(tempmat) * rep(array(kappac), n1) )
    
    Yxc <- Yx
    
    for(i in 1:n2){
      temp.xsrs = Fx_mixtnorm(Xc[i,],probK[i, ],muK,sqrt(sigK),A,B)
      temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
      temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
      Yxc[i,] = qnorm(temp.xsrs) 
    }
    
    #Yxc <- array(qnorm(pnorm(Xc, muX, sd =  sqrt(sigX))))
    
    part4 <- t(Yxc) %*% (IRX-diag(n2)) %*% Yxc
    part3 <- (x-muX)^2 / (sigX)
    
    Umat <- (Wnorm/n2) / matrix(Xnorm, n3, n1)
    part5 <- rowSums((log(Umat)+matrix(Uvarc/2, n3, n1))^2)/Uvarc + n1*log(Uvarc)+2*rowSums(log(Umat))
    
    ret   <- sum(part1) * n1 + sum(part2) + sum(part3)/2 + sum(part4)/2 + sum(part5) / 2
    
    return(ret)
  }
  
  UMu <- function(x, j){
    indj <- which(indK==j)
    trunp <- pnorm(B, x, sd = sqrt(sigK[j])) - pnorm(A, x, sd = sqrt(sigK[j]))
    trunp[trunp<1e-20] <- 1e-20
    part1 <- sum((X[indj] - x)^2 / (2*sigK[j])) + length(indj)*log(trunp)
    
    if(length(indj)==0) {part1 = 0}
    
    ret <- part1 + ((x-muX0)^2) / (2*sigX0)
    
    return(ret)
  }
  
  grad_UMu <- function(x, j){
    indj <- which(indK==j)
    trunp <- pnorm(B, x, sd = sqrt(sigK[j])) - pnorm(A, x, sd = sqrt(sigK[j]))
    trunp[trunp<1e-20] <- 1e-20
    part <- -(dnorm(B, x, sd = sqrt(sigK[j])) - dnorm(A, x, sd = sqrt(sigK[j])))/(trunp*sqrt(sigK[j]))
    partder <- sum(- (X[indj] - x) / (sigK[j])) + length(indj)* part
    if(length(indj)==0) {partder = 0}
    ret <- partder  + (x-muX0) / (sigX0)
    
    return(ret)
  }
  
  Usig <- function(x, j){
    indj <- which(indK==j)
    trunp <- pnorm(B, muK[j], sd = sqrt(x)) - pnorm(A, muK[j], sd = sqrt(x))
    trunp[trunp<1e-20] <- 1e-20
    part1 <- - sum(dnorm(X[indj],muK[j], sd=sqrt(x), log = T)) + length(indj)*log(trunp)
    if(length(indj)==0) {part1 = 0}
    ret <- part1 + 2*log(x) + 1/x #dgamma(1/x, 1,1, log=T)
    # length(indj)*log(x)/2 + sum((X[indj] - muK[j])^2 / (2*x))
    return(ret)
  }
  
  grad_Usig <- function(x, j){
    indj <- which(indK==j)
    trunp <- pnorm(B, muK[j], sd = sqrt(x)) - pnorm(A, muK[j], sd = sqrt(x))
    trunp[trunp<1e-20] <- 1e-20
    part <- -(dnorm(B, muK[j], sd = sqrt(x))*(B-muK[j])/(2*x^(3/2)) - dnorm(A, muK[j], sd = sqrt(x))*(A-muK[j])/(2*x^(3/2)))/(trunp)
    partder <-  sum(- (X[indj] - muK[j])^2 / (2*(x)^2)) + length(indj)*part + length(indj)/(2*x)
    if(length(indj)==0) {partder = 0}
    ret <- partder + (1+1)/x - 1/(x^2)
    
    return(ret)
  }
  
  polrec <- function(angl, r =1){
    len  <- length(angl)
    temp <- exp(cumsum(log(sin(angl[1:(len-1)]))))
    temp <- c(1, temp)
    
    temp[1:len] <- temp[1:len] * cos(angl)
    temp <- c(temp, prod(sin(angl)))
    return(r * temp)
  }
  
  Fx_mixtnorm <- function(x,w,m,s,lwr,upr) 
  {
    y = matrix(0,nrow=length(w),ncol=length(x))
    for(kk in 1:length(w))
      y[kk,] = w[kk]*ptnorm(x,m[kk],s[kk],lwr,upr)
    y = colSums(y)
    return(y)
  }
  
  n1 <- dim(W)[1]
  n2 <- dim(W)[2]
  n3 <- dim(W)[3]
  
  if(is.null(Xini)){
    X <- apply(W, 2:3, mean)
    X <- X - min(X) + A
    X <- min(B,max(abs(X)))*X / max(abs(X)) 
  }
  if(!is.null(Xini)){
    X <- Xini
  }
  
  Xsd <- apply(W, 2:3, sd)
  
  Wnorm <- matrix(0, n3, n1)
  
  for(i in 1:n1){
    Wnorm[, i] <- array(sqrt(colSums(W[i,,]^2)))
  }
  
  Wmat <- t(matrix(array(W), n1))
  Xmat <- matrix(array(X), nrow = n2*n3, ncol = n1)
  
  index <- as.matrix(combinat::combn(1:n2, 2))
  e1    <- rep(0, n2)
  e1[1] <- 1
  
  # if(is.null(Thini)){
  #   theta <- array(exp(TM), dim = c(n2, knot+norder-1))#array(exp(rnorm(n2*(knot+norder-1), sd = 0.2)), dim = c(n2, knot+norder-1)) #
  #   beta  <- rep(TB, knot+norder-1)#(rnorm(knot+norder-1, sd = 0.2))^2
  # }
  theta <- matrix(0, n2, knot+norder-1)
  
  mth <- rep(0, n2)#apply(log(theta),1,mean)
  dth <- rep(0, n2)#apply(log(theta),1,var)
  
  #mbt <- mean(log(beta))
  #dbt <- var(log(beta))
  VB = 1e-1
  
  Xnorm <- sqrt(colSums(X^2))/n2
  BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
  Umat  <- (Wnorm/n2) / matrix(Xnorm, n3, n1)
  lUS   <- (apply(log(Umat), 1, sd))^2
  fit   <- nnls(BSxn, lUS)#lm(lUS~BSxn-1)#
  beta  <- array(fit$x)#fit$coefficients#
  coef <- array(fit$x)
  beta[which(is.na(beta))] <- 0
  beta[which((beta<=0))] <- 1e-6
  
  coefs <- coef[which(coef!=0)]
  mbt <- beta#0#mean(log(coefs))
  dbt <- VB#var(log(coefs))
  #beta <- rep(1e-6, length(beta))
  
  for(j in 1:n2){
    BSxj  <- bsplineS(X[j, ],  breaks=seq(A,B,(B-A)/knot), norder=norder)
    fit   <- nnls(BSxj, Thini[j, ])#lm(Thini[j, ]~BSxj-1)#
    coef <- array(fit$x)
    coefs <- coef[which(coef!=0)]
    #mth[j] <- 0#mean(log(coefs))
    dth[j] <- VB#var(log(coefs))
    #coef[which(is.na(coef))] <- 0
    #coef[which((coef<=0))] <- 1e-6
    theta[j, ] <- coef#array(exp(TM), dim = c(n2, knot+norder-1))
  }
  #theta <- array(exp(TM), dim = c(n2, knot+norder-1))#array(exp(rnorm(n2*(knot+norder-1), sd = 0.2)), dim = c(n2, knot+norder-1)) #
  #beta  <- rep(TB, knot+norder-1)#(rnorm(knot+norder-1, sd = 0.2))^2
  mth <- theta
  TM <- matrix(0, nrow(theta), ncol(theta))#TM <- (mth-50)*(mth>50)#matrix(0, nrow(mth), ncol(mth))#mth / 2
  inds  <- rep(1:n3, each = n2)
  
  Umat <- matrix(1, n3, n1)
  
  kappa <- matrix(0, n2, n3)
  for(i in 1:n2){
    BSXi <- bsplineS(X[i, ],  breaks=seq(A,B,(B-A)/knot), norder=norder)
    kappa[i, ] <- BSXi %*% theta[i, ]
  }
  
  #dth <- matrix(1e12, nrow(theta), ncol(theta))
  #dbt <- rep(1e5, knot+norder-1)
  
  Jth <- knot+norder-1
  
  #beta  <- exp(rnorm(knot+norder-1, sd = 0.2))
  Xnorm <- sqrt(colSums(X^2))/n2
  BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
  Uvar  <- BSxn %*% beta
  
  muX  <- X#matrix(0, n2, n3) #matrix(muK[indK], n2, n3)#
  sigX <- matrix(1, n2, n3) #matrix(sigK[indK], n2, n3)#
  
  if(is.null(uniout)){
    kmean.mean <- kmeans(array(X), Kx)
    
    indK <- matrix(kmean.mean$cluster, n2, n3, byrow = T)
    
    muK   <- rep(0, Kx)#
    sigK  <- rep(1, Kx)#matrix(1, n2, Kx)
    
    for(i in 1:Kx){
      muK[i]  <- mean(X[which(kmean.mean$cluster==i)])
      sigK[i] <- (sd(X[which(kmean.mean$cluster==i)]))^2
    }
    probK <- matrix(table(kmean.mean$cluster)/n3, n2, Kx, byrow = T) 
  }
  
  if(!is.null(uniout)){
    muK   <- uniout$mu
    sigK  <- uniout$sig
    probK <- uniout$prob
    indK  <- uniout$Ind
  }
  muX0  <- mean(X)
  sigX0 <- mean(apply(X[1:n2,],1,var))
  
  for(i in 1:n2){
    muX[i, ]  <- array(muK[indK[i, ]])
    sigX[i, ] <- array(sigK[indK[i, ]])
  }
  
  V     <- diag(n2)
  polth <- 0*diag(n2)
  
  Rx <- cor(t(X))
  V  <- t(chol(Rx))
  
  for(l in 2:n2){
    temp <- rect2polar(V[l, l:1])
    angl <- temp$phi %% (2*pi)
    polth[l, 1:(l-1)] <- array(angl)
  }
  
  # l = 2
  # temp              <- rep(0, l-1)#runif(l-1, 0, pi)
  # polth[l, 1:(l-1)] <- temp
  # V[l, l:1]         <- c(cos(temp), sin(temp))
  # 
  # for(l in 3:n2){
  #   temp              <- rep(0, l-1)#runif(l-1, 0, pi)
  #   polth[l, 1:(l-1)] <- temp
  #   V[l, l:1]         <- polrec(temp)
  # }
  
  indplth <- lower.tri(polth)#which(polth!=0)
  indplth <- which(indplth == T)
  indplth <- indplth
  
  polth[1,1] <- 0
  
  # Rx  <- V %*% t(V)
  IRX <- solve(Rx)
  
  Xlist   <- list()
  betal   <- list()
  thetal  <- list()
  muKl    <- list()
  varKl   <- list()
  RXl     <- list()
  IndKl   <- list()
  probKl  <- list()
  denmaSl <- list()
  Umatl   <- list()
  
  itr <- 0
  noBl <- nbl
  armu   <- rep(0, Kx)
  arsig  <- rep(0, Kx)
  arX    <- 0#matrix(0, n2, noBl)
  arth   <- rep(0, n2)
  arplth <- rep(0, length(indplth))
  arbt   <- 0
  arU    <- 0
  
  co <- rep(0, n2)
  #for(l in 1:n2){co[l] <- log(min(abs(grad_Uth(log(theta[l,]),l))))/log(10)}
  
  bco <- log(max(abs(grad_Ubt(log(beta)))))/log(10)
  
  sdmu   <- rep(1e-4, Kx)
  sdsig  <- rep(1e-8, Kx)
  sdX    <- 1e-4#matrix(1e-6, n2, noBl)#sd(Xini)/6#
  sdth   <- rep(10^(-3), n2)#rep(10^(-2*round(min(co))), n2)
  sdplth <- rep(1e-2, length(indplth))
  sdbt   <- 10^(-2*round(bco))#1e-8#
  sdU    <- 1e-3
  
  Yx <- X
  f  <- 3
  for(j in 1:n2){
    temp.xsrs = Fx_mixtnorm(X[j,],probK[j, ],muK,sqrt(sigK),A,B)
    temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
    temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
    Yx[j,] = qnorm(temp.xsrs)
  }
  
  indexup <- function(i, xl, prbkl){
    pii <- prbkl * dnorm(xl[i], muK, sd = sqrt(sigK)) / trunp# + 1e-20
    ret <- sample(1:Kx, 1, prob = pii)
    return(ret)
  }
  #grid   <- x.grid[1,]#(0:10000)/1000
  grid   <- seq(A, B,length=500)
  denmaS <- matrix(0, n2, length(grid))
  
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Mn <- min(theta)
  tempI <- 0.01*diag(Jth)
  corth <- list()
  for(l in 1:n2) {corth[[l]] <- diag(ncol(theta))}
  
  while(itr < Total_itr){
    itr <- itr + 1
    t1<-Sys.time()
    trunp <- pnorm(B, muK, sd = sqrt(sigK)) - pnorm(A, muK, sd = sqrt(sigK))
    
    trunp[trunp<1e-20] <- 1e-20
    ####HMC part
    
    temp <- HMC(Ubt, grad_Ubt, sdbt, L = 30, beta, arbt)
    beta <- temp$up
    arbt <- temp$arc
    
    Xnorm <- sqrt(colSums(X^2))/n2
    BSxn  <- bsplineS(Xnorm,  breaks=seq(A,B,(B-A)/knot), norder=norder)
    Uvar  <- BSxn %*% beta
    
    #TM <- 50
    for(k in 1:n2){
      
      BSXk  <- bsplineS(X[k, ],  breaks=seq(A,B,(B-A)/knot), norder=norder)
      if(itr > noupT){
        temp <- theta[k, ] + mvtnorm::rmvnorm(1, sigma=sdth[k]*corth[[k]])#HMCth(Uth, grad_Uth, sdth[k], L = 1, theta[k, ], k, arth[k])
        temp <- temp*(temp>TM[k,]) + TM[k,]*(temp<TM[k,])
        D <- Uth(theta[k, ], k) - Uth(temp, k)
        
        if(is.nan(D)){D=-Inf}
        if(is.na(D)){D=-Inf}
        
        if(D > log(runif(1))){
          deltath    <- temp
          theta[k, ] <- deltath
          arth[k]    <- arth[k] + 1
        }
        
        
        al <- 0.1 + Jth/2
        bl <- 0.1 +sum((theta[k, ]-mth[k])^2)/2
        dth[k] <- 1/rgamma(1, al, bl)
        
      }
      
      kappa[k, ] <- BSXk %*% theta[k, ]
      
      temp.xsrs = Fx_mixtnorm(X[k,],probK[k, ],muK,sqrt(sigK),A,B)
      temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
      temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
      Yx[k,] = qnorm(temp.xsrs)
    }
    if(itr > noupX){
      Xc <- rtnorm(n2*n3, mean = array(Xini), sd = sdX, A, B)#rtnorm(n2*n3, mean = array(X), sd = sdX, A, B) #
      Xc <- matrix(Xc, n2, n3) #matrix(rnorm(n2*n3, sd = sdX), n2, n3) #+ matrix(runif(length(indBl), -sdX[k, Bl], sdX[k, Bl]), 1, length(indBl))#
      #Xc <- Xc * (Xc>=A) * (Xc <= B) + A*(Xc<A) + B*(Xc>B)
      
      D  <- UXX(X) - UXX(Xc) + sum((Xc-Xini)^2/sdX^2) / 2 - sum((X-Xini)^2/sdX^2) / 2 #UX(X[k, ], k) - UX(Xc[k, ], k)#
      if(is.nan(D)){D=-Inf}
      if(is.na(D)){D=-Inf}
      
      if(D > log(runif(1))){
        arX <- arX + 1
        X   <- Xc
      }
    }
    # for(k in 1:Jth){
    #   
    # }
    al <- 0.1 + Jth/2
    bl <- 0.1 + sum((beta-mbt)^2)/2
    dbt <- 1/rgamma(1, al, bl)
    
    
    for(j in 1:Kx){
      temp    <- HMCthmu(UMu, grad_UMu, sdmu[j], L = 30, muK[j], j, armu[j])
      muK[j]  <- temp$up
      armu[j] <- temp$arc
      
       
      temp     <- HMCsig(Usig, grad_Usig, sdsig[j], L = 30, sigK[j], j, arsig[j])
      sigK[j]  <- temp$up
      arsig[j] <- temp$arc
    }
    
    for(j in 1:n2){
      temp.xsrs = Fx_mixtnorm(X[j,],probK[j, ],muK,sqrt(sigK),A,B)
      temp.xsrs[temp.xsrs>=1] = 0.999  # Numerical Stability
      temp.xsrs[temp.xsrs<=0] = 0.001  # Numerical Stability
      Yx[j,] = qnorm(temp.xsrs)
    }
    
    Xmat <- matrix(array(X), nrow = n2*n3, ncol = n1)
    
    j <- 0
    ####### M-H part
    for(k in indplth){
      j <- j+1
      polthc <- polth
      polthc[k] <- polth[k] + rnorm(1, sd = sdplth[j]) #runif(1, - sdplth[j], sdplth[j])# 
      
      Vc <- V
      l=2
      temp <- polthc[l, 1:(l-1)]
      #temp <- rtnorm(1, polth[l, 1:(l-1)], sd = sdplth, 0, 2*pi)
      temp <- temp*(temp>0)
      temp[1:(l-2)] <- (temp[1:(l-2)] %% pi)
      temp[l-1]     <- temp[l-1] * (temp[l-1] <= 2*pi) + 2*pi * (temp[l-1] > 2*pi)
      polthc[l, 1:(l-1)] <- temp
      Vc[l, l:1] <- c(cos(temp), sin(temp))
      
      for(l in 3:n2){
        temp <- polthc[l, 1:(l-1)]
        temp <- temp*(temp>0)
        temp[1:(l-2)] <- (temp[1:(l-2)] %% pi)
        temp[l-1]     <- (temp[l-1] %% (2 * pi))
        polthc[l, 1:(l-1)] <- temp
        Vc[l, l:1] <- polrec(temp)
      }
      
      Rxc <- Vc %*% t(Vc)
      
      D <- sum(apply(Yx, 2, mvtnorm::dmvnorm, log=T, sigma = Rxc))-sum(apply(Yx, 2, mvtnorm::dmvnorm, log=T, sigma = Rx))
      if(is.nan(D)){D=-Inf}
      if(is.na(D)){D=-Inf}
      
      if(D > log(runif(1))){
        arplth[j] <- arplth[j] + 1
        V      <- Vc
        Rx     <- Rxc
        #RXei   <- eigen(Rx)
        IRX    <- solve(Rx)#RXei$vectors %*% diag(1/(abs(RXei$values))) %*% t(RXei$vectors)
        polth  <- polthc
      }
    }
    
    
    trunp <- pnorm(B, muK, sd = sqrt(sigK)) - pnorm(A, muK, sd = sqrt(sigK))
    
    trunp[trunp<1e-20] <- 1e-20
    
    
    ####Gibbs part
    for(l in 1:n2){
      pl <- lapply(1:Kx, FUN=function(i){length(which(indK[l, ]==i))})
      probK[l, ] <- rdirichlet(1, unlist(pl)+0.001/Kx)
      
      temp  <- parallel::mcmapply(indexup, 1:n3, MoreArgs = list(xl = X[l, ], prbkl = probK[l, ]))
      indK[l, ] <- unlist(temp)
      muX[l, ]  <- muK[indK[l, ]]
      sigX[l, ] <- sigK[indK[l, ]]
    }
    
   
    
    for(l in 1:n2){
      denma <- matrix(0, Kx, length(grid))
      
      for(i in 1:Kx){
        denma[i, ] <- dnorm(grid, muK[i], sqrt(sigK[i]))
      }
      temp <- probK[l, ]
      denmaS[l, ] <- colSums(denma * matrix(temp, Kx, length(grid))) 
    }
    
    temp <- mcmapply(1:n3, FUN = meanMvMF, MoreArgs = list(kappa=kappa), mc.cores = ncore)
    
    Wmat <- t(matrix(array(W), n1))
    Wmat <- Wmat * rep(array(temp), n1)
    Wc <- array(t(Wmat), dim=c(n1, n2, n3))
    
    for(i in 1:n1){
      Wnorm[, i] <- array(sqrt(colSums(Wc[i,,]^2)))
    }
    
    if(itr %% 100 == 0){
      for(j in 1:Kx){
        C <- 0#max((grad_UMu(muK[j], j))^2)
        ar <- armu[j]/ itr
        if(ar<.60){sdmu[j] <- sdmu[j] * (.1)}
        if(ar>.90){sdmu[j] <- min(sdmu[j] * (10), 2/C)}
        
        C <- 0#max((grad_Usig(sigK[j], j))^2)
        ar <- arsig[j]/ itr
        if(ar<.60){sdsig[j] <- sdsig[j] * (.1)}
        if(ar>.90){sdsig[j] <- min(sdsig[j] * (10), 2/C)}
      }
      
      
      C <- 0#max((grad_Ubt(beta))^2)
      
      ar <- arbt/ itr
      if(ar<.60){sdbt <- sdbt * (.1)}
      if(ar>.90){sdbt <- min(sdbt * (10), 2/C)}
      if(itr > noupT){
        for(j in 1:n2){
          C <- 0#max((grad_Uth(theta[j, ], j)))#max(abs(eigen(C)$values))#
          ar <- arth[j]/ (itr-noupT)
          if(ar<.20){sdth[j] <- sdth[j] * (.1)}
          if(ar>.70){sdth[j] <- min(sdth[j] * (2), 2/C)}
        }
        
        
      }
      
      for(k in 1:length(indplth)){
        ar <- arplth[k]/ (itr-0)
        if(ar<.30){sdplth[k] <- sdplth[k] * (.5)}
        if(ar>.50){sdplth[k] <- sdplth[k] * (5)}
      }
      
      par(mfrow = c(1,n2))
      #for(i in 1:n2){plot(grid,denmaS[i,])}
      
      print(range(kappa))
      
      #print(rowSums(grid[2]*(denmaS-density.x.true)^2))
      
      cat(mean(armu)/itr, "acceptance rate for muK")
      cat(mean(arsig)/itr, "acceptance rate for sigK")
      cat(mean(arX)/(itr-noupX), "acceptance rate for X")
      cat(mean(arth)/(itr-noupT), "acceptance rate for theta")
      cat(mean(arbt)/itr, "acceptance rate for beta")
      cat(mean(arplth)/ itr, "acceptance rate for polth")
    }
    
    if(itr > burn){
      Xlist[[itr-burn]]  <- X
      thetal[[itr-burn]] <- theta
      betal[[itr-burn]]  <- beta
      muKl[[itr-burn]]   <- muK
      varKl[[itr-burn]]  <- sigK
      RXl[[itr-burn]]    <- Rx
      IndKl[[itr-burn]]  <- indK 
      denmaSl[[itr-burn]]<- denmaS
      probKl[[itr-burn]] <- probK
      Umatl[[itr-burn]]  <- Umat
    }
    
    if(itr > burn){
      if(itr %% 100 == 0){
        thetamat <- array(unlist(thetal), dim=c(nrow(theta), ncol(theta), length(thetal)))
        for(l in 1:n2){
          mat <- cor(t(thetamat[l,,]))
          if(sum(is.na(mat))){mat <- diag(ncol(mat))}
          corth[[l]] <- mat
        }
      }
    }
    t2<-Sys.time()
    #Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  out <- list(Xp = Xlist, denmaSp = denmaSl, thetap = thetal, betap = betal, Umatp = Umatl, muKp = muKl, varKp = varKl, RXp = RXl, IndKp = IndKl, probKp=probKl)
  return(out)
}