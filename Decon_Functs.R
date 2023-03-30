#############################################################


newlog <- function(x) { (log(x) > 0)*log(x) + (log(x) < 0)*((x -1) - (x-1)^2/2 + (x-1)^3/3  - (x-1)^4/4) }


#############################################################


P.mat <- function(K)
	{
	# penalty matrix for density
	D <- diag(rep(1,K))
	D <- diff(diff(D))
	P <- t(D)%*%D 
	return(P)
	}


#############################################################


B.basis <- function(x,knots)
	{
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
      		{
		act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
        
		B[act.inds,jj] <- (1/2)*(1-resc.x)^2
		B[act.inds,jj+1] <- -(resc.x^2)+resc.x+1/2
		B[act.inds,jj+2] <- (resc.x^2)/2
		}
	return(B)
	}


#############################################################


B.basis.normalized <- function(x,knots)		# equidistant knots
	{
	B <- B.basis(x,knots)
	K <- length(knots)
	delta <- knots[2]-knots[1]
	
	B[,1] <- B[,1] /(delta/6)
	B[,2] <- B[,2] /(5*delta/6)
	B[,K] <- B[,K] /(5*delta/6)
	B[,K+1] <- B[,K+1] /(delta/6)
	B[,3:(K-1)] <- B[,3:(K-1)]/delta
	
	return(B)
	}


#############################################################


B.basis.local.integral <- function(x,knots)
  {
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
    {
		act.inds <- (1:n)[(x>=knots[jj])&(x<=knots[jj+1])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
        
		B[act.inds,jj] <- 5*delta/6 + (delta/2)*(resc.x-resc.x^2+resc.x^3/3)
		B[act.inds,jj+1] <- delta/6 + delta*(-resc.x^3/3+resc.x^2/2+resc.x/2)
		B[act.inds,jj+2] <- delta*resc.x^3/6
		}
	act.inds <- (1:n)[(x>=knots[1])&(x<=knots[2])]
	B[act.inds,1] <- B[act.inds,1] - 5*delta/6 
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	act.inds <- (1:n)[(x>=knots[2])&(x<=knots[3])]
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	return(B)
  }


B.basis.cum.integral <- function(x,knots)
  {
	delta <- knots[2]-knots[1]
	n <- length(x)
	K <- length(knots)
	B <- matrix(0,n,K+1)
	for (jj in 1:(K-1))
    {
		act.inds <- (1:n)[(x>=knots[jj])]
		act.x <- x[act.inds]
		resc.x <- (act.x-knots[jj])/(knots[jj+1]-knots[jj])
    resc.x[resc.x>1] <- 1
    
		B[act.inds,jj] <- 5*delta/6 + (delta/2)*(resc.x-resc.x^2+resc.x^3/3)
		B[act.inds,jj+1] <- delta/6 + delta*(-resc.x^3/3+resc.x^2/2+resc.x/2)
		B[act.inds,jj+2] <- delta*resc.x^3/6
	  }
	act.inds <- (1:n)[(x>=knots[1])]
	B[act.inds,1] <- B[act.inds,1] - 5*delta/6 
	B[act.inds,2] <- B[act.inds,2] - delta/6 
	return(B)
  }


#############################################################


B.basis.density.coeffs <- function(thetas.density,delta)
  {
  coeffs <- exp(thetas.density)/((delta/6) * (exp(thetas.density[1])+5*exp(thetas.density[2])+6*sum(exp(thetas.density[3:(K.t-1)]))+5*exp(thetas.density[K.t])+exp(thetas.density[K.t+1])))
  return(coeffs)
  }


#############################################################


Fx.norm.B.basis <- function(x,knots,thetas.density)
  {
	delta <- knots[2]-knots[1]
	Fx <- B.basis.cum.integral(x,knots)%*%B.basis.density.coeffs(thetas.density,delta)
	return(Fx)
  }


#############################################################


rnormbspline <- function(n,coeffs,minx,maxx)
	{
	grid <- seq(minx,maxx,length=500)
	delta <- (grid[2]-grid[1])
	knots <- seq(minx,maxx,len=(length(coeffs)-1))
	y <- B.basis.normalized(grid,knots)%*%(coeffs*length(knots))
	y <- y/(sum(y)*delta)
	Fy <- cumsum(y)*delta
	u <- runif(n)	
	x <- numeric(n)
	for(ii in 1:n)
		x[ii] <- grid[min(which(u[ii]<Fy))+1]
	x <- x + runif(n,0,delta)
	return(x)
	}
	

#############################################################


dnormbspline <- function(x,coeffs,minx,maxx)	# x must be equidistant grid
	{
	knots <- seq(minx,maxx,len=(length(coeffs)-1))
	y <- B.basis.normalized(x,knots)%*%(coeffs*length(knots))
	delta <- (x[2]-x[1])
	y <- y/(sum(y)*delta)
	return(y)
	}
	

#############################################################


rskewnorm <- function(n,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	xi = - (sqrt(c) * skewness)/sqrt(1 + skewness^2 * (1-c))
	omega = sqrt(1+xi^2)

	y = delta * abs(rnorm(n,0,1)) + (1-delta^2)^(1/2) * rnorm(n,0,1)
	y = xi + omega * y
	return(mean+sd*y)
	}


dskewnorm <- function(x,mean,sd,skewness)
	{
	c = 2/3.1415926
	delta = skewness/sqrt(1+skewness^2)
	zeta1 = delta * sqrt(c)
	zeta2 = sqrt(1 - c*delta^2)
	y = numeric(length(x))
	xmod = zeta1 + zeta2*(x-mean)/sd
	for(i in 1:length(x))
      	y[i] = (2*zeta2/sd[i]) * dnorm(xmod[i]) * pnorm(skewness*xmod[i])
	return(y)
	}


#############################################################


d.scaled.restricted.mix.norm <- function(x,mean,sd,pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sd.e = sqrt(var.e.fn(pi,params))
	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    	y[kk,] = pi[kk]*(p[kk]*dnorm(x,(mean+sd*mu1[kk])/sd.e,sd*sqrt(sigmasq1[kk])/sd.e) + (1-p[kk])*dnorm(x,(mean+sd*mu2[kk])/sd.e,sd*sqrt(sigmasq2[kk])/sd.e))
	y = colSums(y)
	return(y)
	}



d.restricted.mix.norm <- function(x,mean,sd,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

  y = p*dnorm(x,mean+sd*mu1,sd*sqrt(sigmasq1)) + (1-p)*dnorm(x,mean+sd*mu2,sd*sqrt(sigmasq2))
	return(y)
	}



mixnorm <- function(n,pi,p,mu_curl,sigmasq1,sigmasq2,e.grid,plot=TRUE)
	{
	m = length(pi)
	y = numeric(n)
	density <- numeric(length(e.grid))
	inds = sample(1:m,n,TRUE,prob=pi)
	for(ii in 1:m)
		{
		temp = which(inds==ii)
		y[temp] = r.restricted.mix.norm(length(temp),c(p[ii],mu_curl[ii],sigmasq1[ii],sigmasq2[ii]))
		}
	params = cbind(p,mu_curl,sigmasq1,sigmasq2)
	y = y/sqrt(var.e.fn(pi,params))
	density <- d.scaled.restricted.mix.norm(e.grid,mean=0,sd=1,pi,params)
	if(plot)
		{
		hist(y,xlim=c(min(e.grid),max(e.grid)),breaks=30,freq=FALSE)
		points(e.grid,density,type="l")
		}

	return(list(es=y,density=density))
	}


r.restricted.mix.norm <- function(n,params)
	{
	p = params[1]
	mu_curl = params[2]
	sigmasq1 = params[3]
	sigmasq2 = params[4]

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)

  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

  inds = sample(0:1,n,TRUE,prob=c(p,(1-p)))
	y = numeric(n)
	temp = which(inds==0)
	y[temp] = rnorm(length(temp),mu1,sqrt(sigmasq1)) 
	temp = which(inds==1)
	y[temp] = rnorm(length(temp),mu2,sqrt(sigmasq2))
	return(y)
	}



var.e.fn <- function(pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)


  mu1 = c1*mu_curl
	mu2 = c2*mu_curl


  y = p*(mu1^2+sigmasq1) + (1-p)*(mu2^2+sigmasq2)
	y = sum(pi*y)
	return(y)
	}



gamma.e.fn <- function(pi,params)
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)


  mu1 = c1*mu_curl
	mu2 = c2*mu_curl


  m2 = p*(mu1^{2}+sigmasq1) + (1-p)*(mu2^{2}+sigmasq2)
	m2 = sum(pi*m2)

  m3 = p*(mu1^{3}+3*mu1*sigmasq1) + (1-p)*(mu2^{3}+3*mu2*sigmasq2)
	m3 = sum(pi*m3)

  m4 = p*(mu1^{4}+6*mu1^{2}*sigmasq1+3*sigmasq1^{2}) + (1-p)*(mu2^{3}+6*mu2^{2}*sigmasq2+3*sigmasq2^{2})
	m4 = sum(pi*m4)

	gamma1 = m3/m2^{3/2}
	gamma2 = m4/m2^{2}-3

	return(c(gamma1,gamma2))
	}


r.proposal.params.restricted.mix.norm <- function(p.a,p.b,sigmasq.mu_curl,s11,s12,s21,s22)
	{
	p = rbeta(1,p.a,p.b)

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
	mu_curl = rnorm(1,0,sqrt(sigmasq.mu_curl))
	sigmasq1 = 1/rgamma(1,s11,s12)
	sigmasq2 = 1/rgamma(1,s21,s22)

	y = c(p,mu_curl,sigmasq1,sigmasq2)
	return(y)
	}    


r.tnorm.proposal.params.restricted.mix.norm <- function(params)
	{
	current.p <- params[1]
	current.mu_curl <- params[2]
	current.sigmasq1 <- params[3]
	current.sigmasq2 <- params[4]

	p = rtnorm(1,current.p,0.01,lower=0,upper=1)
	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
    
	sigmasq1 = rtnorm(1,current.sigmasq1,0.1,lower=0,upper=Inf)
	sigmasq2 = rtnorm(1,current.sigmasq2,0.1,lower=0,upper=Inf)

  mu_curl = rnorm(1,current.mu_curl,0.1)

	y = c(p,mu_curl,sigmasq1,sigmasq2)
	return(y)
	}    




#############################################################

dlaplace <- function(x,mu,b)
	{
	y = exp(-abs(x-mu)/b)/(2*b)
	return(y)
	}

plaplace <- function(x,mu,b)
	{
	y = numeric(length(x))
	for(ii in 1:length(y))
		{
		if(x[ii]<mu) 
			y[ii] = (1/2)*exp(-abs(x[ii]-mu)/b)
		else
			y[ii] = 1-(1/2)*exp(-abs(x[ii]-mu)/b)
		}		
	return(y)
	}


rlaplace <- function(n,mu,b)
	{
	u = runif(n)
	y = mu - b*sign(u-0.5)*log(1-2*abs(u-0.5))
	return(y)
	}



mixlaplace <- function(n,pi,mu,b,e.grid,plot=TRUE)		# Produces scaled mixtures of Laplace
	{
	m = length(pi)
	y = numeric(n)
	density = numeric(length(e.grid))
	inds = sample(1:m,n,TRUE,pi)
    mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(ii in 1:m)
		{
		temp = which(inds==ii)
		y[temp] = rlaplace(length(temp),mu[ii],b[ii])
		density = density + pi[ii]*dlaplace(e.grid,(mu[ii]-mu.e)/sd.e,b[ii]/sd.e)
		}
	y = (y-mu.e)/sd.e
	if(plot)
		{
		hist(y,xlim=c(min(e.grid),max(e.grid)),freq=FALSE,breaks=30)
		points(e.grid,density,type="l")
		}
	return(list(es=y,density=density))
	}


d.scaled.restricted.mixlaplace <- function(x,pi,mu,b)		# Produces scaled mixtures of Laplace
	{
	y = matrix(0,nrow=length(pi),ncol=length(x))
    mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*dlaplace(x,(mu[kk]-mu.e)/sd.e,b[kk]/sd.e)
	y = colSums(y)
	return(y)
	}



gamma.e.Laplace.fn <- function(pi,mu,b)
	{
	m1.dash = sum(pi*mu) 
	m2.dash = sum(pi*(2*b^{2}+mu^{2}))
	m3.dash = sum(pi*(mu^{3}+6*b^{2}*mu))
	m4.dash = sum(pi*(mu^{4}+12*b^{2}*mu^{2}+24*b^{4}))

	m2 = m2.dash-m1.dash^{2}
  m3 = m3.dash - 3*m2.dash*m1.dash + 2*m1.dash^{3}
	m4 = m4.dash - 4*m3.dash*m1.dash + 6*m2.dash*m1.dash^{2} - 3*m1.dash^{4}

  gamma1 = m3/m2^{3/2}
	gamma2 = m4/m2^{2}-3

	return(c(gamma1,gamma2))
	}



#############################################################


fr <- function(thetas,xs,mis,knots.t,P.t,s2t,us)			# function to be maximized
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum(rep(log(vars),times=mis))/2 - sum(us^2/rep(vars,times=mis))/2
	return(-y)
  }

fr.density <- function(thetas,xs,knots.t,P.t,s2t,density.est)			# function to be maximized
	{
  delta = knots.t[2]-knots.t[1]
	density.est2 = B.basis(xs,knots.t)%*%B.basis.density.coeffs(thetas,delta)
	y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum((density.est-density.est2)^2)
	return(-y)
	}

fr.prob.consumption <- function(thetas,xs,knots.t,P.t,s2t,probs.emp)			# function to be maximized
{
  probs.est = as.vector(pnorm(B.basis(xs,knots.t)%*%thetas))
  y = - t(thetas) %*% P.t %*% thetas / (2*s2t) - sum((probs.emp-probs.est)^2)
  return(-y)
}


############################################################


gr <- function(thetas,xs,mis,knots.t,P.t,s2t,us)			# gradient function of fr
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)
	y = - P.t %*% thetas / s2t
	B.basis.components = matrix(0,nrow=K.t+1,ncol=n)
	for(kk in 1:(K.t+1))
		{
   		thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
		}
	for(kk in 1:(K.t+1))
		for(ii in 1:n)
			for(jj in (sum(mis[1:ii-1])+1):sum(mis[1:ii]))
				y[kk] = y[kk] - (1 - (us[jj]^2)/vars[ii]) * B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
	return(-y)
	}


############################################################


prop.sig.thetas.fn <- function(thetas,xs,mis,us,s2t,K.t,P.t,knots.t,n)
	{
	vars = B.basis(xs,knots.t)%*%exp(thetas)

	B.basis.components = matrix(0,nrow=K.t+1,ncol=n);
	for(kk in 1:(K.t+1))
		{
   	thetas.new = rep(0,K.t+1)
		thetas.new[kk] = 1
		B.basis.components[kk,] = B.basis(xs,knots.t)%*%thetas.new
	  }

	prop.sig.thetas = matrix(0,nrow=K.t+1,ncol=K.t+1)
	for(kk in 1:(K.t+1))
		{
		for(ll in kk:(K.t+1))
			{
			if(kk==ll)
			  {
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:(ii-1)])+1):sum(mis[1:ii]))
	        				prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2) + (1-(us[jj]^2)/vars[ii])*B.basis.components[kk,ii]*exp(thetas[kk])/(2*vars[ii])
			  }
			else
			  {
	  			for(ii in 1:n)
	     				for(jj in (sum(mis[1:(ii-1)])+1):sum(mis[1:ii]))
	       					prop.sig.thetas[kk,ll] = prop.sig.thetas[kk,ll] + ((us[jj]^2)/vars[ii] - 1/2)*(B.basis.components[kk,ii]*B.basis.components[ll,ii])*exp(thetas[kk]+thetas[ll])/(vars[ii]^2)
			  }
		  } 
		}
	for(kk in 2:(K.t+1))
   		for(ll in 1:(kk-1))
			    prop.sig.thetas[kk,ll] = prop.sig.thetas[ll,kk]
	prop.sig.thetas = prop.sig.thetas + P.t/s2t
	prop.sig.thetas = round(solve(prop.sig.thetas),4)
	
	return(prop.sig.thetas)
	}


#############################################################


pi.fn <- function(x,alpha,beta,xstar,K)
	{
	w <- pnorm(alpha-beta*abs(x-xstar)^2)
	pi <- numeric(K)
	pi[1] <- w[1]
	for(k in 2:(K-1))
		pi[k] <- w[k]*prod(1-w[1:(k-1)])
	pi[K] <- prod(w[1:(K-1)])
	return(pi)
	}


#############################################################
# Folded Normal Distribution 

dfldnorm <- function(x,mu,sd) 
	{
	dens = ifelse(x<0,0,dnorm(-x,mu,sd) + dnorm(x,mu,sd))
	return(dens)
	}



#############################################################


fx_mixnorm <- function(x,w,m,s) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*dnorm(x,m[kk],s[kk])
	y = colSums(y)
	return(y)
	}
	
Fx_mixnorm <- function(x,w,m,s) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*pnorm(x,m[kk],s[kk])
	y = colSums(y)
	return(y)
	}
	
Fx_mixnorm_inv <- function(p,w,m,s,br=c(-1000,1000))	# Currently accepts only scalar p
	{
	G <- function(x) {Fx_mixnorm(x,w,m,s) - p}
	return(uniroot(G,br)$root) 
	}
	
	
	
fx_mixtnorm <- function(x,w,m,s,lwr,upr) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*dtnorm(x,m[kk],s[kk],lwr,upr)
	y = colSums(y)
	return(y)
	}
	
Fx_mixtnorm <- function(x,w,m,s,lwr,upr) 
	{
	y = matrix(0,nrow=length(w),ncol=length(x))
	for(kk in 1:length(w))
	    y[kk,] = w[kk]*ptnorm(x,m[kk],s[kk],lwr,upr)
	y = colSums(y)
	return(y)
	}
	
Fx_mixtnorm_inv <- function(p,w,m,s,lwr,upr,br=c(-1000,1000))	# Currently accepts only scalar p
	{
	G <- function(x) {Fx_mixtnorm(x,w,m,s,lwr,upr) - p}
	return(uniroot(G,br)$root) 
	}
	
	

Fe_scaled_mixnorm <- function(x,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	sd.e = sqrt(var.e.fn(pi,params))
	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    y[kk,] = pi[kk]*(p[kk]*pnorm(x,mu1[kk]/sd.e,sqrt(sigmasq1[kk])/sd.e) + (1-p[kk])*pnorm(x,mu2[kk]/sd.e,sqrt(sigmasq2[kk])/sd.e))
	y = colSums(y)
	return(y)
	}
		
Fe_scaled_mixnorm_inv <- function(p,pi,params,br=c(-100,100))	# Currently accepts only scalar p
	{
	G <- function(x) {Fe_scaled_mixnorm(x,pi,params) - p}
	return(uniroot(G,br)$root) 
	}


Fe_scaled_laplace <- function(x,pi,mu,b) 
	{
	y = matrix(0,nrow=length(pi),ncol=length(x))
  mu.e = sum(pi*mu)
	sd.e = sqrt(sum(pi*(2*b^{2}+mu^{2}))-mu.e^{2})
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*plaplace(x,(mu[kk]-mu.e)/sd.e,b[kk]/sd.e)
	y = colSums(y)
	return(y)
	}
		
Fe_scaled_laplace_inv <- function(p,pi,mu,b,br=c(-100,100))
	{
	G <- function(x) {Fe_scaled_laplace(x,pi,mu,b) - p}
	return(uniroot(G,br)$root) 
	}


fu_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	  y[kk,] = pi[kk]*(p[kk]*dnorm(x,mean+sd*mu1[kk],sd*sqrt(sigmasq1[kk])) + (1-p[kk])*dnorm(x,mean+sd*mu2[kk],sd*sqrt(sigmasq2[kk])))
	y = colSums(y)
	return(y)
	}


Fu_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		p = params[,1]
		mu_curl = params[,2]
		sigmasq1 = params[,3]
		sigmasq2 = params[,4]
		}
	else
		{
		p = params[1]
		mu_curl = params[2]
		sigmasq1 = params[3]
		sigmasq2 = params[4]
		}

	c1 = (1-p)/(p^2+(1-p)^2)
	c2 = -p/(p^2+(1-p)^2)
    
  mu1 = c1*mu_curl
	mu2 = c2*mu_curl

	y = matrix(0,nrow=length(p),ncol=length(x))
	for(kk in 1:length(p))
	    y[kk,] = pi[kk]*(p[kk]*pnorm(x,mean+sd*mu1[kk],sd*sqrt(sigmasq1[kk])) + (1-p[kk])*pnorm(x,mean+sd*mu2[kk],sd*sqrt(sigmasq2[kk])))
	y = colSums(y)
	return(y)
	}


fu_eps_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		mu = params[,1]
		sigmasq = params[,2]
		}
	else
		{
		mu = params[1]
		sigmasq = params[2]
		}
    
	y = matrix(0,nrow=length(pi),ncol=length(x))
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*dnorm(x,mean+sd*mu[kk],sd*sqrt(sigmasq[kk])) 
	y = colSums(y)
	return(y)
	}


Fu_eps_mixnorm <- function(x,mean,sd,pi,params) 
	{
	if(is.matrix(params))
		{
		mu = params[,1]
		sigmasq = params[,2]
		}
	else
		{
		mu = params[1]
		sigmasq = params[2]
		}
    
	y = matrix(0,nrow=length(pi),ncol=length(x))
	for(kk in 1:length(pi))
	    y[kk,] = pi[kk]*pnorm(x,mean+sd*mu[kk],sd*sqrt(sigmasq[kk]))
	y = colSums(y)
	return(y)
	}



#############################################################


formCorr.xs <- function(b,theta)
	{
	p = length(b)+1	
	V = matrix(0,p,p)
	V[1,1] = 1
	V[2,1] = b[1]
	V[2,2] = sqrt(1-b[1]*b[1])
	for(jj in 3:p)
		{
		q1jj = 1+(jj-3)*(jj-2)/2
		q2jj = (jj-2)*(jj-1)/2
		V[jj,1] = b[jj-1]*sin(theta[q1jj])
		if(jj>3)
			for(kk in 2:(jj-2))
				V[jj,kk] = b[jj-1]*prod(cos(theta[q1jj:(q1jj+kk-2)]))*sin(theta[q1jj+kk-1])
		V[jj,jj-1] = b[jj-1]*prod(cos(theta[q1jj:q2jj]))	
		V[jj,jj] = sqrt(1-b[jj-1]*b[jj-1])	
		}
	R = V%*%t(V)
	return(R)	
	}
	

formCorr.es <- function(b,theta)
	{
	p = length(b)+1	
	V = matrix(0,p,p)
	V[1,1] = 1
	V[2,1] = b[1]
	V[2,2] = sqrt(1-b[1]*b[1])
	for(jj in 3:p)
		{
		q1jj = 1+(jj-3)*(jj-2)/2
		q2jj = (jj-2)*(jj-1)/2
		V[jj,1] = b[jj-1]*sin(theta[q1jj])
		if(jj>3)
			for(kk in 2:(jj-2))
				V[jj,kk] = b[jj-1]*prod(cos(theta[q1jj:(q1jj+kk-2)]))*sin(theta[q1jj+kk-1])
		V[jj,jj-1] = b[jj-1]*prod(cos(theta[q1jj:q2jj]))	
		V[jj,jj] = sqrt(1-b[jj-1]*b[jj-1])	
		}
	R = V%*%t(V)
	return(R)	
	}


































#########################################
### Bayesian Univariate Deconvolution ###
#########################################

# The function is for univariate density deconvolution for variables with strictly continuously measured surrogates 
# as described in the paper "Bayesian Copula Density Deconvolution for Zero-Inflated Data with Applications in Nutritional Epidemiology" by Sarkar, Pati, Mallick and Carroll.

# The method uses mixtures of truncated-normals with shared atoms to model the density of interest, 
# mixtures of moment-restricted normals to model the density of the measurement errors,
# and mixtures of B-spines to model the conditional variability of the measurement errors.
# See paper for additional details.



#############
### Input ###
#############

# While running from within the file 'Bayes_Copula_Decon_MVT.R' that implements the multivariate method, these arguments are read from the original file.
# The univariate method can also be independently implemented using the current file.  

# ws <- strictly continuously measured surrogate values
# xs.lwr <- lower limit of the range of the variable of interest
# xs.upr <- upper limit of the range of the variable of interest
# mis <- no of surrogates for each subject, must be greater than or equal to 3
# z.xs.max <- number of mixture components allowed in the model for the density of interest
# z.us.max <- number of mixture components allowed in the model for the density of the measurement errors
# K.t <- number of B-spline knots for the variance functions modeling conditional variability of the measurement errors
# simsize <- total num of MCMC iterations
# burnin <- burnin for the MCMC iterations
# show_progress <- if TRUE, shows progress by printing every 100th iteartion number, MUST be set at FALSE while running in parrellel from within 'Bayes_Copula_Decon_MVT.R'
# plot_results <- if TRUE, plots the estimated density of interest, the estimated density of measurement errors, the estimated variance function etc., MUST be set at FALSE while running in parrellel from within 'Bayes_Copula_Decon_MVT.R'



##############
### Output ###
##############

# Output comprises a list of the following variables. 
# While running from within the file 'Bayes_Copula_Decon_MVT.R' that implements the multivariate method, these variables are used as.

# knots <- knot-points for constructing the B-splines bases that model the conditional variability of the measurement errors 
# thetas <- estimated coefficients of B-splines bases that model the conditional variability of the measurement errors 
# xs <- estimated subject-specific values of the variable of interest
# us <-  estimated subject and replicate-specific values of the measurement errors
# z.xs <- mixture component labels for the mixture model for the density of interest
# pi.xs <- mixture component probabilities for the mixture model for the density of interest
# params.xs <- mixture component-specific parameters for the mixture model for the density of interest
# z.us <- mixture component labels for the mixture model for the density of the measurement errors
# pi.us <- mixture component probabilities for the mixture model for the density of the measurement errors
# params.us <- mixture component-specific parameters for the mixture model for the density of the measurement errors



UNIV_DECON_REGULAR = function(ws, xs.lwr, xs.upr, mis, z.xs.max, z.us.max, K.t, simsize, burnin, show_progress=TRUE, plot_results=TRUE)
{
  #################################
  ### Priors and Initial Values ###
  #################################
  
  ### Initialization and prior of xs and us
  n = length(mis)
  N = sum(mis)
  inds = rep(1:n,times=mis)
  inds1 = inds2 = numeric(n)
  inds1[1] = 1
  inds2[1] = inds1[1]+mis[1]-1
  for(ii in 2:n)
  {
    inds1[ii] = inds1[ii-1]+mis[ii-1]
    inds2[ii] = inds1[ii] + mis[ii]-1
  }
  wbars = tapply(ws,inds,"mean")
  xs = as.vector(wbars)
  xs[xs <= xs.lwr] = xs.lwr+0.1
  xs[xs >= xs.upr] = xs.upr-0.1
  current.xs = start.xs = xs
  us = ws - rep(xs,times=mis)
  range.start.xs = diff(range(xs))
  s2is = as.vector(tapply(ws,inds,var))
  xs.grid = seq(xs.lwr,xs.upr,length=500)
  xs.grid.length = length(xs.grid)
  
  alpha.xs = 1
  
  # Normal 
  mu0.xs = mean(xs)
  sigmasq0.xs = var(xs)
  
  
  # Inverse-Gamma (Independnt from mu - independence is important)
  a.sigmasq.xs = 1
  b.sigmasq.xs = 1
  
  pi.xs = rep(1/z.xs.max,z.xs.max)
  clusters.xs = kmeans(xs,z.xs.max) 
  mu.xs = clusters.xs$center
  z.xs = clusters.xs$cluster
  sigmasq.xs = rep(var(xs)/5,z.xs.max)
  
  d.ordinates.xs = matrix(0,nrow=n,ncol=z.xs.max)
  
  ### Prior and initialization of s2t and thetas
  alpha.t = 100
  beta.t = 1
  s2t = 0.01
  P.t = P.mat(K.t+1) 	# penalty matrix
  knots.t = seq(xs.lwr,xs.upr,length=K.t)
  
  optim_results = optim(rep(1,K.t+1), fr, NULL, xs, mis, knots.t, P.t, s2t, us, method = "BFGS")
  thetas = current.thetas = start.thetas = optim_results$par
  prop.sig.thetas = make.positive.definite(prop.sig.thetas.fn(thetas,xs,mis,us,s2t,K.t,P.t,knots.t,n))
  
  var.grid = seq(xs.lwr,xs.upr,length=100)
  vars = current.vars = t(B.basis(xs,knots.t)%*%exp(current.thetas))
  B.basis.var.grid.knots.t = B.basis(var.grid,knots.t)
  B.basis.store = B.basis(xs.grid,knots.t)
  close.ind = rep(0,n)
  
  ### Prior and initialization for mixture
  simsize.mh.us = 10
  z.us = rep(1,N)
  alpha.us = 0.1
  params.us = matrix(c(0.5,0,1,1),nrow=z.us.max,ncol=4,byrow=T)		# unique values
  pi.us = rep(1/z.us.max,z.us.max)
  d.ordinates.us = matrix(0,nrow=N,ncol=z.us.max)
  
  
  
  #########################
  ### Tuning Parameters ###
  #########################
  
  
  sig.tune.thetas.1 = 0
  sig.tune.thetas.2 = 0.1
  
  
  
  
  
  ###############################
  ### Storage for MCMC Output ###
  ###############################
  
  
  es.grid = seq(-3,3,length=500)
  density.xs.est = numeric(xs.grid.length)
  var.es = numeric(1)
  var.est = numeric(length(var.grid))
  density.es.est = numeric(length(es.grid))
  prob.consumption.est = numeric(xs.grid.length)
  
  proposed.xs = current.xs = xs 
  proposed.us = current.us = us
  proposed.vars = current.vars = vars 
  
  current.likelihood = proposed.likelihood = matrix(1,2,n)
  temp.proposed.us.likelihood = temp.current.us.likelihood = matrix(1,2,N)
  
  thetas.est  = numeric(length(thetas))
  thetas.MCMC = matrix(0,nrow=simsize,ncol=length(thetas))
  
  
  
  ##################
  ### Start MCMC ###
  ##################
  
  
  for (iii in 1:simsize)
  {
    if((show_progress==TRUE)&&(iii%%10==0))
      print(iii)
    
    
    ### Updating z.xs
    for(kk in 1:z.xs.max)
      d.ordinates.xs[,kk] = dtnorm(xs,mu.xs[kk],sqrt(sigmasq.xs[kk]),lower=xs.lwr,upper=xs.upr)
    d.ordinates.xs[is.nan(d.ordinates.xs)] = 0
    d.ordinates.xs[is.infinite(d.ordinates.xs)] = 0
    for(ii in 1:n)
      z.xs[ii] = sample(z.xs.max,1,prob=pi.xs*d.ordinates.xs[ii,])
    
    
    ### Updating cluster probabilities
    n.kk.xs = tabulate(z.xs,nbins=z.xs.max)
    pi.xs = rdirichlet(1,alpha.xs/z.xs.max+n.kk.xs)	
    
    
    ### Updating mu.xs, sigmasq.xs
    
    xs.trans = mu.xs[z.xs]+sqrt(sigmasq.xs[z.xs])*qnorm((pnorm((xs-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))-pnorm((xs.lwr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs])))/(pnorm((xs.upr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))-pnorm((xs.lwr-mu.xs[z.xs])/sqrt(sigmasq.xs[z.xs]))))
    xs.trans[xs.trans < xs.lwr - 10] = xs.lwr - 10
    xs.trans[xs.trans > xs.upr + 10] = xs.upr + 10
    
    for(kk in 1:z.xs.max)
    {
      temp = which(z.xs==kk)
      xspool = xs.trans[temp]
      
      sigmasq.temp = 1/(n.kk.xs[kk]/sigmasq.xs[kk] + 1/sigmasq0.xs)
      mu.temp = (sum(xspool)/sigmasq.xs[kk] + mu0.xs/sigmasq0.xs) * sigmasq.temp
      mu.xs[kk] = rnorm(1,mu.temp,sqrt(sigmasq.temp))	
      
      post.a.sigmasq.xs = a.sigmasq.xs + length(xspool)/2
      post.b.sigmasq.xs = b.sigmasq.xs + sum((xspool-mu.xs[kk])^2)/2
      sigmasq.xs[kk] = 1/rgamma(1,shape=post.a.sigmasq.xs,rate=post.b.sigmasq.xs)
    }
    
    
    
    ### Updating xs (and us)
    
    proposed.xs = rtnorm(n,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)	
    TempMat = abs(matrix(rep(proposed.xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(thetas)
    
    proposed.prior = dtnorm(proposed.xs,mu.xs[z.xs],sqrt(sigmasq.xs[z.xs]),lower=xs.lwr,upper=xs.upr)
    current.prior = dtnorm(current.xs,mu.xs[z.xs],sqrt(sigmasq.xs[z.xs]),lower=xs.lwr,upper=xs.upr)
    
    proposed.us = ws-rep(proposed.xs,times=mis)
    
    k.us = max(z.us)	
    temp.current.us.likelihood = fu_mixnorm(current.us,mean=0,sd=rep(sqrt(current.vars),times=mis),pi.us[1:k.us],params.us[1:k.us,])
    temp.proposed.us.likelihood = fu_mixnorm(proposed.us,mean=0,sd=rep(sqrt(proposed.vars),times=mis),pi.us[1:k.us],params.us[1:k.us,])
    current.likelihood = tapply(temp.current.us.likelihood,inds,"prod")
    proposed.likelihood = tapply(temp.proposed.us.likelihood,inds,"prod")
    
    mh.ratio = (proposed.prior * proposed.likelihood * dtnorm(current.xs,mean=proposed.xs,sd=0.1,lower=xs.lwr,upper=xs.upr)) / (current.prior * current.likelihood * dtnorm(proposed.xs,mean=current.xs,sd=0.1,lower=xs.lwr,upper=xs.upr))
    
    mh.ratio[is.nan(mh.ratio)] = 0
    
    u = runif(n)
    inds.to.replace = (1:n)[u<mh.ratio]
    xs[inds.to.replace] = current.xs[inds.to.replace] = proposed.xs[inds.to.replace]
    vars[inds.to.replace] = current.vars[inds.to.replace] = proposed.vars[inds.to.replace]
    
    us = current.us = ws - rep(xs,times=mis)
    
    
    
    ### Updating thetas
    
    proposed.thetas = t(rmvnorm(1,current.thetas,(diag(rep(sig.tune.thetas.1,(K.t+1)))+sig.tune.thetas.2*prop.sig.thetas)))
    TempMat = abs(matrix(rep(xs,xs.grid.length),n,xs.grid.length)-matrix(rep(xs.grid,n),n,xs.grid.length,byrow=T))
    close.ind = apply(TempMat,1,which.min) 
    proposed.vars = B.basis.store[close.ind,]%*%exp(proposed.thetas)
    
    current.log.prior = - t(current.thetas)%*%P.t%*%current.thetas/(2*s2t)
    proposed.log.prior = - t(proposed.thetas)%*%P.t%*%proposed.thetas/(2*s2t)
    
    temp.current.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=rep(sqrt(current.vars),times=mis),params.us[z.us,])
    temp.proposed.likelihood = d.restricted.mix.norm(us,mean=rep(0,times=N),sd=rep(sqrt(proposed.vars),times=mis),params.us[z.us,])
    
    current.log.likelihood = sum(log(temp.current.likelihood))
    proposed.log.likelihood = sum(log(temp.proposed.likelihood))
    
    log.mh.ratio = proposed.log.prior + proposed.log.likelihood - current.log.likelihood - current.log.prior
    if(is.nan(log.mh.ratio)) log.mh.ratio = -Inf
    if(log(runif(1))<log.mh.ratio)
    {
      thetas = current.thetas = proposed.thetas
      vars = current.vars = proposed.vars
    }
    
    ### Updating s2t
    
    s2t = 1/rgamma(1,shape=alpha.t+(K.t+1)/2,rate=beta.t+t(thetas)%*%P.t%*%thetas)
    
    
    ### Updating z.us
    
    for(ii in 1:N)
    {			
      prob.us = pi.us * d.restricted.mix.norm(us[ii],mean=0,sd=sqrt(vars[inds[ii]]), params.us)
      if(sum(prob.us)==0) {prob.us=rep(1/z.us.max,z.us.max)}
      z.us[ii] = sample(1:z.us.max,1,TRUE,prob.us)   # New z.us[ii] drawn
    }
    
    ### Updating cluster probabilities
    
    n.kk.us = tabulate(z.us,nbins=z.us.max)
    pi.us = rdirichlet(1,alpha.us/z.us.max+n.kk.us)
    
    ### Updating params.us
    
    k.us = max(z.us)                # Number of clusters
    if(iii>2000) simsize.mh.us = 1
    for(rr in 1:simsize.mh.us)
    {
      for(kk in 1:k.us)
      {
        temp = which(z.us==kk)
        uspool = us[temp]
        varspool = vars[inds[temp]]
        
        proposed.params.us = r.tnorm.proposal.params.restricted.mix.norm(params.us[kk,])
        temp.proposed.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),proposed.params.us))
        temp.current.log.likelihood = log(d.restricted.mix.norm(uspool,mean=0,sd=sqrt(varspool),params.us[kk,]))
        temp.proposed.log.likelihood[is.infinite(temp.proposed.log.likelihood)] = 0
        temp.current.log.likelihood[is.infinite(temp.current.log.likelihood)] = 0
        proposed.log.likelihood = sum(temp.proposed.log.likelihood)
        current.log.likelihood = sum(temp.current.log.likelihood)
        
        log.acc.prob = proposed.log.likelihood-current.log.likelihood
        if(log(runif(1))<log.acc.prob)
          params.us[kk,] = proposed.params.us
      }
    }
    if(k.us<z.us.max)
      for(kk in (k.us+1):z.us.max)
        params.us[kk,] = r.proposal.params.restricted.mix.norm(1,1,3,3,3,3,3)	
    
    var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
    params.us[1:k.us,2] = params.us[1:k.us,2]/sqrt(var.es)
    params.us[1:k.us,3:4] = params.us[1:k.us,3:4]/var.es
    
    
    
    thetas.MCMC[iii,] = thetas
    
    if(iii>burnin)
    {		
      for(kk in 1:z.xs.max)
        density.xs.est = density.xs.est + pi.xs[kk]*dtnorm(xs.grid,mu.xs[kk],sqrt(sigmasq.xs[kk]),lower=xs.lwr,upper=xs.upr)
      
      k.us = max(z.us)
      var.es = var.e.fn(pi.us[1:k.us],params.us[1:k.us,])
      density.es.est = density.es.est + d.scaled.restricted.mix.norm(es.grid,0,1,pi.us[1:k.us],params.us[1:k.us,])
      var.est = var.est + B.basis.var.grid.knots.t %*% exp(thetas) * var.es
      thetas.est = thetas.est + log(var.es) + thetas
    }
    
  }
  
  
  
  density.xs.est = density.xs.est/(simsize-burnin)
  density.es.est = density.es.est/(simsize-burnin)
  var.est = var.est/(simsize-burnin)
  thetas.est = thetas.est/(simsize-burnin)
  
  thetas.final = thetas.est
  xs.final = xs
  var.final = sqrt(B.basis(xs.final,knots.t)%*%exp(thetas.final))
  us.final = (ws-rep(xs.final,times=mis))
  
  
  if(plot_results==TRUE)
  {
    dev.new()
    par(mfrow=c(2,2))
    plot(xs.grid,density.xs.est,xlab="x",ylab="f(x)",type="l",lty=1,col="green3",lwd=3)
    
    plot(es.grid,density.es.est,xlab="e",ylab="f(e)",type="l",lty=1,col="green3",lwd=3)
    points(es.grid,dnorm(es.grid),type="l",lty=1)
    
    plot(xs,s2is,pch="*",xlab="x",ylab="v(x)")
    points(var.grid,var.est,type="l",lty=1,col="blue",lwd=2)
    points(var.grid,B.basis.var.grid.knots.t%*%exp(thetas.final),type="l",lty=1,col="green3",lwd=2)
    par(mfrow=c(1,1))
  }
  
  
  params.xs = cbind(mu.xs,sigmasq.xs)
  return(list(knots=knots.t, thetas=thetas.final, 
              xs=xs.final, us=us.final, 
              z.xs=z.xs, pi.xs=pi.xs, params.xs=params.xs, 
              z.us=z.us, pi.us=pi.us, params.us=params.us))
}


