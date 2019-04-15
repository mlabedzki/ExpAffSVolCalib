library(polynom)
###########################################
## SABR Implied Volatility
###########################################
sigma_SABR <- function(f,K,T,beta,alpha,rho,nu,ref = FALSE){
	if(ref){
		x <- log(f/K)
		if(x == 0) I.B <- alpha*K^(beta-1)
		else if(nu == 0) I.B <- (x*alpha*(1-beta))/(f^(1-beta)-K^(1-beta))
		else if(beta == 1){
			z <- nu*x/alpha
			I.B <- (nu*x)/log((sqrt(1-2*rho*z+z^2)+z-rho)/(1-rho))
		}
		else if(beta < 1){
			z <- (nu*(f^(1-beta)-K^(1-beta)))/(alpha*(1-beta))
			I.B <- (nu*x)/log((sqrt(1-2*rho*z+z^2)+z-rho)/(1-rho))
		}
		else{}
		I.H <- 1 + ((1-beta)^2/24*alpha^2/(f*K)^(1-beta) + 1/4 * (rho * beta * nu * alpha)/(f*K)^((1-beta)/2) + (2-3*rho^2)/24*nu^2)*T
		sigma.B <- I.B * I.H
	}
	else{
		if(abs(f-K) >= 0.0001){
			z <- nu/alpha*(f*K)^((1-beta)/2)*log(f/K)
			x <- log((sqrt(1-2*rho*z+z^2)+ z - rho)/(1-rho))
			A <- alpha/((f*K)^((1-beta)/2)*1+(1-beta)^2/24*log(f/K)^2 + (1-beta)^4/1920*log(f/K))
			B <- z/x
			C <- 1 + ((1-beta)^2/24*alpha^2/(f*K)^(1-beta) + 1/4 * (rho * beta * nu * alpha)/(f*K)^((1-beta)/2) + (2-3*rho^2)/24*nu^2)*T
		}
		else{
			A <- alpha/((f)^(1-beta))
			B <- 1
			C <- 1 + ((1-beta)^2/24*alpha^2/(f)^(2-2*beta) + 1/4 * (rho * beta * nu * alpha)/(f)^(1-beta) + (2-3*rho^2)/24*nu^2)*T
		}
		sigma.B <- A*B*C
	}
	return(sigma.B)
}

###########################################
## Alpha Estimator
###########################################
getAlpha <- function(f,K,T,sigma.ATM,beta,rho,nu){
	C0 <- -sigma.ATM * f^(1-beta)
	C1 <- 1 + 1/24 * (2-3*rho^2) * nu^2 * T
	C2 <- 1/4 * rho * beta * nu * f^{beta-1} * T
	C3 <- 1/24 * (1-beta)^2 * f^(2*(beta-1)) * T
	alpha.vec <- solve(as.polynomial(c(C0,C1,C2,C3)))
	index <- which(Im(alpha.vec) == 0 & Re(alpha.vec) > 0)
	alpha <- alpha.vec[index]
	if(length(alpha) == 0) alpha <- 0
	return(min(Re(alpha)))
}

###########################################
## Volatility Fitter
###########################################
fit.vol <- function(f,K,T,beta,alpha,rho,nu,ref=FALSE){
	if(is.null(dim(K))){
		vol <- numeric(length(K))
		for (i in seq_along(K)){
			vol[i] <- sigma_SABR(f=f,K=K[i],T=T,beta=beta,alpha=alpha,rho=rho,nu=nu,ref)
		}
	}
	else{
		vol <- numeric(dim(K)[2]*length(T))
		dim(vol) <- c(length(T),dim(K)[2])
		for(i in seq_along(T)){
			for(j in 1:dim(K)[2]){
				vol[i,j] <- sigma_SABR(f[i],K[i,j],T[i],beta,alpha,rho,nu,ref)
			}
		}
	}
	return(vol)
}

###########################################
## Sum of Sqared Errors (SSE)
###########################################
SSE_SABR <- function(params,sigma.Mkt,f,K,T,beta,ref=FALSE){
	atm <- FALSE
	if(length(params) > 2){
		alpha <- params[1]
		rho <- params[2]
		nu <- params[3]
	}
	else{
		atm <- TRUE
		rho <- params[1]
		nu <- params[2]
	}
	if(is.null(dim(K))){
		error <- numeric(length(K))
		if(atm){
			sigma.ATM <- sigma.Mkt[which(K == f)]
			alpha <- getAlpha(f,f,T,sigma.ATM,beta,rho,nu)
		}
		for (i in seq_along(K)){
			error[i] <- sigma.Mkt[i] - sigma_SABR(f=f,K=K[i],T=T,beta=beta,alpha=alpha,rho=rho,nu=nu,ref)
		}
	}
	else{
		error <- numeric(dim(K)[2]*length(T))
		dim(error) <- c(length(T),dim(K)[2])
		for(i in seq_along(T)){
			if(atm){
				sigma.ATM <- sigma.Mkt[which(K == f[i])]
				alpha <- getAlpha(f[i],f[i],T[i],sigma.ATM,beta,rho,nu)
			}
			for (j in 1:dim(K)[2]){
				error[i,j] <- sigma.Mkt[i,j] - sigma_SABR(f=f[i],K=K[i,j],T=T[i],beta=beta,alpha=alpha,rho=rho,nu=nu,ref)
			}
		}
	}
	SSE <- sum(error^2,na.rm=TRUE)
	#cat(SSE,"\n")
	if(abs(rho) > 1) SSE <- Inf
	return(SSE)
}
SSE_SABRlim <- function(params,sigma.Mkt,f,K,T,beta,ref=FALSE){
	SSE_SABR(c(exp(params[1]),tanh(params[2]),exp(params[3])),sigma.Mkt,f,K,T,beta,ref=FALSE)
}
