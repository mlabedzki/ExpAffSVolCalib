source("calibration.R")

#################################
# Modified function to fitsmile #
#################################
fitsmile <- function(params,S,data,rd,rf,t,T,type="Heston"){
	theta <- params[1]
	sigma <- params[2]
	rho <- params[3]
	kappa <- params[4]
	nu0 <- params[5]
	lambda=0
	if (type=="Bates") {
		lambda <- params[6]
		muJ <- params[7]
		vJ <- params[8]
	}		
	price <- matrix(0,length(T),5)
	vol <- matrix(0,length(T),5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:length(T)){
		for(j in 1:5){
			if (type=="Heston")
				price[i,j] <- hestonc(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i] )
			else if (type=="Bates")
				price[i,j] <- callCF(cf=cfBates, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, lambda=lambda, muJ=muJ, vJ=vJ, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	vol
}
###########################################
## Implied BSM Vola
###########################################
BSMFX <- function(S,K,vol,rd,rf,t=0,T,type=1){
	d1 <- (log(S/K) + (rd - rf + vol^2/2)*(T-t))/(vol*sqrt(T-t))
	d2 <- d1 - vol*sqrt(T-t)
	type * (exp(-rf*(T-t))*S*pnorm(type*d1)-exp(-rd*(T-t))*K*pnorm(type*d2))
}
#w=vol^2*(T-t)
#f=log(S*exp((rd - rf)*(T-t)))
#exp(-rf*(T-t))*S = exp(x) * exp(-rd*t) 
B76 <- function(f,K,w,rd,t=0,T,type=1){
	d1 <- (f - log(K) + 1/2*w)/sqrt(w)
	d2 <- d1 - sqrt(w)
	type * (exp(-rd*(T-t)+f)*pnorm(type*d1)-exp(-rd*(T-t))*K*pnorm(type*d2))
}
BSMFXsolve <- function(vol, price,S,K,rd,rf,t=0,T,type=1){
	d1 <- (log(S/K) + (rd - rf + vol^2/2)*(T-t))/(vol*sqrt(T-t))
	d2 <- d1 - vol*sqrt(T-t)
	(type * (exp(-rf*(T-t))*S*pnorm(type*d1)-exp(-rd*(T-t))*K*pnorm(type*d2)) - price)^2
}
impVol <- function(data,t0,t,S,K,rd,rf,type=1){
	if(length(data)==1){ 
		nlminb(0.4, BSMFXsolve, price=data, S=S,K=K[1],rd=rd[1],rf=rf[1],t=t0,T=t[1],type=type,lower=0)$par
	}else{
		IV <- matrix(0,length(t),length(K))
		for(i in 1:length(t)){
			for(j in 1:length(K)){
				start = data/(0.4*S*exp((rd-rf)*(t-t0))*sqrt(t-t0))
				IV[i,j] <- nlminb(start, BSMFXsolve, price=data[i,j], S=S,K=K[j],rd=rd[i],rf=rf[i],t=t0,T=t[i],type=type,lower=0)$par
				#IV[i,j] <- optim(0.4, fn=BSMFXsolve, price=data[i,j], S=S,K=K[j],rd=rd[i],rf=rf[i],t=t0,T=t[i],type=type)$par
			}			
		}
		IV
	}
}
###########################################
## Implied BSM Delta Strike
###########################################
BSMDeltasolve <- function(K, vol,S,rd,rf,t=0,T,Delta,type=1){
	#Kfwd <- K*exp((rd-rf)*(T-t)) # FWD Delta
	Kfwd <- K #Spot delta must be here to have FWD outside
	d1 <- (log(S/Kfwd) + (rd - rf + vol^2/2)*(T-t))/(vol*sqrt(T-t))
	( (type*exp(-rf*(T-t))*pnorm(type*d1)) - Delta )^2
}
impStrike <- function(data, t0,t,S,rd,rf, Deltafix=""){
	if(length(data)==1){
		nlminb(S, BSMDeltasolve, vol=data, S=2,rd=rd[1],rf=rf[1],t=t0,T=t[1],Deltafix, type=1, lower=0)$par
	}else{
		if(dim(data)[2] != 5) stop("Data format is not correct!")
		Delta <- c(1-0.1,1-0.25,0.5,0.25,0.1)
		K <- matrix(0,length(t),5)
		for(i in 1:length(t)){
			for(j in 1:5){
				options(warn=-1)
				#K[i,j] <- nlminb(S, BSMDeltasolve, vol=data[i,j], S=S,rd=rd[i],rf=rf[i],t=t0,T=t[i],Delta=Delta[j],type=1, lower=0,control = list( eval.max=400,iter.max =300 ))$par
				K[i,j] <- optim(method="Brent",lower=S/3,upper=S*2, S, fn=BSMDeltasolve, vol=data[i,j], S=S,rd=rd[i],rf=rf[i], t=t0,T=t[i],Delta=Delta[j],type=1)$par
				#K[i,j] <- optim(S, fn=BSMDeltasolve, vol=data[i,j], S=S,rd=rd[i],rf=rf[i],t=t0,T=t[i],Delta=Delta[j],type=1)$par
				options(warn=0)
			}
		}
		K
	}
}
###########################################
## Characteristic Function
###########################################
# Here is version with fixed missing sigma in (b-rho*sigma*u*i+d)*(T-t)
CF <- function(u, S,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type="both"){
	if(type=="both") zeta <- c(1,-1)
	else{
		if(type==1) zeta <- 1
		if(type==2) zeta <- -1
	}
	i <- complex(real=0,imaginary=1)
	X <- log(S)
	b <- kappa + lambda - (1+zeta)/2*rho*sigma
	d <- sqrt( (rho*sigma*i*u-b)^2 - sigma^2*(zeta*i*u-u^2) ) # floor 0 is here for weird cases
	g <- (b-rho*sigma*u*i+d)/(b-rho*sigma*u*i-d)
	if(Re(g[1])==Inf){
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2* (b-rho*sigma*u*i+d)*(T-t)
	}else{
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2*( (b-rho*sigma*u*i+d)*(T-t)- 2*log((1-g*exp(d*(T-t)))/(1-g)) ) 
	}
	D <- (b-rho*sigma*u*i+d)/sigma^2 *( (1-exp(d*(T-t)))/(1-g*exp(d*(T-t))) )
	exp(C+D*nu0+i*u*X)
}
CFa <- function(u, S,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type="both"){
	#Albrecher et al. (2007)
	if(type=="both") zeta <- c(1,-1)
	else{
		if(type==1) zeta <- 1
		if(type==2) zeta <- -1
	}
	i <- complex(real=0,imaginary=1)
	X <- log(S)
	b <- kappa + lambda - (1+zeta)/2*rho*sigma
	d <- sqrt( (rho*sigma*i*u-b)^2 - sigma^2*(zeta*i*u-u^2) ) # floor 0 is here for weird cases
	c <- (b-rho*sigma*u*i-d)/(b-rho*sigma*u*i+d) #==1/g
	if(Re(c[1])==Inf){
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2* (b-rho*sigma*u*i-d)*(T-t)
	}else{
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2*( (b-rho*sigma*u*i-d)*(T-t)- 2*log((1-c*exp(-d*(T-t)))/(1-c)) ) 
	}
	D <- (b-rho*sigma*u*i-d)/sigma^2 *( (1-exp(-d*(T-t)))/(1-c*exp(-d*(T-t))) )
	exp(C+D*nu0+i*u*X)
}
CFaa <- function(u,kappa,theta,sigma,rho,lambda,nu0,t,T,type=2){
	#Albrecher et al. (2007) but without log(F) element
	if(type=="both") zeta <- c(1,-1)
	else{
		if(type==1) zeta <- 1
		if(type==2) zeta <- -1
	}
	i <- complex(real=0,imaginary=1)
	b <- kappa + lambda - (1+zeta)/2*rho*sigma
	d <- sqrt( (rho*sigma*i*u-b)^2 - sigma^2*(zeta*i*u-u^2) ) # floor 0 is here for weird cases
	c <- (b-rho*sigma*u*i-d)/(b-rho*sigma*u*i+d) #==1/g
	if(Re(c[1])==Inf){
		C <- kappa*theta/sigma^2* (b-rho*sigma*u*i-d)*(T-t)
	}else{
		C <- kappa*theta/sigma^2*( (b-rho*sigma*u*i-d)*(T-t)- 2*log((1-c*exp(-d*(T-t)))/(1-c)) ) 
	}
	D <- (b-rho*sigma*u*i-d)/sigma^2 *( (1-exp(-d*(T-t)))/(1-c*exp(-d*(T-t))) )
	exp(C+D*nu0)
}
CF_old <- function(u, S,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type="both"){
	if(type=="both") zeta <- c(1,-1)
	else{
		if(type==1) zeta <- 1
		if(type==2) zeta <- -1
	}
	i <- complex(real=0,imaginary=1)
	X <- log(S)
	b <- kappa + lambda - (1+zeta)/2*rho*sigma
	d <- sqrt( (rho*sigma*i*u-b)^2 - sigma^2*(zeta*i*u-u^2) )
	g <- (b-rho*sigma*u*i+d)/(b-rho*sigma*u*i-d)
	if(Re(g[1])==Inf){
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2* (b-rho*u*i+d)*(T-t)
	}else{
		C <- (rd-rf)*u*i*(T-t) + kappa*theta/sigma^2*( (b-rho*u*i+d)*(T-t)- 2*log((1-g*exp(d*(T-t)))/(1-g)) ) 
	}
	D <- (b-rho*sigma*u*i+d)/sigma^2 *( (1-exp(d*(T-t)))/(1-g*exp(d*(T-t))) )
	exp(C+D*nu0+i*u*X)
}
CF.Bates<-function(u,S,tau,r,q,v0,vT,rho,k,sigma,lambda,muJ,vJ) {
	i <- complex(real=0,imaginary=1)
	d = sqrt (( rho*sigma*i*u - k)^2 + sigma^2 * (i*u + u^2))
	g2 = (k - rho * sigma * i * u - d) / (k - rho * sigma * i * u + d)

	cf1 = i * u * ( log (S) + (r - q) * tau )
	cf2 = vT * k / ( sigma ^2) * ((k - rho * sigma * i * u - d) * tau - 2 * log ((1 - g2 * exp (-d * tau )) / (1 - g2)))
	cf3 = v0 / sigma ^2 * (k - rho * sigma * i * u - d) * (1 - exp (-d * tau )) / (1 - g2 * exp (-d * tau ))
	#jump
	cf4 = -lambda * muJ * i * tau * u + lambda * tau * ( (1+ muJ )^(i*u) * exp ( vJ *(i*u /2) * (i*u-1) ) -1 )
	exp( cf1 + cf2 + cf3 + cf4 )
}
CF.Double<-function(u,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type="both") {
	i <- complex(real=0,imaginary=1)
	X <- log(S)
	exp(-i*u*X)*CFa(u,S,kappa=kappa1,theta=theta1,sigma=sigma1,rho=rho1,lambda=lambda1,nu0=nu1,rd,rf,t,T,type)*CFa(u,S,kappa=kappa2,theta=theta2,sigma=sigma2,rho=rho2,lambda=lambda2,nu0=nu2,rd,rf,t,T,type)
}
###########################################
## Heston Option Pricer
###########################################
cfHeston<-function(u,S,K,T,r,v0,theta,rho,kappa,sigma){
	# In order to avoid discontinuities in the complex plane,
	# the term g2 has been modified accordingly to Kahl and Jackel (2006).

	d = sqrt((kappa - rho * sigma * 1i*u)^2 + sigma^2 * (1i*u + u ^ 2));
	g2  = (kappa - rho*sigma*1i*u - d) / (kappa - rho*sigma*1i*u + d);

	cf1 = 1i*u * (log(S) + (r) * T);
	cf2 = theta * kappa / (sigma^2) * ((kappa - rho*sigma*1i*u - d) * T - 2 * log((1 - g2 * exp(-d * T)) / (1 - g2)));
	cf3 = v0 / sigma^2 * (kappa - rho*sigma*1i*u - d) * (1 - exp(-d * T)) / (1 - g2 * exp(-d * T));
	exp(cf1 + cf2 + cf3);
}
cosmethod <- function(S,K,T,r,sigma,v0,rho,theta,kappa){
	n = 17;
	N = 2^n;
	L = 2;
	c1 = r*T + (1-exp(-kappa*T))*(theta-v0)/(2*kappa) - 0.5*theta*T;

	c2 = 1/(8*kappa^3)*(sigma*T*exp(-kappa*T)*(v0-theta)*(8*kappa*rho - 4*sigma)+kappa*rho*sigma*(1-exp(-kappa*T))*(16*theta - 8*v0) + 2*theta*kappa*T*(-4*kappa*rho*sigma + sigma^2 + 4*kappa^2) + sigma^2*((theta-2*v0)*exp(-2*kappa*T) + theta*(6*exp(-kappa*T) - 7) + 2 * v0) + 8*kappa^2*(v0 - theta)*(1 - exp(-kappa*T)));

	a = c1 - L*sqrt(abs(c2));
	b = c1 + L*sqrt(abs(c2));

	x = log(S/K);
	k = (0:(N-1));

	chipsi = coeff_h(0, a, b, 1,k);
	chi=chipsi[1]
	psi=chipsi[2]
	U = 2/(b-a)*(chi - psi);         
	unit = 0.5 * rep(1,N-1);          

	MAT = cfHeston(k*pi/(b-a),S,K,T,r,v0,theta,rho,kappa,sigma)*rep(1,length(K))*exp(1i*k*pi*(x-a)/(b-a))*(U*rep(1,length(K))); # N x K Matrix
	dim(Re(MAT))
	#ret = unit%*%(MAT); # Sum the matrix over the rows (N) per each column (K).
	#ret=K*exp(-r*T)*Re(ret);
}
#cosmethod(1.27,1.29,1,0.005,0.2,0.01,-0.3,0.01,2)
coeff_h<-function(c, a, d, cp,k) {
	# chose:
	# cp = 1 if Call
	# cp = -1 if Put

	tmp = k*pi/(b-a);
	x1 = (d-a)*tmp;
	x2 = (c-a)*tmp;
	exp_c=exp(c);
	exp_d=exp(d);
	chi = (cos(x1)*exp_d - cos(x2)*exp_c + tmp*(sin(x1)*exp_d - sin(x2)* exp_c) ) / (1+tmp^2);
	psi = (sin(x1) - sin(x2))/tmp;
	if(length(k) == 1) psi = d-c;

	psi = cp*psi;
	chi = cp*chi;
	
	c(chi,psi)
}
hestonc <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,alpha=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	Y <- log(K)
	#alpha = 5
	exp(-alpha*Y)/pi* integrate(function(u) Re(exp(-i*u*Y)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*CF(u-(alpha+1)*i,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	#lower limit of integral is problematic, originaly was 1e-16 but does not work in some cases, 1e-2 works in all cases
	#upper limit was not checked
}
hestong <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	k <- log(K)
	P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*CFa(u, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=1)),1e-16,100)$value
	P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*CFa(u, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	p1 <- (1-type)/2 + type*P1
	p2 <- (1-type)/2 + type*P2
	type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
}
hestongg <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	k <- log(K)
	stale <- CFa(-i, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)
	P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*CFa(u-i,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)/stale),1e-16,100)$value
	P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*CFa(u,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	p1 <- (1-type)/2 + type*P1
	p2 <- (1-type)/2 + type*P2
	type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
}
hestonggPre <- function(S,K,rd,rf,t,T,u,I1,I2,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	k <- log(K)
	P1 <- 1/2 + 1/pi * trapz(u,Re(exp(-i*u*k)*I1))
	P2 <- 1/2 + 1/pi * trapz(u,Re(exp(-i*u*k)*I2))
	p1 <- (1-type)/2 + type*P1
	p2 <- (1-type)/2 + type*P2
	type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
}
hestona <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	#Attari(2004) equation 12 used for P1
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	l <- log(exp(-(rd-rf)*(T-t))*K/S)
	P1 <- 1 + exp(l)*1/pi * integrate(function(u) Re(exp(-i*u*l)/(i*u-1)*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*CF(u, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*l)/(i*u)*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*CF(u, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	p1 <- (1-type)/2 + type*P1
	p2 <- (1-type)/2 + type*P2
	type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
}
hestonaa <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	#Attari(2004) last equation - one integrand
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	l <- log(exp(-(rd-rf)*(T-t))*K/S)
	integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*CF(u, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
}
hestonaa1 <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	#Attari(2004) last equation - one integrand
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	F = S*exp((rd-rf)*(T-t))
	l <- log(K/F)
	integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*log(F))*CFa(u,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)),1e-16,100)$value
	#u = (0:300)/3+(1e-16)
	#integrand <- simpson(u,Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*log(F))*CFa(u,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)))
	exp(-rd*(T-t))*(F - K * (1/2 + 1/pi * integrand))
}
hestonaa2 <- function(S,K,kappa,theta,sigma,rho,lambda,nu0,rd,rf,t,T,type=1){
	#Attari(2004) last equation - one integrand
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	F = S*exp((rd-rf)*(T-t))
	l <- log(K/F)
	#u = (0:300)/3+(1e-16)
	u = c(1e-16,(1:20)*0.1,2+(1:10)*0.2,4+(1:30)/3,14+(1:86))
	phi = exp(-i*u*log(F))*CFa(u,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=2)
	integrand <- trapz(u,((Re(phi)+Im(phi)/u)*cos(u*l)+(Im(phi)-Re(phi)/u)*sin(u*l))/(1+u*u))
	exp(-rd*(T-t))*(F - K * (1/2 + 1/pi * integrand))
}
#hestonc(S=1.27,K=1.31,kappa=2,theta=0.1^2,sigma=0.3,rho=-0.5,lambda=0,nu0=0.1^2,rd=0.05,rf=0.02,t=0,T=2,alpha=1)
#heston1(S=1.27,K=1.31,kappa=2,theta=0.1^2,sigma=0.3,rho=-0.5,lambda=0,nu0=0.1^2,rd=0.05,rf=0.02,t=0,T=2)
#hestona(S=1.27,K=1.31,kappa=2,theta=0.1^2,sigma=0.3,rho=-0.5,lambda=0,nu0=0.1^2,rd=0.05,rf=0.02,t=0,T=2,type=1) # 0.5132753 0.4601121
#hestong(S=1.27,K=1.31,kappa=2,theta=0.1^2,sigma=0.3,rho=-0.5,lambda=0,nu0=0.1^2,rd=0.05,rf=0.02,t=0,T=2,type=1) #0.5195925,  0.4690809
heston2 <- function(S,K,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,rd,rf,t,T,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	Y <- log(K)
	P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*CF.Double(u,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=1)),1e-16,100)$value
	P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*CF.Double(u,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=2)),1e-16,100)$value
	stale <- CF.Double(-i,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=2)
	P12 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*CF.Double(u-i,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=2)/stale),1e-16,100)$value
	#p1 <- (1-type)/2 + type*P1
	#p2 <- (1-type)/2 + type*P2
	#type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
	exp(-rf*(T-t))*S*P12 - exp(-rd*(T-t))*K*P2
	#c(P12,P2)
	#k <- log(K)-(rd-rf)*(T-t)
	#integrand <- integrate(function(u) 1/(u^2+1/4)*Re(exp(-i*u*k)*CF.Double(u-i/2,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=2)),1e-16,100)$value
	#S - exp(-rd*(T-t)/2) * sqrt(S*K) * 1/pi * integrand
}
#heston2(S=1.27,K=1.31,kappa1=110,theta1=0.04,sigma1=0.05,rho1=0.9,lambda1=0,nu1=0.04,kappa2=7.5,theta2=0.03,sigma2=0.5,rho2=-0.95,lambda2=0,nu2=0.0003,rd=0.05,rf=0.02,t=0,T=2,type=1)
#exp(-1i*(-1i)*log(1.27))*CFa(-1i,S=1.27,kappa=7.5,theta=0.03,sigma=0.5,rho=-0.95,lambda=0,nu0=0.0003,rd=0.05,rf=0.02,t=0,T=2,type=2)*CFa(-1i,S=1.27,kappa=110,theta=0.04,sigma=0.05,rho=0.9,lambda=0,nu0=0.04,rd=0.05,rf=0.02,t=0,T=2,type=2)
sz2d <- function(S,K,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2,rd,rf,t,T,type=1,alpha=3){
	i <- complex(real=0,imaginary=1)
	Y <- log(K)
	#P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*cfSZd(u,S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2)),1e-16,100)$value
	#P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*cfSZd(u,S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2)),1e-16,100)$value
	#stale <- cfSZd(-i,S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2)
	#P12 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*Y)/(i*u)*cfSZd(u-i,S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2)/stale),1e-16,100)$value
	#p1 <- (1-type)/2 + type*P1
	#p2 <- (1-type)/2 + type*P2
	#type * ( exp(-rf*(T-t))*S*p1 - exp(-rd*(T-t))*K*p2 )
	#exp(-rf*(T-t))*S*P12 - exp(-rd*(T-t))*K*P2
	#c(P12,P2)
	#k <- log(K)-(rd-rf)*(T-t)
	#integrand <- integrate(function(u) 1/(u^2+1/4)*Re(exp(-i*u*k)*CF.Double(u-i/2,S,rd,rf,t,T,kappa1,theta1,sigma1,rho1,lambda1,nu1,kappa2,theta2,sigma2,rho2,lambda2,nu2,type=2)),1e-16,100)$value
	#S - exp(-rd*(T-t)/2) * sqrt(S*K) * 1/pi * integrand
	exp(-alpha*Y)/pi* integrate(function(u) Re(exp(-i*u*Y)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfSZd(u-i*(alpha+1),S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2)),1e-16,100)$value
}
#sz2d(S=1.27,K=1.31,kappa1=10,theta1=0.4,sigma1=0.05,rho1=0.9,nu1=0.4,kappa2=5,theta2=0.3,sigma2=0.5,rho2=-0.95,nu2=0.03,rd=0.05,rf=0.02,t=0,T=2)
#cfSZd(1,S=1.27,kappa1=10,theta1=0.4,sigma1=0.05,rho1=0.9,nu1=0.4,kappa2=5,theta2=0.3,sigma2=0.5,rho2=-0.95,nu2=0.03,rd=0.05,rf=0.02,t=0,T=2)
#exp(-1i*1*log(1.27))*cfSZlk(1,S=1.27,kappa=10,theta=0.4,sigma=0.05,rho=0.9,nu0=0.4,rd=0.05,rf=0.02,t=0,T=2)*cfSZlk(1,S=1.27,kappa=10,theta=0.3,sigma=0.5,rho=-0.9,nu0=0.03,rd=0.05,rf=0.02,t=0,T=2)


szhw <- function(S,K,kappa,ksi,tau,nu0,theta,a,sigma,rhorv,rhoxr,rhoxv,rd,rf,t,T,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	Y <- log(K)
	alpha = 5
	exp(-alpha*Y)/pi* integrate(function(u) Re(exp(-i*u*Y)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfSZHW(u-(alpha+1)*i,S=S,kappa=kappa,ksi=ksi,tau=tau,nu0=nu0,theta=theta,a=a,sigma=sigma,rhorv=rhorv,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,10)$value
}
szhwa <- function(S,K,kappa,ksi,tau,nu0,theta,a,sigma,rhorv,rhoxr,rhoxv,rd,rf,t,T,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	l <- log(exp(-(rd-rf)*(T-t))*K/S)
	integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*cfSZHW(u,S=S,kappa=kappa,ksi=ksi,tau=tau,nu0=nu0,theta=theta,a=a,sigma=sigma,rhorv=rhorv,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,20)$value
	exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
}
#szhw(S=1.27,K=1.29,kappa=2,ksi=0.15,tau=0.3,nu0=0.1,theta=0.04,a=2,sigma=0.8,0,0,-0.6,rd=0.005,rf=0,t=0,T=2)
#szhwa(1.27,1.31,2,0.2,0.4,0.3,0.01,0.5,0.8,0,0,-0.6,0.005,0,0,2)
bond <- function(t,T) {
	(1-exp(-a*(T-t)))/a
}
szhwC <- function(u,t,T,g0,g,a) {
	-u*(1i+u)*((g[3]-g[4]*exp(-2*g0*T))-(g[5]*exp(-a*T)-g[6]*exp(-(2*g0+a)*T))-g[7]*exp(-g0*T))/(g[1]+g[2]*exp(-2*g0*T))
}
szhwD <- function(u,t,T,g0,g) {
	-u*(1i+u)*(1-exp(-2*g0*T))/(g[1]+g[2]*exp(-2*g0*T))
}
cfSZHW<- function(u,S,kappa,ksi,tau,nu0,theta,a,sigma,rhorv,rhoxr,rhoxv,rd,rf,t,T){
	y <- log(S)
	V = sigma^2/a^2 * (T + 2/a * exp(-a*T) - 1/2/a * exp(-2*a*T) - 3/2/a)
	g0 <- sqrt((kappa-rhoxv*tau*1i*u)^2+tau^2*u*(1i+u))
	g1 <- g0 + (kappa - rhoxv*tau*1i*u)
	g2 <- g0 - (kappa - rhoxv*tau*1i*u)
	g3 <- (rhoxr*sigma*g1+kappa*a*ksi+rhorv*sigma*tau*(1+1i*u))/a/g0
	g4 <- (rhoxr*sigma*g2-kappa*a*ksi-rhorv*sigma*tau*(1+1i*u))/a/g0
	g5 <- (rhoxr*sigma*g1+rhorv*sigma*tau*(1+1i*u))/a/(g0-a)
	g6 <- (rhoxr*sigma*g2-rhorv*sigma*tau*(1+1i*u))/a/(g0+a)
	g7 <- (g3-g4)-(g5-g6)
	g <- c(g1,g2,g3,g4,g5,g6,g7)
	B <- 1i*u	
	A <- -0.5*u*(1i+u)*V + integrate(function(s) Re((kappa*ksi+rhorv*(1+1i*u)*tau*sigma*bond(0,s))*szhwC(u,0,s,g0,g,a)+0.5*tau^2*(szhwC(u,0,s,g0,g,a)^2+szhwD(u,0,s,g0,g))),1e-16,T)$value + 1i*integrate(function(s) Im((kappa*ksi+rhorv*(1+1i*u)*tau*sigma*bond(0,s))*szhwC(u,0,s,g0,g,a)+0.5*tau^2*(szhwC(u,0,s,g0,g,a)^2+szhwD(u,0,s,g0,g))),1e-16,T)$value
	exp(A+B*y+szhwC(u,0,T,g0,g,a)*nu0+0.5*szhwD(u,0,T,g0,g)*nu0^2)
}
#cfSZHW(u=1,S=1.27,kappa=2,ksi=0.15,tau=0.3,nu0=0.1,theta=0.04,a=2,sigma=0.8,0,0,-0.6,rd=0.005,rf=0,t=0,T=2)
fv<-function(tau,beta,d,lambda,gamma,kappa,sigma) {
	o = exp(tau/2)
	g = (beta-d)/(beta+d)
	f0 = o^2/(o^(2*d)-g)
	f1 = 16*kappa*sigma/4/gamma^2/d * (beta-d) * sinh(tau*d/4)^2
	f2 = 2/d * ((o^d-1) + g*(o^-d-1))
	f3 = 2 * (o^(d-2*lambda)-1)/(d-2*lambda) - 2*g*(1-o^(2*lambda-d))/(d+2*lambda)
	f4 = 2/(d-2*lambda) - 4/d + 2/(d+2*lambda)
	f5 = o^(-2*lambda-d)*(2/d*o^(2*lambda)*(1+o^(2*d))-2*o^(2*d)/(d-2*lambda)-2/(d+2*lambda))
	c(f1,f2,f3,f4,f5,f0)
}
fBs<-function(u,tau,beta,d,lambda,gamma,kappa,sigma,eta,rhoxr,rhorv) {
	f = fv(tau,beta,d,lambda,gamma,kappa,sigma)
	f[6] * (f[1] + 1/lambda*(1i*u-1)*(eta*rhoxr*1i*u*(f[2]-f[3])+eta*rhorv/2/gamma*(beta-d)*(f[4]+f[5])))
}
fBr<-function(u,tau,lambda) {
	1/lambda * (1i*u-1) * (1-exp(-tau*lambda))
}
frt<-function(t,r0,lambda,theta) { #==ksi
	#in case of time dependent theta:
	#r0*exp(-lambda*t) + integrate(function(s) theta*exp(-lambda*(t-s)),0,t)$value
	r0*exp(-lambda*t) + theta*(1-exp(-lambda*t))
}
frtSEE<-function(pars,data,tau) {
	lambda=exp(pars[2])
	r0=pars[1]
	theta=pars[3]
	n=length(tau)
	err=rep(0,n)
	for (i in 1:n) {
		err[i]=frt(tau[i],r0,lambda,theta)-data[i]
	}
	sqrt( sum( (err)^2 ) )
}
rfit<-function(rdiff,tau){
	n=dim(rdiff)[1]
	r0=rep(0,n);lambda=r0;theta=r0
	for (i in 1:n) {
		opt=optim(c(rdiff[i,1],1,rdiff[i,6]),frtSEE,data=rdiff[i,],tau=tau,method="BFGS")$par
		r0[i]=opt[1];theta[i]=opt[3];lambda[i]=exp(opt[2])
	}
	cbind(r0,lambda,theta)
}
cfGrzelak<- function(u,S,lambda,theta,r0,eta,kappa,sigma,s0,gamma,rhorv,rhoxr,rhoxv,rd,rf,t,T){
	tau = T-t
	o = exp(tau/2)
	x <- log(S)
	r20 = 0
	beta = 2*(kappa-rhoxv*gamma*u*1i)
	alpha = -1/2*u*(1i+u)
	d = sqrt(beta^2-8*alpha*gamma^2)
	g = (beta-d)/(beta+d)
	f6 = 1/4/gamma^2 * (beta-d) * tau
	f7 = (1i*u-1)^2*(3+o^(-4*lambda)-4*o^(-2*lambda)-2*tau*lambda)
	Gint = integrate(function(s) Re(fBs(u,s,beta,d,lambda,gamma,kappa,sigma,eta,rhoxr,rhorv)*(kappa*sigma+1/2*gamma^2*fBs(u,s,beta,d,lambda,gamma,kappa,sigma,eta,rhoxr,rhorv)+eta*rhorv*gamma*fBr(u,s,lambda))),0,tau)$value
	
	Br <- fBr(u,tau,lambda)
	Bs <- fBs(u,tau,beta,d,lambda,gamma,kappa,sigma,eta,rhoxr,rhorv)
	Bx <- 1i*u
	Bv <- 1/4/gamma^2 * ((1-o^(-2*d))/(1-g*o^(-2*d))) * (beta-d)
	A <- f6 - 1/2/gamma^2 * log((g*o^(-2*d)-1)/(g-1)) - 1/2/lambda^3 * f7 + Gint
	#Ab <- eta^2/4/lambda^3 * (1+2*lambda*tau-(exp(-lambda*tau)-2)^2)
	#P(0,tau) = exp(-integrate()$value+Ab)
	Upsilon  = (1i*u-1) * integrate(function(s) frt(s,r0,lambda,theta),t,T)$value
	#Upsilon = -(1i*u-1) * (log(P(0,T)/P(0,t))) + eta^2/2/lambda^2*(tau+2/lambda*(exp(-lambda*tau)-1)-1/2/lambda*(exp(-2*lambda*tau)-1)))
	A2 <- A + Upsilon
	exp(A2 + Bx*x + Br*r20 + Bv*s0^2 + Bs*s0)
}
#cfGrzelak(u=1,S=1.27,lambda=2,theta=0.02,r0=0.01,eta=0.8,kappa=2,sigma=0.1,s0=0.15,gamma=0.3,0,0,-0.6,rd=0.005,rf=0,t=0,T=2)
szhwg <- function(S,K,lambda,theta,r0,eta,kappa,sigma,s0,gamma,rhorv,rhoxr,rhoxv,rd,rf,t,T){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	#l <- log(exp(-(rd-rf)*(T-t))*K/S)
	#integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*cfGrzelak(u,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,s0=s0,gamma=gamma,rhorv=rhorv,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,20)$value
	#exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
	Y <- log(K)
	alpha = 5
	exp(-alpha*Y)/pi* integrate(function(u) Re(exp(-i*u*Y)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfGrzelak(u,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,s0=s0,gamma=gamma,rhorv=rhorv,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,10)$value
}
#szhwg(S=1.27,K=1.29,lambda=2,theta=0.02,r0=0.01,eta=0.8,kappa=2,sigma=0.1,s0=0.15,gamma=0.3,0,0,-0.6,rd=0.005,rf=0,t=0,T=1)
#szhwg(S=1.2701,K=1.29,lambda=2,theta=0.02,r0=0.01,eta=0.8,kappa=2,sigma=0.1,s0=0.15,gamma=0.3,0,0,-0.6,rd=0.005,rf=0,t=0,T=1)

cfHHW <- function(u,S,lambda,theta,r0,eta,kappa,sigma,nu0,gamma,rhoxr,rhoxv,rd,rf,t,T){
	tau = T-t
	o = exp(tau/2)
	x <- log(S)
	r20 = 0
	beta = kappa-gamma*rhoxv*u*1i
	d = sqrt(beta^2+(1i*u+u^2)*gamma^2)
	g = (beta-d)/(beta+d)	
	Bx <- 1i*u
	Br <- fBr(u,tau,lambda)
	Bv <- 1/gamma^2 * ((1-o^(-2*d))/(1-g*o^(-2*d))) * (beta-d)
	
	f1 = -1/2*rhoxr^2*u*(1i+u)
	f2 = kappa*sigma
	f3 = rhoxr*eta*u*1i
	f4 = 1/2*eta^2
	f5 = 1/gamma^2 * (beta-d)
	f6 = 1/lambda*(1i*u-1)
	
	G1 = -f4*f6^2 + 2*f6*(f3+2*f4*f6)*o^(2*lambda)
	G2 = -f6*(2*f3+3*f4*f6)+2*(f1+f2*f5+f6*(f3+f4*f6))*lambda*tau
	G3 = 2*f2*f5*o^(4*lambda)*(g-1)*lambda
	
	if (g[1]==0)
		A <- 0
	else
		A <- o^(-4*lambda)/2/d/g/lambda*(d*g*(G1+o^(4*lambda)*G2)+G3*log((g*o^(-2*d)-1)/(g-1)))
	#Ab <- eta^2/4/lambda^3 * (1+2*lambda*tau-(exp(-lambda*tau)-2)^2)
	#P(0,tau) = exp(-integrate()$value+Ab)
	Upsilon  = (1i*u-1) * integrate(function(s) frt(s,r0,lambda,theta),t,T)$value
	#Upsilon = -(1i*u-1) * (log(P(0,T)/P(0,t))) + eta^2/2/lambda^2*(tau+2/lambda*(exp(-lambda*tau)-1)-1/2/lambda*(exp(-2*lambda*tau)-1)))
	#Upsilon = -(1i*u-1) * (r0*tau + eta^2/2/lambda^2*(tau+2/lambda*(exp(-lambda*tau)-1)-1/2/lambda*(exp(-2*lambda*tau)-1)))
	A2 <- A + Upsilon
	exp(A2 + Bx*x + Br*r0 + Bv*nu0)
}
cfHHW(u=-1i,S=1.27,lambda=2,theta=0.02,r0=0.01,eta=0.8,kappa=2,sigma=0.1^2,nu0=0.15^2,gamma=0.3,0,-0.6,rd=0.005,rf=0,t=0,T=0.1)

hhw <- function(S,K,lambda,theta,r0,eta,kappa,sigma,nu0,gamma,rhoxr,rhoxv,rd,rf,t,T){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	k <- log(K)
	#alpha = 0.001
	l <- log(exp(-(rd-rf)*(T-t))*K/S)
	integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*cfHHW(u,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,nu0=nu0,gamma=gamma,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,100)$value
	#P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfHHW(u,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,nu0=nu0,gamma=gamma,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,100)$value
	#stale=cfHHW(-i,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,nu0=nu0,gamma=gamma,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)
	#P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfHHW(u-i,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,nu0=nu0,gamma=gamma,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)/stale),1e-16,100)$value
	exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
	#exp(-rf*(T-t))*S*P1 - K*exp(-rd*(T-t)) * P2
	#exp(-alpha*k)/pi* integrate(function(u) Re(exp(-i*u*k)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfHHW(u,S=S,lambda=lambda,theta=theta,r0=r0,eta=eta,kappa=kappa,sigma=sigma,nu0=nu0,gamma=gamma,rhoxr=rhoxr,rhoxv=rhoxv,t=t,T=T)),1e-16,100)$value
}
#hhw(S=1.27,K=1.315,lambda=4,theta=0.005,r0=0.005,eta=0.005,kappa=1.33,sigma=0.03,nu0=0.054,gamma=0.1,0,-0.6,rd=0.005,rf=0,t=0,T=1)
#(heston(S=1.27,K=1.315,kappa=1.33,theta=0.03,sigma=0.1,rho=-0.6,lambda=0,nu0=0.054,rd=0.005,rf=0,t=0,T=1)-heston(S=1.2701,K=1.315,kappa=1.33,theta=0.03,sigma=0.1,rho=-0.6,lambda=0,nu0=0.054,rd=0.005,rf=0,t=0,T=1))/0.0001
cfSZ <- function(u,S,theta,sigma,rho,kappa,nu0,rd,rf,t,T,type=2){
	if(type==1) z <- 1
	if(type==2) z <- 0
	tau = T-t
	x <- log(S)
	Bx = 1i*u
	#s1=-1/2*1i*u*(1i*u*(1-rho^2)-1+2*rho*kappa/sigma)
	s1=-1/2*(z+1i*u)^2*(1-rho^2)+(z+1i*u)*(1-2*rho*kappa/sigma)
	s2=(z+1i*u)*rho*kappa*theta/sigma
	s3=(z+1i*u)*rho/2/sigma
	g1=sqrt(2*sigma^2*s1+kappa^2)
	#g2=kappa^2*theta-s2*sigma^2
	#g3=1/g1*(kappa-2*sigma^2*s3)
	g2=1/g1*(kappa-2*sigma^2*s3)
	g3=kappa^2*theta-s2*sigma^2
	g4=cosh(g1*tau)+g2*sinh(g1*tau)
	Bs=1/sigma^2/g1*((kappa*theta*g1-g2*g3+g3*(sinh(g1*tau)+g2*cosh(g1*tau)))/(cosh(g1*tau)+g2*sinh(g1*tau))-kappa*theta*g1)
	Bv=1/sigma^2*(kappa-g1*(sinh(g1*tau)+g2*cosh(g1*tau))/(cosh(g1*tau)+g2*sinh(g1*tau)))
	#Bs=1/g1/g4/sigma^2*((kappa*theta*g1-g2*g3)*(1-cosh(g1*tau))-(kappa*theta*g1*g2-g3)*sinh(g1*tau))
	#Bv=kappa/sigma^2-g1/sigma^2/g4*(sinh(g1*tau)+g2*cosh(g1*tau))
	if (g4[1]==0)
		A = 0
	else
		A = -1/2*log(cosh(g1*tau)+g2*sinh(g1*tau))+1/2*kappa*tau+(kappa^2*theta^2*g1^2-g3^2)/2/sigma^2/g1^3*(sinh(g1*tau)/(cosh(g1*tau)+g2*sinh(g1*tau))-g1*tau)+(kappa*theta*g1-g2*g3)*g3/sigma^2/g1^3*(cosh(g1*tau)-1)/(cosh(g1*tau)+g2*sinh(g1*tau))
		#A=-1/2*log(g4)+sinh(g1*tau)/2/g1^3/g4/sigma^2*((kappa*theta*g1-g2*g3)^2-g3^2*(1-g2^2))+g3/g1^3/g4/sigma^2*(kappa*theta*g1-g2*g3)*(g4-1)+tau/2/g1^2/sigma^2*(kappa*g1^2*(sigma^2-kappa*theta^2)+g3^2)
	exp(Bx*(x+(rd-rf)*tau))*exp(-s3*nu0^2-1/2*(z+1i*u)*rho*sigma*tau)*exp(A + Bs*nu0 + 1/2*Bv*nu0^2)
}
cfSZlk <- function(u,S,theta,sigma,rho,kappa,nu0,rd,rf,t,T){
	tau = T-t
	x <- log(S*exp((rd-rf)*tau))
	alpha = -1/2*u*(u+1i)
	beta = 2*(kappa-1i*sigma*rho*u)
	gamma=2*sigma^2
	d=sqrt(beta^2-4*alpha*gamma)
	g=(beta-d)/(beta+d)
	Bx = 1i*u
	Bs = 2*kappa*theta*(beta-d)*(1-exp(-d/2*tau))^2/d/gamma/(1-g*exp(-d*tau))
	Bv = (beta-d)*(1-exp(-d*tau))/2/gamma/(1-g*exp(-d*tau))
	if (g[1]==0)
		A = 0
	else
		A = (beta-d)*kappa^2*theta^2/2/d^3/sigma^2*(beta*(d*tau-4)+d*(d*tau-2)+(4*exp(-d/2*tau)*((d^2-2*beta^2)/(beta+d)*exp(-d/2*tau)+2*beta))/(1-g*exp(-d*tau))) + 1/4*(beta-d)*tau-1/2*log((g*exp(-d*tau)-1)/(g-1))
		#A = (beta-d)*kappa^2*theta^2/2/d^3/sigma^2*((d^3 *tau *(e^(d *tau)+1)+2 *d^2 *(beta*tau*e^(d*tau)+2)+8*beta^2*(e^((d*tau)/2)-1)+d*beta*(beta*tau*e^(d*tau)+8*e^((d*tau)/2)-beta*tau))/((beta*(e^(d*tau)-1)+d*(e^(d*tau)+1)))) + 1/4*(beta-d)*tau-1/2*log((g*exp(-d*tau)-1)/(g-1))	#alternative formula
	exp(A + Bx*x + Bs*nu0 + Bv*nu0^2)
}
cfSZd<-function(u,S,rd,rf,t,T,theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2) {
	x <- log(S*exp((rd-rf)*(T-t)))
	exp(-1i*u*x)*cfSZlk(u,S,kappa=kappa1,theta=theta1,sigma=sigma1,rho=rho1,nu0=nu1,rd,rf,t,T)*cfSZlk(u,S,kappa=kappa2,theta=theta2,sigma=sigma2,rho=rho2,nu0=nu2,rd,rf,t,T)
}
#hestonaa(S=1.27,K=1.28,kappa=1,theta=VIX[i,1]^2,sigma=1.2*VIX[i,1],rho=-0.6,lambda=0,nu0=VIX[i,1]^2,rd=rd[4],rf=rf[4],t=0,T=tau[4])
sz <- function(S,K,theta,sigma,rho,kappa,nu0,rd,rf,t,T,alpha=5){
	i <- complex(real=0,imaginary=1); k <- log(K); #alpha = 1.5
	#l <- log(exp(-(rd-rf)*(T-t))*K/S)
	#integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*cfSZlk(u,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
	#P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfSZlk(u,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
	#stale = cfSZlk(-i,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)
	#P12 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfSZlk(u-i,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)/stale),1e-16,100)$value
	#P1 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfSZ(u,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T,type=1)),1e-16,100)$value
	#exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
	#exp(-rf*(T-t))*S*P12 - K*exp(-rd*(T-t)) * P2
	exp(-alpha*k)/pi* integrate(function(u) Re(exp(-i*u*k)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfSZlk(u-i*(alpha+1),S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
}
#0.04548
#sz(S=1.27,K=1.28,kappa=1,theta=VIX[i,1],sigma=1.2*VIX[i,1],rho=-0.6,nu0=VIX[i,1],rd=rd[4],rf=rf[4],t=0,T=tau[4],alpha=10)


cfHSC <- function(u,S,theta,sigma,rho0,kappa,nu0,rho2,thetar,kappar,sigmar,rd,rf,t,T){
	tau = T-t
	x <- log(S)
	r = rd-rf
	m = sqrt(theta-sigma^2/8/kappa)
	n = sqrt(nu0)-m
	k = sqrt((nu0*exp(-kappa)-sigma^2*(1-exp(-kappa))/4/kappa)+theta*(1-exp(-kappa))+sigma^2*theta*(1-exp(-kappa))^2/(8*kappa*theta+8*kappa*exp(-kappa)*(nu0-theta)))
	l = -log(1/n*(k-m))
	
	d = sqrt(kappa^2+(1i*u+u^2)*sigma^2)
	g = (kappa-d)/(kappa+d)	
	alpha = kappar*thetar+m*sigmar*rho2*u*1i
	beta = n*sigmar*rho2*u*1i
	l1 = -log((exp(-d)-g*exp(-d))/(1-g*exp(-d)))
	c1 = 1i*u*(kappa-d)/sigma^2
	c2 = (theta-nu0)/(kappa+kappar-l1)*exp(-kappa*T)+(nu0-theta)/(kappa+kappar)*exp(-kappa*T)-theta/kappar+1/(kappar-l1)
	
	h1 = (1i*u-1)*r*tau+kappa*theta/sigma^2*((kappa-d)*tau-2*log((1-g*exp(-d*tau))/(1-g)))
	h2c = -c1*(theta-nu0)*exp(T*(-kappa))/(kappa+kappar-l1)/(kappa-l1)-c1*(nu0-theta)*exp((-T)*(kappa))/(kappa)/(kappa+kappar)-c1*theta/(kappar-l1)/(l1)+c1*c2/(kappar)
	h3c = c1*(theta-nu0)*exp(T*(kappa+l))/(kappa+kappar-l1)/(kappa+l-l1)+c1*(nu0-theta)*exp((-T)*(l+kappa))/(l+kappa)/(kappa+kappar)+c1*theta*exp(-l*T)/kappar/l-c1*theta*exp(-l*T)/(kappar-l1)/(l-l1)+c1*c2*exp(-l*T)/(l-kappar)
	h2 = c1*(theta-nu0)*exp(kappa*(tau-T)-l1*tau)/(kappa+kappar-l1)/(kappa-l1)+c1*(nu0-theta)*exp(kappa*(tau-T))/kappa/(kappa+kappar)+c1*theta*tau/kappar+c1*theta*exp(-l1*tau)/(kappar-l1)/l1-c1*c2*exp(-kappar*tau)/kappar+h2c
	h3 = c1*(theta-nu0)*exp(tau*(kappa+l-l1)-T*(kappa+l))/(kappa+kappar-l1)/(kappa+l-l1)+c1*(nu0-theta)*exp((tau-T)*(l+kappa))/(l+kappa)/(kappa+kappar)+c1*theta*exp(l*(tau-T))/kappar/l-c1*theta*exp(tau*(l-l1)-l*T)/(kappar-l1)/(l-l1)+c1*c2*exp(tau*(l-kappar)-l*T)/(l-kappar)+h3c
	h4c1 = c1^2*(nu0-theta)^2/2/(kappa)/(kappa+kappar)^2
	h4c2 = c1^2*(nu0-theta)^2/2/(kappa-l)/(2*kappa+kappar-l1)^2
	h4c3 = 2*c1^2*(nu0-theta)^2/(l1-2*kappa)/(kappa+kappar-l)/(kappa+kappar)
	h4c4 = 2*c1^2*(theta)^2/(kappar-l1)/(kappa)/l1
	h4c5 = 2*c1^2*c2*(theta)/(kappar^2-l1^2)
	h4c6 = -c1^2*c2^2/2/kappar
	h4c7 = -c1^2*c2*(theta)*2/(kappar)^2
	h4c8 = -c1^2*(theta)^2/2/(l1)/(kappar-l1)^2
	h4c9 = c1^2*c2*(nu0-theta)*2/(kappa-kappar)/(kappa+kappar)
	h4c10 = -2*c1^2*(nu0*theta-theta^2)/(kappa-l1)/(kappa+kappar)/(kappar-l1)
	h4c11 = -2*c1^2*c2*(nu0-theta)/(kappa+kappar-l1)/(kappa-kappar-l1)
	h4c12 = -2*c1^2*(nu0*theta-theta^2)/(kappa-l1)/kappar/(kappa+kappar-l1)
	h4c13 = 2*c1^2*(nu0*theta-theta^2)/(kappa-2*l)/(kappar-l1)/(kappa+kappar-l1)
	h4c14 = 2*c1^2*(nu0*theta-theta^2)/(kappa)/kappar/(kappa+kappar)^2
	h4c = (h4c1+h4c2+h4c3)*exp(-2*kappa*T)+h4c4+h4c5+h4c6+h4c7+h4c8+(h4c9+h4c10+h4c11+h4c12+h4c13+h4c14)*exp(-kappa*T)
	h4 = h4c1*exp(2*kappa*(tau-T)) + h4c4*exp(-tau*l1) + h4c5*exp(-l1-kappar*(tau)) + h4c9*exp((kappa-kappar)*tau-kappa*T) + h4c2*exp(2*(kappa-l1)*tau-2*(kappa)*T) + h4c3*exp((2*kappa-l1)*tau-2*(kappa)*T) + h4c11*exp((kappa-kappar-l1)*tau-(kappa)*T) + h4c12*exp((kappa-l1)*tau-(kappa)*T) + h4c13*exp((kappa-2*l1)*tau-(kappa)*T) + h4c6*exp(1*-2*(kappar)*tau) + h4c7*exp(-kappar*tau) + h4c8*exp(-2*l1*tau) + h4c10*exp((kappa-l1)*tau-(kappa)*T) + h4c14*exp((kappar)*(tau-T)) + c1^2*theta^3*tau/kappar^2 + h4c
	
	Bx <- 1i*u
	Bv <- (kappa-d)/sigma^2*(1-exp(-d*tau))/(1-g*exp(-d*tau))
	Br <- c1*(theta-nu0)/(kappa+kappar-l1)*exp((kappa-l1)*tau-kappa*T)+c1*(nu0-theta)/(kappa+kappar)*exp((kappa)*(tau-T))+c1*theta/kappar-c1*theta/(kappar-l1)*exp(-l1)+c1*c2*exp(-kappar*tau)
	A <- h1+alpha*h2+beta*h3+sigmar^2/2*h4

	exp(-(rd-rf)*tau + A + Bx*x + Bv*nu0 + Br*rho0)
}
#cfHSC(u=1,S=1.27,kappa=1,theta=VIX[i,1]^2,sigma=1.2*VIX[i,1],rho0=-0.6,nu0=VIX[i,1]^2,rho2=0,thetar=-0.6,kappar=3,sigmar=0.01,rd=rd[4],rf=rf[4],t=0,T=tau[4])
HSC <- function(S,K,theta,sigma,rho0,kappa,nu0,rho2,thetar,kappar,sigmar,rd,rf,t,T,alpha=10){
	i <- complex(real=0,imaginary=1); k <- log(K); #alpha = 1.5
	#l <- log(exp(-(rd-rf)*(T-t))*K/S)
	#integrand <- integrate(function(u) Re(exp(-i*u*l)*((1-i/u)/(1+u^2))*exp(-i*u*(log(S)+(rd-rf)*(T-t)))*cfSZlk(u,S=S,theta=theta,sigma=sigma,rho=rho,kappa=kappa,nu0=nu0,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
	P2 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfHSC(u,S=S,theta=theta,sigma=sigma,rho0=rho0,kappa=kappa,nu0=nu0,rho2=rho2,thetar=thetar,kappar=kappar,sigmar=sigmar,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
	stale = cfHSC(-i,S=S,theta=theta,sigma=sigma,rho0=rho0,kappa=kappa,nu0=nu0,rho2=rho2,thetar=thetar,kappar=kappar,sigmar=sigmar,rd=rd,rf=rf,t=t,T=T)
	#P12 <- 1/2 + 1/pi * integrate(function(u) Re(exp(-i*u*k)/(i*u)*cfHSC(u-i,S=S,theta=theta,sigma=sigma,rho0=rho0,kappa=kappa,nu0=nu0,rho2=rho2,thetar=thetar,kappar=kappar,sigmar=sigmar,rd=rd,rf=rf,t=t,T=T)/stale),1e-16,100)$value
	#exp(-rf*(T-t))*S - K*exp(-rd*(T-t)) * (1/2 + 1/pi * integrand)
	#exp(-rf*(T-t))*S*P12 - K*exp(-rd*(T-t)) * P2
	#P12
	exp(-alpha*k)/pi* integrate(function(u) Re(exp(-i*u*k)*exp(-rd*(T-t))/(alpha^2 + alpha - u^2 + i*(2*alpha + 1)*u)*cfHSC(u-i*(alpha+1),S=S,theta=theta,sigma=sigma,rho0=rho0,kappa=kappa,nu0=nu0,rho2=rho2,thetar=thetar,kappar=kappar,sigmar=sigmar,rd=rd,rf=rf,t=t,T=T)),1e-16,100)$value
}
#HSC(S=1.27,K=1.28,kappa=1,theta=VIX[i,1]^2,sigma=1.2*VIX[i,1],rho0=-0.6,nu0=VIX[i,1]^2,rho2=0,thetar=-0.6,kappar=3,sigmar=0.01,rd=rd[4],rf=rf[4],t=0,T=tau[4],alpha=0.5)

###########################################
## Calibration
###########################################
## Sum of Squared Errors Vol
SSE_old <- function(params,data,S,kappa,lambda,nu0,rd,rf,t,T,w){
	theta <- params[1]
	sigma <- params[2]
	rho <- params[3]
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i] )
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	sqrt( sum( (vol - data)^2 ) )
}
## Fitting of the volatility surface
fitsmile_old <- function(params, S,data, kappa,lambda,nu0,rd,rf,t,T){
	theta <- params[1]
	sigma <- params[2]
	rho <- params[3]
	price <- matrix(0,length(T),5)
	vol <- matrix(0,length(T),5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:length(T)){
		for(j in 1:5){
			price[i,j] <- heston(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i] )
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	vol
}

###########################################
## MC - PRICING
###########################################
heston.pather <- function(S0,T,delta,rd,rf, kappa,theta,sigma,nu0,rho){
	# Feller Condition
	if(sigma^2 >= 2*kappa*theta) stop("Variance can get zero!")
	n <- T/delta+1
	S <- numeric(n)
	v <- numeric(n)
	S[1] <- S0
	v[1] <- nu0
	X1 <- rnorm(n)
	X2 <- rnorm(n)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	for( i in 2:n){
		S[i] <- S[i-1] + (rd-rf)*S[i-1]*delta + sqrt(abs(v[i-1]))*S[i-1]*sqrt(delta)*W1[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta +sigma*sqrt(abs(v[i-1]))*sqrt(delta)*W2[i]
	}
	S
}
gdiff.pather <- function(S0,T,delta,rd,rf,kappa,theta,sigma,nu0,rho){
	n <- T/delta+1
	S <- numeric(n)
	v <- numeric(n)
	S[1] <- S0
	v[1] <- nu0
	X1 <- rnorm(n)
	X2 <- rnorm(n)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	for( i in 2:n){
		S[i] <- S[i-1] + (rd-rf)*S[i-1]*delta + sqrt(abs(v[i-1]))*S[i-1]*sqrt(delta)*W1[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta + sigma*abs(v[i-1])*sqrt(delta)*W2[i]
	}
	S
}
bates.pather <- function(S0,T,delta,rd,rf,kappa,theta,sigma,nu0,rho,lambda,muJ,vJ){
	# Feller Condition
	if(sigma^2 >= 2*kappa*theta) stop("Variance can get zero!")
	n <- T/delta+1
	S <- numeric(n)
	v <- numeric(n)
	S[1] <- S0
	v[1] <- nu0
	X1 <- rnorm(n)
	X2 <- rnorm(n)
	W3 <- rnorm(n)
	dq <- rbinom(n, 1, lambda*delta)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	for( i in 2:n){
		S[i] <- S[i-1] + ((rd-rf)-lambda*muJ)*S[i-1]*delta + sqrt(abs(v[i-1]))*S[i-1]*sqrt(delta)*W1[i] + (muJ+sqrt(vJ)*W3[i])*S[i-1]*dq[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta + sigma*sqrt(abs(v[i-1]))*sqrt(delta)*W2[i]
	}
	S
}
heston.pather2 <- function(S0,T,delta,rd,rf, kappa,theta,sigma,nu0,rho){
	# Feller Condition
	if(sigma^2 >= 2*kappa*theta) stop("Variance can get zero!")
	n <- T/delta+1
	S <- numeric(n)
	v <- numeric(n)
	S[1] <- S0
	v[1] <- nu0
	X1 <- rnorm(n)
	X2 <- rnorm(n)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	for( i in 2:n){
		S[i] <- S[i-1] + (rd-rf)*S[i-1]*delta + sqrt(abs(v[i-1]))*S[i-1]*sqrt(delta)*W1[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta +sigma*sqrt(abs(v[i-1]))*sqrt(delta)*W2[i]
	}
	list(S,v)
}
heston2d.pather <- function(S0,T,delta,rd,rf,params,randoms,last=TRUE){
	theta = params[1]
	sigma = params[2]
	rho = params[3]
	kappa = params[4]
	nu0 = params[5]
	theta2 = params[6]
	sigma2 = params[7]
	rho2 = params[8]
	kappa2 = params[9]
	nu2 = params[10]
	# Feller Condition
	#if(sigma^2 >= 2*kappa*theta) stop("Variance can get zero!")
	#if(sigma2^2 >= 2*kappa2*theta2) stop("Variance can get zero!")
	n <- T/delta+1
	s <- numeric(n)
	v <- numeric(n)
	v2 <- numeric(n)
	s[1] <- log(S0)
	v[1] <- nu0
	v2[1] <- nu2
	sb=s;vb=v;v2b=v2
	X1 = randoms[1:n]
	X2 = randoms[(n+1):(2*n)]
	X3 = randoms[(2*n+1):(3*n)]
	X4 = randoms[(3*n+1):(4*n)]
	#X1 <- rnorm(n)
	#X2 <- rnorm(n)
	#X3 <- rnorm(n)
	#X4 <- rnorm(n)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	W3 <- X3
	W4 <- rho*X3 + sqrt(1-rho2^2)*X4
	for( i in 2:n){
		s[i] <- s[i-1] + (rd-rf)*delta + sqrt((v[i-1]))*sqrt(delta)*W1[i] + sqrt((v2[i-1]))*sqrt(delta)*W3[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta +sigma*sqrt((v[i-1]))*sqrt(delta)*W2[i] + delta*(W2[i]^2-1)*sigma^2/4
		v2[i] <- v2[i-1] + kappa2*(theta2-v2[i-1])*delta +sigma2*sqrt((v2[i-1]))*sqrt(delta)*W4[i] + delta*(W2[i]^2-1)*sigma^2/4
		sb[i] <- sb[i-1] + (rd-rf)*delta - sqrt((vb[i-1]))*sqrt(delta)*W1[i] - sqrt((v2b[i-1]))*sqrt(delta)*W3[i]
		vb[i] <- vb[i-1] + kappa*(theta-vb[i-1])*delta - sigma*sqrt((vb[i-1]))*sqrt(delta)*W2[i] + delta*(W2[i]^2-1)*sigma^2/4
		v2b[i] <- v2b[i-1] + kappa2*(theta2-v2b[i-1])*delta - sigma2*sqrt((v2b[i-1]))*sqrt(delta)*W4[i] + delta*(W2[i]^2-1)*sigma^2/4
	}
	if(last==TRUE) exp(c(s[n],sb[n]))
	else exp(cbind(s,sb))
}
ouou.pather <- function(S0,T,delta,rd,rf,params,randoms,last=TRUE){
	theta = params[1]
	sigma = params[2]
	rho = params[3]
	kappa = params[4]
	nu0 = params[5]
	theta2 = params[6]
	sigma2 = params[7]
	rho2 = params[8]
	kappa2 = params[9]
	nu2 = params[10]
	# Feller Condition
	#if(sigma^2 >= 2*kappa*theta) stop("Variance can get zero!")
	#if(sigma2^2 >= 2*kappa2*theta2) stop("Variance can get zero!")
	n <- T/delta+1
	s <- numeric(n)
	v <- numeric(n)
	v2 <- numeric(n)
	s[1] <- log(S0)
	v[1] <- nu0
	v2[1] <- nu2
	sb=s;vb=v;v2b=v2
	X1 = randoms[1:n]
	X2 = randoms[(n+1):(2*n)]
	X3 = randoms[(2*n+1):(3*n)]
	X4 = randoms[(3*n+1):(4*n)]
	#X1 <- rnorm(n)
	#X2 <- rnorm(n)
	#X3 <- rnorm(n)
	#X4 <- rnorm(n)
	W1 <- X1
	W2 <- rho*X1 + sqrt(1-rho^2)*X2
	W3 <- X3
	W4 <- rho*X3 + sqrt(1-rho2^2)*X4
	for( i in 2:n){
		s[i] <- s[i-1] + (rd-rf)*delta + (v[i-1])*sqrt(delta)*W1[i] + (v2[i-1])*sqrt(delta)*W3[i]
		v[i] <- v[i-1] + kappa*(theta-v[i-1])*delta +sigma*sqrt(delta)*W2[i] + delta*(W2[i]^2-1)*sigma^2/4
		v2[i] <- v2[i-1] + kappa2*(theta2-v2[i-1])*delta +sigma2*sqrt(delta)*W4[i] + delta*(W2[i]^2-1)*sigma^2/4
		sb[i] <- sb[i-1] + (rd-rf)*delta - (vb[i-1])*sqrt(delta)*W1[i] - (v2b[i-1])*sqrt(delta)*W3[i]
		vb[i] <- vb[i-1] + kappa*(theta-vb[i-1])*delta - sigma*sqrt(delta)*W2[i] + delta*(W2[i]^2-1)*sigma^2/4
		v2b[i] <- v2b[i-1] + kappa2*(theta2-v2b[i-1])*delta - sigma2*sqrt(delta)*W4[i] + delta*(W2[i]^2-1)*sigma^2/4
	}
	if(last==TRUE) exp(c(s[n],sb[n]))
	else exp(cbind(s,sb))
}
heston2d.pricer<- function(S0,T,delta,rd,rf,params,K){
	m=20000
	n <- (T/delta+1)*4
	randoms = matrix( rnorm(n*m,mean=0,sd=1), m, n) 
	S = as.numeric(apply(randoms,1,function(x) heston2d.pather(S0,T,delta,rd,rf,params,x)))
	all=sum(pmax(S-K,0))
	all/m/2*exp(-rd*T)
}
ouou.pricer<- function(S0,T,delta,rd,rf,params,K){
	m=20000
	n <- (T/delta+1)*4
	randoms = matrix( rnorm(n*m,mean=0,sd=1), m, n) 
	S = as.numeric(apply(randoms,1,function(x) ouou.pather(S0,T,delta,rd,rf,params,x)))
	all=sum(pmax(S-K,0))
	all/m/2*exp(-rd*T)
}
ouou2.pricer<- function(S0,T,delta,rd,rf,params,K){
	sum=0
	m=100000
	for(i in 1:m) {
		S=ouou.pather(S0,T,delta,rd,rf,params)
		sum=sum+max(S-K,0)
	}
	sum/m*exp(-rd*T)
}
###########################################
## OTHERS
###########################################
cdf.price<-function(S,K,T,returns,rd=0,rf=0,type=1,upper=7) {
	xy=density(log(S)+(returns+(rd-rf)/252-mean(returns))*sqrt(252*T))
	f<- approxfun(xy[1]$x, xy[2]$y, yleft=0, yright=0)
	if (type==1) 
		exp(-rd*T)*integrate(function(u) (exp((rd-rf)*T)*exp(u)-K)*f(u),log(K),upper)$value
	else
		exp(-rd*T)*integrate(function(u) (exp((rd-rf)*T)*exp(u)-K)*f(u),log(K),upper)$value - S + K*exp(-rd*T)
}
getK<-function(F,vol,tau,rf) {
		K10p=F*exp(vol[,1]^2*tau/2 + vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		K25p=F*exp(vol[,2]^2*tau/2 + vol[,2]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		Katmf=F*exp(vol[,3]^2*tau/2)
		K25c=F*exp(vol[,4]^2*tau/2 - vol[,4]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		K10c=F*exp(vol[,5]^2*tau/2 - vol[,5]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		cbind(K10p,K25p,Katmf,K25c,K10c)
}
data.load1<-function(datatable,i){
		vol<<-t(matrix(as.numeric(datatable[i,2:31]/100),nrow = 5,ncol = 6))
		S<<-as.numeric(datatable[i,32])
		if (pair=="EURUSD") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			ois=as.numeric(datatable[i,39:44]/100)
			fois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="GBPUSD") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			ois=as.numeric(datatable[i,39:44]/100)
			fois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="USDJPY") {
			fwd=as.numeric(datatable[i,33:38]/10^2)
			fois=as.numeric(datatable[i,39:44]/100)
			ois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="USDPLN") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			fois=as.numeric(datatable[i,39:44]/100)
			ois=as.numeric(datatable[i,45:50]/100)
		}
		#rd <<- -log(1/(1+ois*tau))/tau
		#rf <<- -log((1+fwd/S)*exp(-rd*tau))/tau
		rf <<- -log(1/(1+fois*tau))/tau
		rd <<- log((1+fwd/S)*exp(rf*tau))/tau
		#spot+fwd1m=spot*(1+usd1md*yf)/(1+eur1md*yf)
		#spot+fwd1m=spot*exp((usd1md-eur1md)*yf)
		#eur1md=(spot*(1+usd1md*yf)/(spot+fwd1m)-1)/yf
		#eur1md=-log((1+fwd1m/spot)*exp(-usd1md*yf))/yf
		Fwd<<-S*exp((rd-rf)*tau) #F2=S+fwd
		K<<-getK(Fwd,vol,tau,rf)
}
data.load<-function(datatable,pair){
	n=dim(datatable)[1]
	for (i in 1:n) {
		vol=t(matrix(as.numeric(datatable[i,2:31]/100),nrow = 5,ncol = 6))
		S=as.numeric(datatable[i,32])
		if (pair=="EURUSD") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			ois=as.numeric(datatable[i,39:44]/100)
			fois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="GBPUSD") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			ois=as.numeric(datatable[i,39:44]/100)
			fois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="USDJPY") {
			fwd=as.numeric(datatable[i,33:38]/10^2)
			fois=as.numeric(datatable[i,39:44]/100)
			ois=as.numeric(datatable[i,45:50]/100)
		}
		if (pair=="USDPLN") {
			fwd=as.numeric(datatable[i,33:38]/10^4)
			fois=as.numeric(datatable[i,39:44]/100)
			#ois=as.numeric(datatable[i,45:50]/100)
		}
		#rd=-log(1/(1+ois*tau))/tau
		#rf=-log((1+fwd/S)*exp(-rd*tau))/tau
		rf=-log(1/(1+fois*tau))/tau
		rd=log((1+fwd/S)*exp(rf*tau))/tau
		#spot+fwd1m=spot*(1+usd1md*yf)/(1+eur1md*yf)
		#spot+fwd1m=spot*exp((usd1md-eur1md)*yf)
		#eur1md=(spot*(1+usd1md*yf)/(spot+fwd1m)-1)/yf
		#eur1md=-log((1+fwd1m/spot)*exp(-usd1md*yf))/yf
		F=S*exp((rd-rf)*tau) #F2=S+fwd
		K=getK(F,vol,tau,rf)
		#Kb=cbind(K[,1],(K[,1]+K[,2])/2,K[,2],(K[,2]+K[,3])/2,K[,3],(K[,3]+K[,4])/2,K[,4],(K[,4]+K[,5])/2,K[,5])
		#Kext=cbind(3*K[,1]-2*K[,2],K,3*K[,5]-2*K[,4])
		#volbig=cbind(vol[,1],(vol[,1]+vol[,2])/2,vol[,2],(vol[,2]+vol[,3])/2,vol[,3],(vol[,3]+vol[,4])/2,vol[,4],(vol[,4]+vol[,5])/2,vol[,5])
		#volext=cbind(3*vol[,1]-2*vol[,2],vol,3*vol[,5]-2*vol[,4])
		params2=smileRegLog(S,K[1,],vol[1,])
		S_Leif[i]<<-as.numeric(params2[1])
		C_Leif[i]<<-as.numeric(params2[2]*2)
		params=smileRegVec(S,K,vol,tau)
		SS[i,]<<-as.numeric(params[,1])
		CC[i,]<<-as.numeric(params[,2]*2)
		dK=cbind(K[,2]-K[,1],(K[,3]-K[,1])/2,(K[,4]-K[,2])/2,(K[,5]-K[,3])/2,K[,5]-K[,4])
		#dKext=cbind(Kext[,2]-Kext[,1],(Kext[,3]-Kext[,1])/2,(Kext[,4]-Kext[,2])/2,(Kext[,5]-Kext[,3])/2,(Kext[,6]-Kext[,4])/2,(Kext[,7]-Kext[,5])/2,Kext[,7]-Kext[,6])
		#dKb=cbind(Kb[,2]-Kb[,1],(Kb[,3]-Kb[,1])/2,(Kb[,4]-Kb[,2])/2,(Kb[,5]-Kb[,3])/2,(Kb[,6]-Kb[,4])/2,(Kb[,7]-Kb[,5])/2,(Kb[,8]-Kb[,6])/2,(Kb[,9]-Kb[,7])/2,Kb[,9]-Kb[,8])
		GKp=BSMFX(S,K[,1:3],vol[,1:3],rd,rf,0,tau,type=-1)
		GKc=BSMFX(S,K[,3:5],vol[,3:5],rd,rf,0,tau,type=1)
		GK=cbind(GKp[,-3],(GKp[,3]+GKc[,1])/2,GKc[,-1])
		#GKpe=BSMFX(S,Kext[,1:4],volext[,1:4],rd,rf,0,tau,type=-1)
		#GKce=BSMFX(S,Kext[,4:7],volext[,4:7],rd,rf,0,tau,type=1)
		#GKext=cbind(GKpe[,-4],(GKpe[,4]+GKce[,1])/2,GKce[,-1])
		#GKpb=BSMFX(S,Kb[,1:5],volbig[,1:5],rd,rf,0,tau,type=-1)
		#GKcb=BSMFX(S,Kb[,5:9],volbig[,5:9],rd,rf,0,tau,type=1)
		#GKb=cbind(GKpb[,-5],(GKpb[,5]+GKcb[,1])/2,GKcb[,-1])	
		Vix=sqrt(2/tau*exp((rd-rf)*tau)*rowSums(dK*GK/K^2)-1/tau*(F/K[,3]-1)^2)
		#VIX2=(1/tau*E2-(1/tau*E1^2))^(1/2)
		#VIX3=sqrt(-2/tau*P1)
		#bE1=(-1/1*exp((rd-rf)*tau)*rowSums(dK*GK/K^2) - 1/1*(1+log(F/Katmf)-F/Katmf))
		#bE2=(2/1*exp((rd-rf)*tau)*rowSums((1-log(K/F))*dK*GK/K^2) + 1/1*2*log(Katmf/F)*(F/Katmf-1)+1/2*log(Katmf/F)^2)
		#bE3=3/1*exp((rd-rf)*tau)*rowSums((2*log(K/F)-log(K/F)^2)*dK*GK/K^2) + 1/1*3*log(Katmf/F)^2*(1/3*log(Katmf/F)-1+F/Katmf)
		#bE4=4/1*exp((rd-rf)*tau)*rowSums((3*log(K/F)^2-log(K/F)^3)*dK*GK/K^2) + 1/1*3*log(Katmf/F)^2*(1/3*log(Katmf/F)-1+F/Katmf)
		#E1=-1/1*exp((rd)*tau)*rowSums(dK*GK/K^2)
		E2=2*exp((rd)*tau)*rowSums((1-log(K/F))*dK*GK/K^2)
		#E2e=2*exp((rd)*tau)*rowSums((1-log(Kext/F))*dKext*GKext/Kext^2)
		#E2b=2*exp((rd)*tau)*rowSums((1-log(Kb/F))*dKb*GKb/Kb^2)
		E3=3*exp((rd)*tau)*rowSums((2*log(K/F)-log(K/F)^2)*dK*GK/K^2)
		#E3e=3*exp((rd)*tau)*rowSums((2*log(Kext/F)-log(Kext/F)^2)*dKext*GKext/Kext^2)
		#E3b=3*exp((rd)*tau)*rowSums((2*log(Kb/F)-log(Kb/F)^2)*dKb*GKb/Kb^2)
		E4=4*exp((rd)*tau)*rowSums((3*log(K/F)^2-log(K/F)^3)*dK*GK/K^2) 
		#E4e=4*exp((rd)*tau)*rowSums((3*log(Kext/F)^2-log(Kext/F)^3)*dKext*GKext/Kext^2) 
		E5=5*exp((rd)*tau)*rowSums((4*log(K/F)^3-log(K/F)^4)*dK*GK/K^2)
		E6=6*exp((rd)*tau)*rowSums((5*log(K/F)^4-log(K/F)^5)*dK*GK/K^2)
		u = exp((rd-rf)*tau) -1 -1/2*E2-1/6*E3-1/24*E4
		E1 = u 
		#E1e = exp((rd-rf)*tau) -1 -1/2*E2e-1/6*E3e-1/24*E4e
		#P3[i]<<-mean(E3)
		mu[i,]<<-u
		P1[i,]<<-E1
		P2[i,]<<-E2
		P3[i,]<<-E3
		P4[i,]<<-E4
		P5[i,]<<-E5
		P6[i,]<<-E6
		Q1[i,]<<- E1 - (rd-rf)*tau
		Q2[i,]<<- E2 - 2*E1*(rd-rf)*tau + ((rd-rf)*tau)^2
		Q3[i,]<<- E3 - 3*E2*(rd-rf)*tau + 3*E1*((rd-rf)*tau)^2 - ((rd-rf)*tau)^3
		Q4[i,]<<- E4 - 4*E3*(rd-rf)*tau + 6*E2*((rd-rf)*tau)^2 - 4*E1*((rd-rf)*tau)^3 + ((rd-rf)*tau)^4
		VIX[i,]<<-Vix
		#VIXm[i]<<-mean(Vix)
		paramsZhang=smileRegVecZhang(F,K,vol,tau,Vix)
		#level[i,]<<-as.numeric(paramsZhang[,1])
		slope <- as.numeric(paramsZhang[,1])#/as.numeric(paramsZhang[,1])
		curv <- as.numeric(paramsZhang[,2])#/as.numeric(paramsZhang[,1])
		skew[i,] = (vol[,5]-vol[,1])/vol[,3]
		kurt[i,] = (vol[,5]+vol[,1]-2*vol[,3])/2/vol[,3]
		#n0=(1-l2/24)*sigma; #n1=l1/6; #n2=1/24*(1-l2/16)=1/24-l2/16/24
		l1 <<- 6*slope;
		l2 <<- (1/24-curv)*24*16; #l2=24*curv ?
	}
}
int2date <- function(x) {
	if (x<10) paste("0",x,sep="")
	else paste(x,sep="")
}
volofvol <- function(x,m,kappa=0) {
	n=length(x)
	varr=(x^2)
	#if (kappa!=0) e = kappa/(exp(kappa*tau1)-1)
	#else e=1
	varret=c(0,varr[-1]-varr[-n])
	#sqrt(EMA(varret^2,m)*252/(x^2*e))
	sqrt(EMA(varret^2,m)*252/(x^2))
}
hcorr <- function(x,y,m,kappa=0) {
	varret=c(0,y[-1]-y[-n])
	logret=c(0,x[-1]-x[-n])
	EMA(logret*varret,m)/sqrt(EMA(logret^2,m))/sqrt(EMA(varret^2,m))
	#*VIX[,1]
}
volatility <- function(x,k) {
	n=length(x)
	logspot=log(x)
	logret=c(0,logspot[-1]-logspot[-n])
	sqrt(EMA(logret^2,k)*252)
}
TurnbullWakemanAsian<-function(S,SA,x,T,t2,r,B,v,type) {
    # // SA= Realized average so far
    # // t1 = Time to start of average period in years
    # // T =  Time to maturity in years of option  T
    # // T2 = Original time in average period in years, constant over life of option   
    #//tau: reminding time of average period

    t1 = max(0, T - t2)
    tau = t2 - T
   
    if (B==0) {
		M1 = 1
    }
    else {
        M1 = (exp(B * T) - exp(B * t1)) / (B * (T - t1))
    }
    #//Take into account when option wil  be exercised
    if (tau > 0) {
        if (t2 / T * x - tau / T * SA < 0) {
            if (type == 1) {
                SA = SA * (t2 - T) / t2 + S * M1 * T / t2 
				 #//expected average at maturity
                max(0, SA - x) * exp(-r * T)
			}
            else {
                0
            }
            quit
        }
    }
	if (B == 0) {  #//   Extended to hold for options on futures 16 May 1999 Espen Haug
       M2 = 2 * exp(v * v * T) / (v ^ 4 * (T - t1) ^ 2) - 2 * exp(v * v * t1) * (1 + v * v * (T - t1)) / (v ^ 4 * (T - t1) ^ 2)
	}
	else {
        M2 = 2 * exp((2 * B + v * v) * T) / ((B + v * v) * (2 * B + v * v) * (T - t1) ^ 2) + 2 * exp((2 * B + v * v) * t1) / (B * (T - t1) ^ 2) * (1 / (2 * B + v * v) - exp(B * (T - t1)) / (B + v * v))
    }
    bA = log(M1) / T
    vA = sqrt(log(M2) / T - 2 * bA)
    if (tau > 0) {
		x = t2 / T * x - tau / T * SA
        GBSM(S, x, T, r, bA, vA, type) * T / t2
	}
    else {
        GBSM(S, x, T, r, bA, vA, type)
    }
}
GBSM <- function(S,K,tau,r,b,vol,type=1){
	d1 <- (log(S/K) + (b + vol^2/2)*tau)/(vol*sqrt(tau))
	d2 <- d1 - vol*sqrt(tau)
	type * (exp((b-r)*tau)*S*pnorm(type*d1)-exp(-r*tau)*K*pnorm(type*d2))
}
VORST <- function(S,E,vol,r,q,tau,n,type=1) {
	h <- tau/n
	mu <- log(S) + (r - q - vol^2/2)*(tau+h)/2
	sig <- vol*sqrt(h*(2*n+1)*(n+1)/6/n)
	eg <- exp(mu + sig^2/2)
	ea <- S/n*exp((r-q)*h)*(1-exp((r-q)*n*h))/(1-exp((r-q)*h))
	GEOM_BS(S,E-(ea-eg),vol,r,q,tau,n,type)
}
GEOM_BS <- function(S,E,vol,r,q,tau,n,type=1){
	h <- tau/n
	mu <- log(S) + (r - q - vol^2/2)*(tau+h)/2
	sig <- vol*sqrt(h*(2*n+1)*(n+1)/6/n)
	d1 <- (mu - log(E) + sig^2)/sig
	d2 <- d1 - sig
	type * exp(-r*tau)*(exp(mu+sig^2/2)*pnorm(type*d1)-E*pnorm(type*d2))
}
rowMedians <- function(X) {
	n = dim(X)[1]
	y = numeric(n)
	for (i in 1:n) {
		y[i] = median(X[i,])
	}
	y
}
colMedians <- function(X) {
	n = dim(X)[2]
	y = numeric(n)
	for (i in 1:n) {
		y[i] = median(X[,i])
	}
	y
}
int2date <- function(x) {
	if (x<10) paste("0",x,sep="")
	else paste(x,sep="")
}
plot2 <- function(x,y,a=0,b=length(x),ylim=0) {
	if (ylim==0) ylim=c(min(c(x[a:b],y[a:b]),na.rm=TRUE),max(c(x[a:b],y[a:b]),na.rm=TRUE))
	plot(ts(x),ylim=ylim,xlim=c(a,b))
	lines(y,col=2)
}
trapz <- function(x,y) {
	idx = 2:length(x)
	return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1]) / 2))
}
plotLines <- function(y,x,j) {
	plot(y=y[j,],x=y[j,]); lines(y=y[j,],x=x[j,])
	for (i in 1:5) 
	{
		points(y=y[i,],x=K[i,],col=i+1); lines(y=y[i,],x=x[i,],col=i+1)
	}
}
simpson <- function(x,y) {
	sum=0
	for(i in seq(from=3, to=length(x), by=2))
		sum = sum + (x[i] - x[i-2]) * (y[i] + y[i-2] + 4*y[i-1]) / 6	
	sum
}
linspacing <- function(a,b,spacing) {
	N = floor((b-a)/spacing)+1
	a+(0:N)*spacing
}
hestonExpansion <- function(F,K,kappa,theta,sigma,rho,lambda,nu0,rd,t,T,type=1){}
last <- function(x) {
  x[length(x)]
}
EMA <-function (x, n) {
  #ratio = 1/n
  ratio= 2/(n+1)
  x[1:n]= mean(x[1:n])
  #x = coredata(x_xts)
  if(is.na(x[1])) {
    s = which(is.na(x))
    s = s[length(s)]
    y = c(rep(NA,s),filter(x[-(1:s)] * ratio, 1 - ratio, "recursive", init = x[s+1]))
  }
  else {
    y = c(filter(x * ratio, 1 - ratio, "recursive", init = x[1]))
  }
  #xts(y, order.by=index(x_xts))
  y[1:(n-1)]=NA
  y
}
summarypy <- function(x){
  tmp=summary(x)
  y=c(tmp[4],sd(x,na.rm=T),tmp[1:3],tmp[5:6])
  names(y)[2]="S. D."
  y
}
