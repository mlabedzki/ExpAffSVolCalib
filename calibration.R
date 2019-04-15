#################################
# Modified function to optimize #
#################################
feller <- function(params) {
	theta = params[,1]
	sigma = params[,2]
	kappa = params[,4]
	as.numeric(2*kappa*theta-sigma^2)
}
  #Implementation here relies on standard optim() or independent Nelder Mead implementation in C++
	#Outside optim(), but still in the default-install collection, nlm() and nlminb() contain optimization.
	#Following can be also used:
	#fminsearch(fun = NULL, x0 = NULL, options = NULL, verbose=FALSE)
	#from library(neldermead)
	#OR
	#library(adagio)
	#nelmin(fn, x0, tol = 1e-10, maxfeval = 1e4, step = rep(1.0, length(x0)))
	#nelminb(fn, x0, lower, upper, tol = 1e-10, maxfeval = 10000, step = rep(1, length(x0)))
	#OR
	#library("DEoptim")
	#DEoptim(fn=E4,nu0=nu0,lower=c(-5,-5),upper=c(5,5),control=list(storepopfrom=1, trace=FALSE))
  #OR
	#Levenberg–Marquardt
	#also known as the damped least-squares (DLS) methodmethod,
	#The LMA interpolates between the Gauss–Newton algorithm (GNA) and the method of gradient descent
	#library(minpack.lm)
	#nlsLM(yeps ~ p1 / (1 + exp(p2 - x)) * exp(p4 * x), start=list(p1=10,p2=18,p4=-.03))
opter <- function(pars,modType,opType,a,b,i,pair,costfun="MSE",optmet="Nelder-Mead",parOpt1=NA,parOpt2=NA,parOpt3=NA) {
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
		ois=as.numeric(datatable[i,45:50]/100)
	}
	#rd=-log(1/(1+ois*tau))/tau
	#rf=-log((1+fwd/S)*exp(-rd*tau))/tau
	rf=-log(1/(1+fois*tau))/tau
	rd=log((1+fwd/S)*exp(rf*tau))/tau
	F=fwd+S
	K10p=F*exp(vol[,1]^2*tau/2)*exp(vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
	K25p=F*exp(vol[,2]^2*tau/2)*exp(vol[,2]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
	Katmf=F*exp(vol[,3]^2*tau/2)
	K25c=F*exp(vol[,4]^2*tau/2)*exp(-vol[,4]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
	K10c=F*exp(vol[,5]^2*tau/2)*exp(-vol[,5]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
	K=cbind(K10p,K25p,Katmf,K25c,K10c)
	#K <- impStrike(data=vol,t0=0,t=tau,S=S,rd=rd,rf=rf)
	n=length(tau)
	mktPrice=matrix(0,n,5)
	vega=matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			d1 <- (log(S/K[i,j]) + (rd[i] - rf[i] + vol[i,j]^2/2)*tau[i])/(vol[i,j]*sqrt(tau[i]))
			d2 <- d1 - vol[i,j]*sqrt(tau[i])
			vega[i,j] <- K[i,j]*exp(-rd[i]*tau[i])*sqrt(tau[i])*dnorm(d2)
			if (costfun=="vol")
				mktPrice[i,j] = vol[i,j]
			else
				mktPrice[i,j] = 1 * (exp(-rf[i]*tau[i])*S*pnorm(1*d1)-exp(-rd[i]*tau[i])*K[i,j]*pnorm(1*d2))/vega[i,j]
		}
	}
	if (modType=="Heston") {theta=pars[1];sigma=pars[2];rho=pars[3];kappa=pars[4];nu0=pars[5]}
	if (modType=="SZ") {theta=pars[1];sigma=pars[2];rho=pars[3];kappa=pars[4];nu0=pars[5]}
	if (modType=="Bates") {theta=pars[1];sigma=pars[2];rho=pars[3];kappa=pars[4];nu0=pars[5];lambda=pars[6];muJ=pars[7];vJ=pars[8]}
	if (modType=="Heston2d") {theta1=pars[1];sigma1=pars[2];rho1=pars[3];kappa1=pars[4];nu1=pars[5];theta2=pars[6];sigma2=pars[7];rho2=pars[8];kappa2=pars[9];nu2=pars[10]}
	if (modType=="SZ2d") {theta1=pars[1];sigma1=pars[2];rho1=pars[3];kappa1=pars[4];nu1=pars[5];theta2=pars[6];sigma2=pars[7];rho2=pars[8];kappa2=pars[9];nu2=pars[10]}
	if (modType=="HHW") {theta=pars[1];sigma=pars[2];rho=pars[3];kappa=pars[4];nu0=pars[5];eta=pars[6];rhoxr=pars[7]}
	
	Ts=Sys.time()
	if (opType=="E5") {
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa),log(nu0)),E5,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun,method=optmet)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]))
	}
	if (opType=="E4") {
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa)),E4,nu0=nu0,data=vol[a:b,],S=S,rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),nu0)
	}
	if (opType=="E2") {
		opt=optim(c(log(sigma),atanh(rho)),E2,theta=theta,kappa=kappa,nu0=nu0,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(theta,exp(opt$par[1]),tanh(opt$par[2]),kappa,nu0)
	}
	if (opType=="E3my") {
		opt=optim(c(log(theta),log(sigma),log(kappa)),E3my,nu0=nu0,P3=rho,data=vol[a:b,],S=S,rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),exp(opt$par[2]),rhoFixP3(nu0,rho,mean(rd-rf),exp(opt$par[1]),exp(opt$par[2])),exp(opt$par[3]),nu0)
	}
	if (opType=="E3guil") {
		opt=optim(c(log(sigma),atanh(rho),log(kappa)),E3guil,theta=theta,nu0=nu0,data=vol[a:b,],S=S,rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(theta,exp(opt$par[1]),tanh(opt$par[2]),exp(opt$par[3]),nu0)
	}
	if (opType=="E3smile") {		
		#init=optim(c(1,1),smileError,data=vol[a:b,],S=S,K=K,t=0,T=tau[a:b])$par
		params=smileReg(S,K,vol,tau,1);init=c(params[1],params[2]*2)
		#opt=optim(c(log(theta),log(kappa)),E3smile,nu0=nu0,init=init,data=vol[a:b,],S=S,rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		opt=optim(c(log(theta),log(kappa)),E3smile,nu0=nu0,init=init,data=vol[a:b,],S=S,rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),smile2omega(sqrt(nu0),init[1],init[2]),smile2rho(sqrt(nu0),init[1],init[2]),exp(opt$par[2]),nu0)
	}
	if (opType=="EB8") {
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa),log(nu0),log(lambda),muJ,log(vJ)),EB8,data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),opt$par[7],exp(opt$par[8]))
	}
	if (opType=="EB3+5")	{
		init=optim(c(log(lambda),muJ,log(vJ)),EM3,data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa),log(nu0)),EB5,lambda=exp(init$par[1]),muJ=init$par[2],vJ=exp(init$par[3]),data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(init$par[1]),init$par[2],exp(init$par[3]))
		#result=c(exp(init$par[1]),init$par[2],exp(init$par[3]))	
	}
	if (opType=="E2H10") {
		opt=optim(c(log(theta1),log(sigma1),atanh(rho1),log(kappa1),log(nu1),log(theta2),log(sigma2),atanh(rho2),log(kappa2),log(nu2)),E2H10,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),exp(opt$par[7]),tanh(opt$par[8]),exp(opt$par[9]),exp(opt$par[10]))
	}
	if (opType=="E2H10f") {
		opt=optim(c(log(theta1),log(sigma1),atanh(rho1),log(kappa1),log(nu1),log(theta2),log(sigma2),atanh(rho2),log(kappa2),log(nu2)),E2H10f,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),exp(opt$par[7]),tanh(opt$par[8]),exp(opt$par[9]),exp(opt$par[10]))
	}
	if (opType=="E2H8") {
		opt=optim(c(log(theta1),log(sigma1),log(kappa1),log(nu1),log(theta2),log(sigma2),log(kappa2),log(nu2)),EE2H8,rho1=rho1,rho2=rho2,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),rho1,exp(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),rho2,exp(opt$par[7]),exp(opt$par[8]))
	}
	if (opType=="E2H8+2")	{
		opt=optim(c(log(theta1),log(sigma1),atanh(rho1),log(kappa1),log(theta2),log(sigma2),atanh(rho2),log(kappa2)),E2H8F,nu1=nu1,nu2=nu2,data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		vopt=optim(c(log(nu1),log(nu2)),E2H2F,theta1=exp(opt$par[1]),sigma1=exp(opt$par[2]),rho1=tanh(opt$par[3]),kappa1=exp(opt$par[4]),theta2=exp(opt$par[5]),sigma2=exp(opt$par[6]),rho2=tanh(opt$par[7]),kappa2=exp(opt$par[8]),data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		oopt=optim(c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5],opt$par[6],opt$par[7],opt$par[8]),E2H8F,nu1=exp(vopt$par[1]),nu2=exp(vopt$par[2]),data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		opt=optim(c(vopt$par[1],vopt$par[2]),E2H2F,theta1=exp(oopt$par[1]),sigma1=exp(oopt$par[2]),rho1=tanh(oopt$par[3]),kappa1=exp(oopt$par[4]),theta2=exp(oopt$par[5]),sigma2=exp(oopt$par[6]),rho2=tanh(oopt$par[7]),kappa2=exp(oopt$par[8]),data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(oopt$par[1]),exp(oopt$par[2]),tanh(oopt$par[3]),exp(oopt$par[4]),exp(opt$par[1]),exp(oopt$par[5]),exp(oopt$par[6]),tanh(oopt$par[7]),exp(oopt$par[8]),exp(opt$par[2]))
	}
	if (opType=="E2H8F2")	{
		opt=optim(c(qnorm(nu2),log(sigma1),atanh(rho1),log(kappa1),log(sigma2),atanh(rho2),log(kappa2)),E2H8F2,nu0=nu1,theta=theta1,data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(theta1*pnorm(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),nu1*pnorm(opt$par[1]),theta1*(1-pnorm(opt$par[1])),exp(opt$par[5]),tanh(opt$par[6]),exp(opt$par[7]),nu1*(1-pnorm(opt$par[1])))
	}
	if (opType=="EHHW7") {
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa),log(nu0),log(eta),rhoxr),EHHW7,rlambda=rlambda[i],rtheta=rtheta[i],r0=r0[i],data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b])
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),opt$par[7])
	}	
	if (opType=="ESZ5") {
		opt=optim(c(log(theta),log(sigma),atanh(rho),log(kappa),log(nu0)),ESZ5,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]))
	}
	if (opType=="ESZ2") {
		opt=optim(c(log(sigma),atanh(rho)),E2,model=modType,theta=theta,kappa=kappa,nu0=nu0,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(theta,exp(opt$par[1]),tanh(opt$par[2]),kappa,nu0)
	}
	if (opType=="ESZ2D10") {
		opt=optim(c(log(theta1),log(sigma1),atanh(rho1),log(kappa1),log(nu1),log(theta2),log(sigma2),atanh(rho2),log(kappa2),log(nu2)),ESZ2D10,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),exp(opt$par[5]),exp(opt$par[6]),exp(opt$par[7]),tanh(opt$par[8]),exp(opt$par[9]),exp(opt$par[10]))
	}
	if (opType=="E2H7IM") {
		opt=optim(c(log(theta1),log(sigma1),atanh(rho1),log(kappa1),qnorm(nu1),log(theta2),log(kappa2)),E2H7IM,omegas=parOpt1,rhos=parOpt2,nu0=nu2,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),exp(opt$par[2]),tanh(opt$par[3]),exp(opt$par[4]),pnorm(opt$par[5]),exp(opt$par[6]),exp(opt$par[7]))
	}
	if (opType=="E2H2IM") {
		opt=optim(c(log(sigma1),atanh(rho1)),E2H2IM,omegas=parOpt1,rhos=parOpt2,nuR=parOpt3,nu0=nu2,nu0R=nu1,theta1=theta1,theta2=theta2,kappa1=kappa1,kappa2=kappa2,data=mktPrice[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun=costfun)
		result=c(exp(opt$par[1]),tanh(opt$par[2]))
	}
	t=Sys.time()-Ts;
	if (costfun=="MSE") opt_value=sqrt(opt$value)
	else {if (costfun=="MAE") opt_value=opt$value
	else {if (costfun=="MSPE") opt_value=sqrt(opt$value)
	else {if (costfun=="MAPE") opt_value=opt$value
	else opt_value=sqrt(opt$value) }}}
	c(result,opt_value,t)
}
#5 params to opt: kappa, rho, lambda/sigma, nu0, theta
E5old <- function(params,paramsContr,data,S,K,rd,rf,t,T,vega,costfun){
	#if(length(data)<6) data=rbind(data,data)
	if (paramsContr[1]==0) theta <- exp(params[1]) else theta <- paramsContr[1]
	if (paramsContr[2]==0) sigma <- exp(params[2]) else sigma <- paramsContr[2]
	if (paramsContr[3]==0) rho <- tanh(params[3]) else rho <- rhoFix(VIX[2],Skew[2],theta,sigma,rd[2])
	if (paramsContr[4]==0) kappa <- exp(params[4]) else kappa <- paramsContr[4]
	if (paramsContr[5]==0) nu0 <- exp(params[5]) else nu0 <- paramsContr[5]
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	#vol <- matrix(0,n,5)
	if(n==1) {
		K=t(as.matrix(K))
		data=t(as.matrix(data))
	}
	for(i in 1:n){
		for(j in 1:5){
			if (is.na(K[i,j])) {
				price[i,j] <- NA
			}
			else {
				price[i,j] <- hestongg(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])/vega[i,j]
				#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
				#price[i,j] <- callCF(cf=cfHeston, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
				#vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
		}
	}
	#if (2*kappa*theta>sigma^2) {
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)}}}
		#sqrt( sum( (vol - data)^2 ) )
	#}
	#else
	#	NA
}

E5 <- function(params,data,S,K,rd,rf,t,T,vega,costfun){
	#if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	nu0 <- exp(params[5])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	Call <- matrix(0,n,5)
	Put <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	u = (0:100)+(1e-16)
	if(n==1) {
		K=t(as.matrix(K))
		data=t(as.matrix(data))
	}
	nu0=nu0[i];theta=theta[i];kappa=kappa[i];rho=rhoLmed[i];sigma=omegaLmed[i];lambda=0
	for(i in 1:n){
		stale <- CFa(-1i, S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i],type=2)
		I1 <- CFa(u-1i,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i],type=2)/(1i*u)/stale
		I2 <- CFa(u,S=S,kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=lambda,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i],type=2)/(1i*u)
		for(j in 1:5){
			if (is.na(K[i,j])) {
				price[i,j] <- NA
			}
			else {
				price[i,j] <- hestonggPre(S=S,K=K[i,j],rd=rd[i],rf=rf[i],t=0,T=T[i],u=u,I1=I1,I2=I2)#/vega[i,j]
				#vol[i,j] <- hestonggVol(S,K[i,j],rd[i],rf[i],t=0,tau[i],u,I2,dtype=-1)
				Call <- hestonggPre(S=S,K=K[i,j],rd=rd[i],rf=rf[i],t=0,T=T[i],u=u,I1=I1,I2=I2)
				Put <- hestonggPre(S=S,K=K[i,j],rd=rd[i],rf=rf[i],t=0,T=T[i],u=u,I1=I1,I2=I2,type=-1)
				#vol[i,j] <- impvol_approx(S,K[i,j],rd[i],rf[i],tau[i],Call,Put)
				#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
				#price[i,j] <- callCF(cf=cfHeston, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
		}
	}
	
	#if (2*kappa*theta>sigma^2) {
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)}}}
		#sqrt( sum( (vol - data)^2 ) )
	#}
	#else
	#	NA
}

E2 <- function(params,theta,kappa,nu0,data,S,K,rd,rf,t,T,vega,costfun){
	#if(length(data)<6) data=rbind(data,data)
	sigma <- exp(params[1])
	rho <- tanh(params[2])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	#vol <- matrix(0,n,5)
	if(n==1) {
		K=t(as.matrix(K))
		data=t(as.matrix(data))
	}
	for(i in 1:n){
		for(j in 1:5){
			if (is.na(K[i,j])) {
				price[i,j] <- NA
			}
			else {
				price[i,j] <- hestongg(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])/vega[i,j]
				#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
				#price[i,j] <- callCF(cf=cfHeston, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
				#vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
		}
	}
	#if (2*kappa*theta>sigma^2) {
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)}}}
		#sqrt( sum( (vol - data)^2 ) )
	#}
	#else
	#	NA
}
E2 <- function(params,theta,kappa,nu0,data,S,K,rd,rf,t,T,vega,costfun,model="Heston"){
	#if(length(data)<6) data=rbind(data,data)
	sigma <- exp(params[1])
	rho <- tanh(params[2])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	#vol <- matrix(0,n,5)
	if(n==1) {
		K=t(as.matrix(K))
		data=t(as.matrix(data))
	}
	for(i in 1:n){
		for(j in 1:5){
			if (is.na(K[i,j])) {
				price[i,j] <- NA
			}
			else {
				if (model=="Heston")
					price[i,j] <- hestongg(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])/vega[i,j]
				else 
					price[i,j] <- sz(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])/vega[i,j]
				#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
				#price[i,j] <- callCF(cf=cfHeston, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
				#vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
		}
	}
	#if (2*kappa*theta>sigma^2) {
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)}}}
		#sqrt( sum( (vol - data)^2 ) )
	#}
	#else
	#	NA
}

###################################
# 1 Fixed parameters optimization #
###################################
E4 <- function(params,nu0,data,S,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- hestonaa(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}

ESZ5 <- function(params,data,S,K,rd,rf,t,T,vega,costfun){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	nu0 <- exp(params[5])
	n <- length(T)
	price <- matrix(0,n,5)
	#vol <- matrix(0,n,5)
	#K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- sz(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])/vega[i,j]
			#vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)}}}
		#sqrt( sum( (vol - data)^2 ) )
}

#######################
# Function to fix Rho #
#######################
rhoFixP3<-function(nu0,P3,mu,theta,sigma){
	min(max((P3-mu*nu0)/(nu0*sigma),-1),1)
}
rhoFixS<-function(nu0,Skew,mu,theta,sigma){
	(Skew*nu0^(3/2)-mu*theta+3*mu*nu0+mu^3)/(nu0*sigma)
}
rhoFix<-function(VIX,Skew,theta,sigma,rd){
	(-rd*theta+Skew*VIX^3+3*rd*VIX^2-2*rd^3)/VIX/sqrt(sigma)
}
rhoFixC<-function(VIX,Skew,theta,sigma,rd){
	(-rd*theta+Skew*VIX^3+3*rd*VIX^2)/VIX/sigma
}
#dV = 
#
#
smile2rho<-function(VIX,SS,CC){
	result=3*VIX*CC+3*VIX*SS+10*SS^2
	w=which(result>0)
	result[w] = 2*SS*result^(-1/2)
	result[!w] = sign(SS)
	result
}
smile2rhom<-function(omega,VIX,SS,CC){
	4*SS*VIX/omega
}
smile2omega<-function(VIX,SS,CC){
	result=3*VIX*CC+3*VIX*SS+10*SS^2
	w=which(result>0)
	result[w] = 2*VIX*result^(1/2)
	result[!w] = 0.1 #2*VIX*(3*VIX*CC+3*VIX*0.01+10*SS^2)^(-1/2)
	result
}
smile2om2<-function(VIX,SS,CC){
	4*VIX^2*(3*VIX*CC+3*VIX*SS+10*SS^2)
}
smile<-function(atm,strike,spot,SS,CC) {
	atm+(strike/spot-1)*SS+(strike/spot-1)^2*CC/2
}
smileError<- function(params,data,S,K,t,T){
	SS <- params[1]
	CC <- params[2]
	n <- length(T)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			if (n>1)
				vol[i,j] <- smile(data[i,3],K[i,j],S,SS,CC)
			else
				vol[i,j] <- smile(data[3],K[j],S,SS,CC)
		}
	}
	sqrt( sum( (vol - data)^2 ) )
}
smileRegLog<-function(S,K,vol){
	sk = as.numeric(log(K/S))
	sk2 = sk^2
	volSpread = as.numeric(vol-vol[3])
	model = lm(volSpread~sk+sk2-1)
	model$coef
}
smileReg<-function(S,K,vol,tau,i){
	sk = (K[i,])/S-1
	sk2 = ((K[i,])/S-1)^2
	volSpread = as.numeric(vol[i,]-vol[i,3])
	model = lm(volSpread~sk+sk2-1)
	model$coef
}
smileRegVec<-function(S,K,vol,tau){
	n <- length(tau)
	coefs <- matrix(0,n,2)
	for(i in 1:n){
		coefs[i,]=smileReg(S,K,vol,tau,i)
	}
	coefs
}
smileRegAll<-function(S,K,vol,tau){
	sk = as.numeric(K[,])/S-1
	sk2 = (as.numeric(K[,])/S-1)^2
	volSpread = as.numeric(vol[,]-vol[,3])
	model = lm(volSpread~sk+sk2-1)
	model$coef
}
smileZhang<-function(level,slope,curv,K,F,tau,sigma) {
	xi = log(K/F)/sqrt(tau)/sigma
	level*(1+xi*slope+xi^2*curv)
}
smileRegZhang<-function(F,K,vol,tau,sigma,i){
	xi = log(K[i,]/F[i])/sqrt(tau[i])/mean(vol[,3])
	xi1 = xi
	xi2 = xi^2
	volSpread = as.numeric(vol[i,]/vol[i,3]-1)
	model = lm(volSpread~xi1+xi2-1)
	model$coef
}
smileRegVecZhang<-function(F,K,vol,tau,sigma){
	n <- length(tau)
	coefs <- matrix(0,n,2)
	for(i in 1:n){
		coefs[i,]=smileRegZhang(F,K,vol,tau,sigma,i)
	}
	coefs
}
lmuvFit <- function(par,S,K,rd,rf,tau,i) {
	fit=rep(0,5)
	for(j in 1:5){
		bs1=BSMFX(S,K[i,j],par[1],rd[i],rf[i],t=0,T=tau[i],type=1)
		bs2=BSMFX(S,K[i,j],par[2],rd[i],rf[i],t=0,T=tau[i],type=1)
		fit[j]=nlminb(0.2, BSMFXsolve, price=par[3]*bs1+(1-par[3])*bs2, S=S,K=K[i,j],rd=rd[i],rf=rf[i],t=0,T=tau[i],type=1, lower=0)$par
	}
	fit
}
lmuvErr <- function(params,data,S,K,rd,rf,tau){
	v1 <- params[1]
	v2 <- params[2]
	l <-params[3]
	#n <- length(T)
	#vol <- matrix(0,n,5)
	n=1
	vol = rep(0,5)
	for(i in 1:n){
		for(j in 1:5){
			if (n<=1) {
				bs1=BSMFX(S,K[i,j],v1,rd[i],rf[i],t=0,tau[i],type=1)
				bs2=BSMFX(S,K[i,j],v2,rd[i],rf[i],t=0,tau[i],type=1)
				vol[j] <- nlminb(0.2, BSMFXsolve, price=l*bs1+(1-l)*bs2, S=S,K=K[i,j],rd=rd[i],rf=rf[i],t=0,T=tau[i],type=1, lower=0)$par
			}
			else
				vol[i,j] <- smile(data[3],K[j],S,SS,CC)
		}
	}
	sqrt( sum( (vol - data)^2 ) )
}

###################################
# 2 Fixed parameters optimization #
###################################
E3smile <- function(params,nu0,init,data,S,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	kappa <- exp(params[2])
	#sigma <- exp(params[3])
	sigma <- smile2omega(sqrt(nu0),init[1],init[2])
	#rho <- tanh(params[3])
	rho <- smile2rho(sqrt(nu0),init[1],init[2])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- hestonaa(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
E3my <- function(params,nu0,P3,data,S,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- rhoFixP3(nu0,P3,mean(rd-rf),theta,sigma)
	kappa <- exp(params[3])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- hestonaa(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
E3guil <- function(params,theta,nu0,data,S,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	sigma <- exp(params[1])
	rho <- tanh(params[2])
	kappa <- exp(params[3])
	lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	K <- impStrike(data=data,t0=0,t=T,S=S,rd=rd,rf=rf)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- hestonaa(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}

ESZ2D10 <- function(params,data,S,K,rd,rf,t,T,vega,costfun){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])
	kappa1 <- exp(params[4])
	nu1 <- exp(params[5])
	theta2 <- exp(params[6])
	sigma2 <- exp(params[7])
	rho2 <- tanh(params[8])
	kappa2 <- exp(params[9])
	nu2 <- exp(params[10])
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- sz2d(S=S,K=K[i,j],theta1,sigma1,rho1,kappa1,nu1,theta2,sigma2,rho2,kappa2,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol")
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			else
				price[i,j] <- price[i,j] / vega[i,j]
		}
	}
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)
		else {if (costfun=="vol") sum( (vol - data)^2 ,na.rm=TRUE) }}}}
}

E2H10 <- function(params,data,S,K,rd,rf,t,T,vega,costfun){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])
	kappa1 <- exp(params[4])
	nu1 <- exp(params[5])
	theta2 <- exp(params[6])
	sigma2 <- exp(params[7])
	rho2 <- tanh(params[8])
	kappa2 <- exp(params[9])
	nu2 <- exp(params[10])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol")
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			else
				price[i,j] <- price[i,j] / vega[i,j]

		}
	}
	if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)
		else {if (costfun=="vol") sum( (vol - data)^2 ,na.rm=TRUE) }}}}
	else
		NA
}
E2H10f <- function(params,data,S,K,rd,rf,t,T,vega,costfun){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])
	kappa1 <- exp(params[4])
	nu1 <- exp(params[5])
	theta2 <- exp(params[6])
	sigma2 <- exp(params[7])
	rho2 <- tanh(params[8])
	kappa2 <- exp(params[9])
	nu2 <- exp(params[10])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol")
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			else
				price[i,j] <- price[i,j] / vega[i,j]

		}
	}
	#if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		if (costfun=="MSE") sum( (price - data)^2 ,na.rm=TRUE)
		else {if (costfun=="MAE") sum( abs(price - data) ,na.rm=TRUE)
		else {if (costfun=="MSPE") sum( (price/data-1)^2 ,na.rm=TRUE)
		else {if (costfun=="MAPE") sum( abs(price/data-1) ,na.rm=TRUE)
		else {if (costfun=="vol") sum( (vol - data)^2 ,na.rm=TRUE) }}}}
	#else
	#	NA
}
E2H8 <- function(params,rho1,rho2,data,S,K,rd,rf,t,T,vega,costfun){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	kappa1 <- exp(params[3])
	nu1 <- exp(params[4])
	theta2 <- exp(params[5])
	sigma2 <- exp(params[6])
	kappa2 <- exp(params[7])
	nu2 <- exp(params[8])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol")
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			else
				price[i,j] <- price[i,j] / vega[i,j]

		}
	}
	if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		if (costfun=="MSE") sum( (price - data)^2 )
		else {if (costfun=="MAE") sum( abs(price - data) )
		else {if (costfun=="MSPE") sum( (price/data-1)^2 )
		else {if (costfun=="MAPE") sum( abs(price/data-1) )
		else {if (costfun=="vol") sum( (vol - data)^2 ) }}}}
	else
		NA
}
E2H8F2 <- function(params,nu0,theta,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	share <- pnorm(params[1])
	theta1 <- theta*share
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])/
	kappa1 <- exp(params[4])
	nu1 <- nu0*share
	theta2 <- theta*(1-share)
	sigma2 <- exp(params[5])
	rho2 <- tanh(params[6])
	kappa2 <- exp(params[7])
	nu2 <- nu0*(1-share)
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
E2H8F <- function(params,nu1,nu2,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])
	kappa1 <- exp(params[4])
	theta2 <- exp(params[5])
	sigma2 <- exp(params[6])
	rho2 <- tanh(params[7])
	kappa2 <- exp(params[8])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	#if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		sqrt( sum( (vol - data)^2 ) )
	#else
	#	NA
}
E2H2F <- function(params,theta1,sigma1,rho1,kappa1,theta2,sigma2,rho2,kappa2,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	nu1 <- exp(params[1])
	nu2 <- exp(params[2])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1,theta1,sigma1,rho1,lambda1=0,nu1,kappa2,theta2,sigma2,rho2,lambda2=0,nu2,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	#if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		sqrt( sum( (vol - data)^2 ) )
	#else
	#	NA
}
EB8 <- function(params,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	nu0 <- exp(params[5])
	lambda <- exp(params[6])
	muJ <- params[7]
	vJ <- exp(params[8])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			#price[i,j] <- heston(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
			price[i,j] <- callCF(cf=cfBates, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, lambda=lambda, muJ=muJ, vJ=vJ, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
EB5 <- function(params,lambda,muJ,vJ,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	nu0 <- exp(params[5])
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			#price[i,j] <- heston(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
			price[i,j] <- callCF(cf=cfBates, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, lambda=lambda, muJ=muJ, vJ=vJ, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
EM3 <- function(params,data,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	lambda <- exp(params[1])
	muJ <- params[2]
	vJ <- exp(params[3])
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			#price[i,j] <- heston(S=S,K=K[i,j], kappa=kappa,theta=theta,sigma=sigma,rho=rho,lambda=0,nu0=nu0,rd=rd[i],rf=rf[i],t=0,T=T[i])
			#price[i,j] <- callHestoncf(S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v0=nu0, vT=theta, rho=rho, k=kappa, sigma=sigma, implVol = FALSE)#$impliedVol
			price[i,j] <- callCF(cf=cfMerton, S=S, X=K[i,j], tau=T[i], r=rd[i], q=rf[i], v=data[i,3]^2, lambda=lambda, muJ=muJ, vJ=vJ, implVol = FALSE, uniroot.control = list(), uniroot.info = FALSE)
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	sqrt( sum( (vol - data)^2 ) )
}

EHHW7 <- function(params,data,rlambda,rtheta,r0,S,K,rd,rf,t,T){
	if(length(data)<6) data=rbind(data,data)
	theta <- exp(params[1])
	sigma <- exp(params[2])
	rho <- tanh(params[3])
	kappa <- exp(params[4])
	nu0 <- exp(params[5])
	eta <- exp(params[6])
	rhoxr <- params[7]
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- hhw(S=S,K=K[i,j],lambda=rlambda,theta=rtheta,r0=r0,eta=eta,kappa=kappa,sigma=theta,nu0=nu0,gamma=sigma,rhoxr=rhoxr,rhoxv=rho,rd=rd[i],rf=rf[i],t=0,T=T[i])
			vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
		}
	}
	if (2*kappa*theta>sigma^2)
		sqrt( sum( (vol - data)^2 ) )
	else
		NA
}
szTSvar<-function(params){
	k=dim(params)[1]
	variance=matrix(0,k,6)
	for(i in 1:k){
		nu0 = params[i,1]
		theta = params[i,2]
		kappa = params[i,3]
		x = params[i,4]^2/2/kappa
		variance[i,] = (exp(-2*kappa*tau)*(exp(2*kappa*tau)*(theta^2*(2*kappa*tau-3)+2*theta*nu0+x*(2*kappa*tau-1)+nu0^2)+4*theta*(theta-nu0)*exp(kappa*tau)-(theta-nu0)^2+x))/(2*kappa*tau)
	}
	variance
}
szTS<-function(params,data,tau,omega){
	nu0 = exp(params[1])
	theta = exp(params[2])
	kappa = exp(params[3])
	x = omega^2/2/kappa		
	k=length(tau)
	variance=rep(0,k)
	for(i in 1:k){
		variance[i] = (exp(-2*kappa*tau[i])*(exp(2*kappa*tau[i])*(theta^2*(2*kappa*tau[i]-3)+2*theta*nu0+x*(2*kappa*tau[i]-1)+nu0^2)+4*theta*(theta-nu0)*exp(kappa*tau[i])-(theta-nu0)^2+x))/(2*kappa*tau[i])
	}
	sum( abs(variance - data) )
}
szTSfit<-function(var,tau,kappa=2,a=2,b=6,omega=NA){
	n = dim(vol)[1]
	TSpar=matrix(0,n,4)
	for(i in 1:n){
		params=c(log(sqrt(var[i,1])),log(sqrt(var[i,6])),log(kappa))
		opt=optim(params,szTS,data=var[i,a:b],tau=tau[a:b],omega=omega[i],method="BFGS")
		TSpar[i,1]=exp(opt$par[1])
		TSpar[i,2]=exp(opt$par[2])
		TSpar[i,3]=exp(opt$par[3])
		TSpar[i,4]=opt$value
	}
	TSpar
}
hestonTS<-function(params,data,tau){
	nu0 = exp(params[1])
	theta = exp(params[2])
	kappa = exp(params[3])
	k=length(tau)
	variance=rep(0,k)
	for(i in 1:k){
		variance[i] <- theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i])
	}
	sqrt( sum( (variance - data)^2 ) )
}
hestonTS2<-function(params,kappa,data,tau){
	nu0 = exp(params[1])
	theta = exp(params[2])
	k=length(tau)
	variance=rep(0,k)
	for(i in 1:k){
		variance[i] <- theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i])
	}
	sqrt( sum( (variance - data)^2 ) )
}
hestonTSfull<-function(params,data,tau,omega,rho){
	if (is.na(omega)) omega=0
	if (is.na(rho)) rho=0
	nu0 = exp(params[1])
	theta = exp(params[2])
	kappa = exp(params[3])
	k=length(tau)
	XY=rep(0,k);Y2=XY
	variance=XY
	for(i in 1:k){
		kt = kappa*tau[i]
		ey2int=((exp(-kt)*(((2*kt-3)*exp(kt)+4*kt)*theta-4*kt*nu0))-((exp(-2*kt))*((2*exp(2*kt)-4*exp(kt)-1)*theta-2*nu0*exp(2*kt)+2*nu0)))/(2*kappa^3)
		exyint = exyint = kappa^(-2)*(exp(-kt)*(2*theta-nu0-kt*(nu0-theta))+(nu0+theta*(kt-2)))
		variance[i] <- theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i]) - omega*rho*exyint + omega^2*ey2int/4
	}
	sqrt( sum( (variance - data)^2 ) )
}
hestonTSfit<-function(vol,tau,last=FALSE,kappa=2,a=2,b=5,full=FALSE,omega=NA,rho=NA){
	n = dim(vol)[1]
	hestonTSpar=matrix(0,n,4)
	for(i in 1:n){
		if (is.na(vol[i,6])) {
			params=c(log(vol[i,1]^2),log(vol[i,3]^2),log(2))
			opt=optim(params,hestonTS,data=vol[i,2]^2,tau=tau[2],method="BFGS");
		}
		else { if (is.na(vol[i,1])) {
			params=c(log(vol[i,2]^2),log(vol[i,6]^2),log(2))
			opt=optim(params,hestonTS,data=vol[i,2]^2,tau=tau[4:5],method="BFGS");		
		}
		else {
			if (i>1 & last) kappa.init = hestonTSpar[i-1,3]
			else kappa.init = kappa
			params=c(log(vol[i,1]^2),log(vol[i,6]^2),log(kappa.init))
			if (full) opt=optim(params,hestonTSfull,omega=omega[i],rho=rho[i],data=vol[i,a:b]^2,tau=tau[a:b])#,method="BFGS")
			else opt=optim(params,hestonTS,data=vol[i,a:b]^2,tau=tau[a:b])#,method="BFGS")
		}}
		hestonTSpar[i,1]=exp(opt$par[1])
		hestonTSpar[i,2]=exp(opt$par[2])
		hestonTSpar[i,3]=exp(opt$par[3])
		hestonTSpar[i,4]=opt$value
	}
	hestonTSpar
}
hestonTSfit2<-function(vol,tau,kappa,a,b){
	n = dim(vol)[1]
	hestonTSpar=matrix(0,n,3)
	for(i in 1:n){
		params=c(log(vol[i,1]^2),log(vol[i,6]^2))
		opt=optim(params,hestonTS2,kappa=kappa,data=vol[i,a:b]^2,tau=tau[a:b],method="BFGS")
		hestonTSpar[i,1]=exp(opt$par[1])
		hestonTSpar[i,2]=exp(opt$par[2])
		hestonTSpar[i,3]=opt$value
	}
	hestonTSpar
}
kappaErr<-function(param,nu0,theta,data,tau){
	kappa = exp(param)
	k=length(tau)
	variance=rep(0,k)
	for(i in 1:k){
		variance[i] <- theta + (nu0-theta) * (1-exp(-kappa*tau[i]))/(kappa*tau[i])
	}
	sqrt( sum( (variance - data)^2 ) )
}
kappaFit<-function(vol,tau){
	n = dim(vol)[1]
	kappa0=rep(0,n)
	for(i in 1:n){
		if (is.na(vol[i,6])) {
			param=log(2)
			opt=optimize(kappaErr,nu0=log(vol[i,1]^2),theta=log(vol[i,6]^2),data=vol[i,1:3]^2,tau=tau[1:3],interval=c(-10,15));
		}
		else {
			param=log(2)
			opt=optimize(kappaErr,nu0=log(vol[i,1]^2),theta=log(vol[i,6]^2),data=vol[i,]^2,tau=tau,interval=c(-10,15));
		}
		kappa0[i]=opt$minimum
	}
	exp(kappa0)
}
kappaReg<-function(variance) {
	n=length(variance)
	k=numeric(n)
	k[1:30] = rep(NA,30)
	for (i in 30:n) {
		y = variance[1:i]
		x = c(NA,variance[1:(i-1)])
		lm1 = lm(y~x)
		k[i] = coef(lm1)[1]/(1-coef(lm1)[2])*252
	}
	k
}
tester1d <- function(paramsPos,F,K,rd,mktPrice,vega,tau,model=1){
  n=length(tau)
  errMatrix <- matrix(0,n,5)
  
    price <- numeric(n*5)
    for(i in 1:n){
      for(j in 1:5){
        ij = (i-1)*5+j
        if (is.na(K[ij])) {
          price[ij] <- NA
          #ivol[ij] <- NA
        }
        else {
          if (model==1)#Heston
            price[ij] <- hestonAttari(paramsPos, F[i], K[ij], rd[i], tau[i])
          else if (model==2)#Schoebel-Zhu
            price[ij] <- szAttari(paramsPos, F[i], K[ij], rd[i], tau[i])
          else if (model==3)#Bates
            price[ij] <- batesAttari(paramsPos, F[i], K[ij], rd[i], tau[i])
          else if (model==4)#OUOU
            price[ij] <- ououAttari(paramsPos, F[i], K[ij], rd[i], tau[i])
          #if (errType=="IV") ivol[ij] <- impVol(data=price[i,j],t0=0,t=tau[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
          price[ij] <- price[ij]/vega[ij]
          errMatrix[i,j] = errMatrix[i,j]  + (price[ij]-mktPrice[ij])^2
          #errMatrix[i,j]  + (ivol[ij] - vol[k,ij])^2
        }
      }
    }
    #if (errType=="IV") err[k]=sqrt( sum( (ivol - vol)^2,na.rm=TRUE) )
    err=sqrt(mean((price - mktPrice)^2,na.rm=TRUE))
  
  list(err,errMatrix)
}
tester <- function(paramsPos,F,K,rd,mktPrice,vega,tau,model=1){
	m=dim(K)[1]
	n=length(tau)
	err=rep(NA,m)
	errMatrix <- matrix(0,n,5)
	for(k in 1:m){
		price <- numeric(n*5)
		for(i in 1:n){
			for(j in 1:5){
				ij = (i-1)*5+j
				if (is.na(K[k,ij])) {
					price[ij] <- NA
					#ivol[ij] <- NA
				}
				else {
					if (model==1)#Heston
						price[ij] <- hestonAttari(paramsPos[k,], F[k,i], K[k,ij], rd[k,i], tau[i])
					else if (model==2)#Schoebel-Zhu
						price[ij] <- szAttari(paramsPos[k,], F[k,i], K[k,ij], rd[k,i], tau[i])
					else if (model==3)#Bates
						price[ij] <- batesAttari(paramsPos[k,], F[k,i], K[k,ij], rd[k,i], tau[i])
					else if (model==4)#OUOU
						price[ij] <- ououAttari(paramsPos[k,], F[k,i], K[k,ij], rd[k,i], tau[i])
					#if (errType=="IV") ivol[ij] <- impVol(data=price[i,j],t0=0,t=tau[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
					price[ij] <- price[ij]/vega[k,ij]
					errMatrix[i,j] = errMatrix[i,j]  + (price[ij]-mktPrice[k,ij])^2
					#errMatrix[i,j]  + (ivol[ij] - vol[k,ij])^2
				}
			}
		}
		#if (errType=="IV") err[k]=sqrt( sum( (ivol - vol)^2,na.rm=TRUE) )
		err[k]=sqrt(mean((price - mktPrice[k,])^2,na.rm=TRUE))
	}
	list(err,sqrt(errMatrix/m))
}
calibRisk <- function(k,optMSE,optMAE,optMAPE) {
#max(c(abs(optMSE[k]-optMAE[k]-1),abs(optMAE[k]-optMAPE[k]),abs(optMSE[k]-optMAPE[k])))
	apply(cbind(abs(optMSE[,k]-optMAE[,k]),abs(optMAE[,k]-optMAPE[,k]),abs(optMSE[,k]-optMAPE[,k])),1,max)
}
testRMSE <- function(params,m) {
		vol=t(matrix(as.numeric(datatable[m,2:31]/100),nrow = 5,ncol = 6))
		S=as.numeric(datatable[m,32])
		fwd=as.numeric(datatable[m,33:38]/10^4)
		eois=as.numeric(datatable[m,39:44]/100)
		rf=-log(1/(1+eois*tau))/tau
		rd=-log((1+fwd/S)*exp(-rf*tau))/tau
		F=fwd+S
		K10p=F*exp(vol[,1]^2*tau/2)*exp(vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		K25p=F*exp(vol[,2]^2*tau/2)*exp(vol[,2]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		Katmf=F*exp(vol[,3]^2*tau/2)
		K25c=F*exp(vol[,4]^2*tau/2)*exp(-vol[,4]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		K10c=F*exp(vol[,5]^2*tau/2)*exp(-vol[,5]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		K=cbind(K10p,K25p,Katmf,K25c,K10c)
		n=length(tau)
		price <- matrix(0,n,5)
		ivol <- matrix(0,n,5)
		for(i in 1:n){
			for(j in 1:5){
				price[i,j] <- hestongg(S=S,K=K[i,j], kappa=params[4],theta=params[1],sigma=params[2],rho=params[3],lambda=0,nu0=params[5],rd=rd[i],rf=rf[i],t=0,T=tau[i])
				ivol[i,j] <- impVol(data=price[i,j],t0=0,t=tau[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
		}
		sqrt( sum( (ivol - vol)^2 ) )
}
Px <- function(f,K,w,rd,T,type) {
	(B76(f+0.00001,K,w,rd,t=0,T,type)-B76(f-0.00001,K,w,rd,t=0,T,type))/0.00002
}
Pxx <- function(f,K,w,rd,T,type) {
	(B76(f+0.00001,K,w,rd,t=0,T,type)-2*B76(f,K,w,rd,t=0,T,type)+B76(f-0.00001,K,w,rd,t=0,T,type))/0.00002^2
}
Pww <- function(f,K,w,rd,T,type) {
	(B76(f,K,w+0.00001,rd,t=0,T,type)-2*B76(f,K,w,rd,t=0,T,type)+B76(f,K,w-0.00001,rd,t=0,T,type))/0.00002^2
}
gauthier <- function(nu0,theta,kappa,S,K1,K2,vol1,vol2,rd,rf,tau)
{
	rho=rep(NA,6);omega=rho;
	for (j in 1:6) {
		T = tau[j]
		r0 = 1/4 * 1/kappa^3 * (-4*exp(-kappa*T)*kappa*T+2-2*exp(-2*kappa*T))
		r1 = 1/4 * 1/kappa^3 * (4*exp(-kappa*T)*(kappa*T+1)+(2*kappa*T-5)+exp(-2*kappa*T))
		p0 = 1/kappa^2 * (-exp(-kappa*T)*kappa*T+1-exp(-kappa*T))
		p1 = 1/kappa^2 * (exp(-kappa*T)*kappa*T+(kappa*T-2)+2*exp(-kappa*T))
		#^na kappaoncu + czy - ?
		q0 = 1/2 * 1/kappa^3 * (-exp(-kappa*T)*kappa*T*(kappa*T+2)+2-2*exp(-kappa*T))
		q1 = 1/2 * 1/kappa^3 * (2*(kappa*T-3)+exp(-kappa*T)*kappa*T*(kappa*T+4)+6*exp(-kappa*T))
		#w = theta + (nu0 - theta)*(1-exp(kappa*T))/kappa
		w = nu0 * (1-exp(kappa*T))/kappa + theta *(T- (1-exp(kappa*T))/kappa)

		#P1 + B*sigma^2 + C*rho*sigma + D*sigma^2*rho^2 - P1mkt = 0
		#P2 + E*sigma^2 + F*rho*sigma + G*sigma^2*rho^2 - P2mkt = 0
		
		f=log(S)+(rd[j]-rf[j])*(T)
		A = B76(f,K1[j],w,rd[j],t=0,T,-1) - BSMFX(S,K1[j],vol1[j],rd[j],rf[j],t=0,T,-1)
		E = B76(f,K2[j],w,rd[j],t=0,T,1) - BSMFX(S,K2[j],vol2[j],rd[j],rf[j],t=0,T,1)
		P1ww = Pww(f,K1[j],w,rd[j],T,-1)
		P2ww = Pww(f,K2[j],w,rd[j],T,1)
		P1xw = (Px(f,K1[j],w+0.00001,rd[j],T,-1)-Px(f,K1[j],w-0.00001,rd[j],T,-1))/0.00002
		P2xw = (Px(f,K2[j],w+0.000001,rd[j],T,1)-Px(f,K2[j],w-0.00001,rd[j],T,1))/0.000002
		P1xxw = (Pxx(f,K1[j],w+0.00001,rd[j],T,-1)-Pxx(f,K1[j],w-0.00001,rd[j],T,-1))/0.00002
		P2xxw = (Pxx(f,K2[j],w+0.00001,rd[j],T,1)-Pxx(f,K2[j],w-0.00001,rd[j],T,1))/0.00002
		P1xxww = (Pxx(f,K1[j],w+0.00001,rd[j],T,-1)-2*Pxx(f,K1[j],w,rd[j],T,-1)+Pxx(f,K1[j],w-0.00001,rd[j],T,-1))/0.00002^2
		P2xxww = (Pxx(f,K2[j],w+0.00001,rd[j],T,1)-2*Pxx(f,K2[j],w,rd[j],T,1)+Pxx(f,K2[j],w-0.00001,rd[j],T,1))/0.00002^2
		B = (nu0*r0+theta*r1)*P1ww
		C = (nu0*p0+theta*p1)*P1xw #== a_1
		#a_1 = rho*sigma *exp(-kappa*T)/kappa^2*(nu0*(-kappa*T+exp(kappa*T)-1)+theta*(kappa*T+exp(kappa*T)*(kappa*T-2)+2))
		D = (nu0*q0+theta*q1)*P1xxw + 1/2*(nu0*p0+theta*p1)^2*P1xxww #==b_2==a_1^2/2
		#D = (nu0*q0+theta*q1)*P1xxw + 1/2
		
		F = (nu0*r0+theta*r1)*P2ww
		G = (nu0*p0+theta*p1)*P2xw
		H = (nu0*q0+theta*q1)*P2xxw + 1/2*(nu0*p0+theta*p1)^2*P2xxww
		prim = D*E/H-A
		bis = D*F/H-B
		tert = C - D*G/H
		omega[j] = sqrt(tert^2/bis^2/2/D*(-B-C*bis/tert-2*D*prim*bis/tert-sqrt((B+C*bis/tert+2*D*prim*bis/tert)^2-4*D*bis^2/tert^2*(A+C*prim/tert+D*prim^2/tert^2))))
		rho[j] = prim/tert/omega[j] + bis/tert*omega[j]
	}
	cbind(omega,rho)
}
gauthierAprox <- function(nu0,theta,kappa,sigma,rho,S,K,rd,rf,T,type=1) {
	r0 = 1/4 * 1/kappa^3 * (-4*exp(-kappa*T)*kappa*T+2-2*exp(-2*kappa*T))
	#r0 = 1/4 * 1/kappa^3 * exp(-2*kappa*T)*(2*exp(2*kappa*T)-4*exp(kappa*T)*kappa*T-2)
	
	r1 = 1/4 * 1/kappa^3 * (4*exp(-kappa*T)*(kappa*T+1)+(2*kappa*T-5)+exp(-2*kappa*T))
	#r1 = 1/4 * 1/kappa^3 * exp(-2*kappa*T)*(exp(2*kappa*T)*(2*kappa*T-5)+4*exp(kappa*T)*(kappa*T+1)+1)
	
	p0 = 1/kappa^2 * (-exp(-kappa*T)*kappa*T+1-exp(-kappa*T))
	#p0 = 1/kappa^2 * exp(-kappa*T) * (exp(kappa*T)-kappa*T-1)
	
	p1 = 1/kappa^2 * (exp(-kappa*T)*kappa*T+(kappa*T-2)+2*exp(-kappa*T))
	#p1 = 1/kappa^2 * exp(-kappa*T) * (exp(kappa*T)*(kappa*T-2)+kappa*T-2)
	
	#^na koncu + czy - ? odp: +
	
	q0 = 1/2 * 1/kappa^3 * (-exp(-kappa*T)*kappa*T*(kappa*T+2)+2-2*exp(-kappa*T))
	#q0 = 1/2 * 1/kappa^3 * exp(-kappa*T) * (2*exp(kappa*T)-kappa*T*(kappa*T+2)-2)
	
	q1 = 1/2 * 1/kappa^3 * (2*(kappa*T-3)+exp(-kappa*T)*kappa*T*(kappa*T+4)+6*exp(-kappa*T))
	#q1 = 1/2 * 1/kappa^3 * exp(-kappa*T) * (2*exp(kappa*T)*(kappa*T-3)+kappa*T*(kappa*T+4)+6)
	
	#w = theta + (nu0 - theta)*(1-exp(-kappa*T))/kappa
	w = nu0 * (1-exp(-kappa*T))/kappa + theta *(T- (1-exp(-kappa*T))/kappa)

	f=log(S)+(rd-rf)*(T)
	A = B76(f,K,w,rd,t=0,T,type)
	P1ww = Pww(f,K,w,rd,T,type)
	P1xw = (Px(f,K,w+0.00001,rd,T,type)-Px(f,K,w-0.00001,rd,T,type))/0.00002
	P1xxw = (Pxx(f,K,w+0.00001,rd,T,type)-Pxx(f,K,w-0.00001,rd,T,type))/0.00002
	P1xxww = (Pxx(f,K,w+0.00001,rd,T,type)-2*Pxx(f,K,w,rd,T,type)+Pxx(f,K,w-0.00001,rd,T,type))/0.00002^2
	B = (nu0*r0+theta*r1)*P1ww
	C = (nu0*p0+theta*p1)*P1xw #== a_1
	#a_1 = rho*sigma *exp(-kappa*T)/kappa^2*(nu0*(-kappa*T+exp(kappa*T)-1)+theta*(kappa*T+exp(kappa*T)*(kappa*T-2)+2))
	D = (nu0*q0+theta*q1)*P1xxw + 1/2*(nu0*p0+theta*p1)^2*P1xxww 
	
	A + B*sigma^2 + C*rho*sigma + D*sigma^2*rho^2
}
gauthierOmega <- function(omega,nu0,theta,kappa,S,K1,K2,vol1,vol2,rd,rf,T) {
	r0 = 1/4 * 1/kappa^3 * (-4*exp(-kappa*T)*kappa*T+2-2*exp(-2*kappa*T))
	r1 = 1/4 * 1/kappa^3 * (4*exp(-kappa*T)*(kappa*T+1)+(2*kappa*T-5)+exp(-2*kappa*T))
	p0 = 1/kappa^2 * (-exp(-kappa*T)*kappa*T+1-exp(-kappa*T))
	p1 = 1/kappa^2 * (exp(-kappa*T)*kappa*T+(kappa*T-2)+2*exp(-kappa*T))
	#^na koncu + czy - ?
	q0 = 1/2 * 1/kappa^3 * (-exp(-kappa*T)*kappa*T*(kappa*T+2)+2-2*exp(-kappa*T))
	q1 = 1/2 * 1/kappa^3 * (2*(kappa*T-3)+exp(-kappa*T)*kappa*T*(kappa*T+4)+6*exp(-kappa*T))
	#w = theta + (nu0 - theta)*(1-exp(kappa*T))/kappa
	w = nu0 * (1-exp(kappa*T))/kappa + theta *(T- (1-exp(kappa*T))/kappa)
	
	f=log(S)+(rd-rf)*(T)
	A = B76(f,K1,w,rd,t=0,T,-1) - BSMFX(S,K1,vol1,rd,rf,t=0,T,-1)
	E = B76(f,K2,w,rd,t=0,T,1) - BSMFX(S,K2,vol2,rd,rf,t=0,T,1)
	P1ww = Pww(f,K1,w,rd,T,-1)
	P2ww = Pww(f,K2,w,rd,T,1)
	P1xw = (Px(f,K1,w+0.00001,rd,T,-1)-Px(f,K1,w-0.00001,rd,T,-1))/0.00002
	P2xw = (Px(f,K2,w+0.000001,rd,T,1)-Px(f,K2,w-0.00001,rd,T,1))/0.000002
	P1xxw = (Pxx(f,K1,w+0.00001,rd,T,-1)-Pxx(f,K1,w-0.00001,rd,T,-1))/0.00002
	P2xxw = (Pxx(f,K2,w+0.00001,rd,T,1)-Pxx(f,K2,w-0.00001,rd,T,1))/0.00002
	P1xxww = (Pxx(f,K1,w+0.00001,rd,T,-1)-2*Pxx(f,K1,w,rd,T,-1)+Pxx(f,K1,w-0.00001,rd,T,-1))/0.00002^2
	P2xxww = (Pxx(f,K2,w+0.00001,rd,T,1)-2*Pxx(f,K2,w,rd,T,1)+Pxx(f,K2,w-0.00001,rd,T,1))/0.00002^2
	B = (nu0*r0+theta*r1)*P1ww
	C = (nu0*p0+theta*p1)*P1xw #== a_1
	#a_1 = rho*sigma *exp(-kappa*T)/kappa^2*(nu0*(-kappa*T+exp(kappa*T)-1)+theta*(kappa*T+exp(kappa*T)*(kappa*T-2)+2))
	D = (nu0*q0+theta*q1)*P1xxw + 1/2*(nu0*p0+theta*p1)^2*P1xxww #==b_2==a_1^2/2
	
	F = (nu0*r0+theta*r1)*P2ww
	G = (nu0*p0+theta*p1)*P2xw
	H = (nu0*q0+theta*q1)*P2xxw + 1/2*(nu0*p0+theta*p1)^2*P2xxww
	prim = D*E/H-A
	bis = D*F/H-B
	tert = C - D*G/H
	rho = prim/tert/omega + bis/tert*omega
	
	fun1 = A + B*omega^2 + C*rho*omega + D*omega^2*rho^2
	fun2 = E + E*omega^2 + F*rho*omega + G*omega^2*rho^2
	
	abs(fun1)+abs(fun2)
}
hestonTS2d<-function(params,nu0,data,tau){
	nuR = pnorm(params[1])
	theta1 = exp(params[2])
	theta2 = exp(params[3])
	kappa1 = exp(params[4])
	kappa2 = exp(params[5])
	variance=rep(0,5)
	for(i in 1:6){
		variance[i] <- theta1 + (nu0*nuR-theta1) * (1-exp(-kappa1*tau[i]))/(kappa1*tau[i]) + theta2 + (nu0*(1-nuR)-theta2) * (1-exp(-kappa2*tau[i]))/(kappa2*tau[i])
	}
	sqrt( sum( (variance - data)^2 ) )
}
hestonTS2dfit<-function(vol,tau){
	n = dim(vol)[1]
	hestonTSpar=matrix(0,n,5)
	for(i in 1:n){
		params=c(qnorm(0.5),log(vol[i,6]^2)/2,log(vol[i,6]^2)/2,log(2),log(2))
		opt=optim(params,hestonTS2d,nu0=vol[i,1]^2,data=vol[i,]^2,tau=tau);
		hestonTSpar[i,1]=opt$par[1]
		hestonTSpar[i,2]=opt$par[2]
		hestonTSpar[i,3]=opt$par[3]
		hestonTSpar[i,4]=opt$par[4]
		hestonTSpar[i,5]=opt$par[5]
	}
	hestonTSpar[,1] = pnorm(hestonTSpar[,1])
	hestonTSpar[,2:5]=exp(hestonTSpar[,2:5])
	
	nuR=matrix(0,n,6)
	for(i in 1:6){
		nuR[,i] = (hestonTSpar[,2] + (VIX[,1]^2*hestonTSpar[,1]-hestonTSpar[,2]) * (1-exp(-hestonTSpar[,4]*tau[i]))/(hestonTSpar[,4]*tau[i])) / (hestonTSpar[,2] + (VIX[,1]^2*hestonTSpar[,1]-hestonTSpar[,2]) * (1-exp(-hestonTSpar[,4]*tau[i]))/(hestonTSpar[,4]*tau[i]) + hestonTSpar[,3] + (VIX[,1]^2*(1-hestonTSpar[,1])-hestonTSpar[,3]) * (1-exp(-hestonTSpar[,5]*tau[i]))/(hestonTSpar[,5]*tau[i]))
	}
	cbind(nuR,hestonTSpar)
}
E2H2IM <- function(params,omegas,rhos,nuR,nu0,nu0R,theta1,theta2,kappa1,kappa2,data,S,K,rd,rf,t,T,vega,costfun="MSE"){
	if(length(data)<6) data=rbind(data,data)
	sigma1 <- exp(params[1])
	rho1 <- tanh(params[2])
	
	if (mean((nuR*sigma1^2-omegas^2)/(nuR-1)) < 0) {
		NA
	}
	else {
	sigma2 = sqrt((mean((nuR*sigma1^2-omegas^2)/(nuR-1))))
	rho2 = min(1,max(-1,mean((nuR*sigma1*rho1-omegas*rhos)/(nuR-1))/sigma2))
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1=kappa1,theta1=theta1,sigma1=sigma1,rho1=rho1,lambda1=0,nu1=nu0*nu0R,kappa2=kappa2,theta2=theta2,sigma2=sigma2,rho2=rho2,lambda2=0,nu2=nu0*(1-nu0R),rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol") {
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			}
			else {
				price[i,j] <- price[i,j] / vega[i,j]
			}
		}
	}
	if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		if (costfun=="MSE") sum( (price - data)^2 )
		else {if (costfun=="MAE") sum( abs(price - data) )
		else {if (costfun=="MSPE") sum( (price/data-1)^2 )
		else {if (costfun=="MAPE") sum( abs(price/data-1) )
		else {if (costfun=="vol")  sum( (vol - data)^2 ) }}}}
	else
		NA
	}
}
tester2d <- function(datatable,tau,omegas,rhos,nuR,nu0,nu0R,theta1,theta2,kappa1,kappa2){
	a=1;b=6
	lambda=0
	m=dim(datatable)[1]
	err=rep(NA,m)
	for	(k in 1:m) {
		vol<-t(matrix(as.numeric(datatable[k,2:31]/100),nrow = 5,ncol = 6))
		S<-as.numeric(datatable[k,32])
		fwd<-as.numeric(datatable[k,33:38]/10^4)
		eois<-as.numeric(datatable[k,39:44]/100)
		rf<--log(1/(1+eois*tau))/tau
		rd<--log((1+fwd/S)*exp(-rf*tau))/tau
		F<-S*exp((rd-rf)*tau) #F2=S+fwd
		K10p=F*exp(vol[,1]^2*tau/2)*exp(vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		K25p=F*exp(vol[,2]^2*tau/2)*exp(vol[,2]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		Katmf=F*exp(vol[,3]^2*tau/2)
		K25c=F*exp(vol[,4]^2*tau/2)*exp(-vol[,4]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
		K10c=F*exp(vol[,5]^2*tau/2)*exp(-vol[,5]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
		K<-cbind(K10p,K25p,Katmf,K25c,K10c)
		n <- length(tau)
		mktPrice=matrix(0,n,5)
		vega=matrix(0,n,5)
		for(i in 1:n){
			for(j in 1:5){
				d1 <- (log(S/K[i,j]) + (rd[i] - rf[i] + vol[i,j]^2/2)*tau[i])/(vol[i,j]*sqrt(tau[i]))
				d2 <- d1 - vol[i,j]*sqrt(tau[i])
				vega[i,j] <- K[i,j]*exp(-rd[i]*tau[i])*sqrt(tau[i])*dnorm(d2)
				mktPrice[i,j] = 1 * (exp(-rf[i]*tau[i])*S*pnorm(1*d1)-exp(-rd[i]*tau[i])*K[i,j]*pnorm(1*d2))/vega[i,j]
			}
		}
		opt=optim(c(log(mean(omegas[k,])),atanh(mean(rhos[k,]))),E2H2IM,omegas=omegas[k,],rhos=rhos[k,],nuR=nuR[k,],nu0=nu0[k],nu0R=nu0R[k],theta1=theta1[k],theta2=theta2[k],kappa1=kappa1[k],kappa2=kappa2[k],data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun="vol")
		sigma1 = exp(opt$par[1])
		rho1 = tanh(opt$par[2])
		sigma2 = sqrt(mean((nuR[k,]*sigma1^2-omegas[k,]^2)/(nuR[k,]-1)))
		rho2 = min(1,max(-1,mean((nuR[k,]*sigma1*rho1-omegas[k,]*rhos[k,])/(nuR[k,]-1))/sigma2))
		err[k] = E2H2IM(params=c(log(sigma1),atanh(rho1)),omegas=omegas[k,],rhos=rhos[k,],nuR=nuR[k,],nu0=nu0[k],nu0R=nu0R[k],theta1=theta1[k],theta2=theta2[k],kappa1=kappa1[k],kappa2=kappa2[k],data=vol[a:b,],S=S,K=K[a:b,],rd=rd[a:b],rf=rf[a:b],t=0,T=tau[a:b],vega=vega,costfun="vol")
	}
	err
}
E2H7IM <- function(params,omegas,rhos,nu0,data,S,K,rd,rf,t,T,vega,costfun="MSE"){
	if(length(data)<6) data=rbind(data,data)
	theta1 <- exp(params[1])
	sigma1 <- exp(params[2])
	rho1 <- tanh(params[3])
	kappa1 <- exp(params[4])
	nu0R <- pnorm(params[5])
	theta2 <- exp(params[6])
	kappa2 <-exp(params[7])
	
	nuR=matrix(0,n,6)
	for(i in 1:6){
		nuR[,i] = (theta1 + (nu0*nu0R-theta1) * (1-exp(-kappa1*tau[i]))/(kappa1*tau[i]))/(theta1 + (nu0*nu0R-theta1) * (1-exp(-kappa1*tau[i]))/(kappa1*tau[i]) + theta2 + (nu0*(1-nu0R)-theta2) * (1-exp(-kappa2*tau[i]))/(kappa2*tau[i]))
	}

	if (mean((nuR*sigma1^2-omegas^2)/(nuR-1)) < 0) {
		NA
	}
	else {	
	sigma2 = sqrt((mean((nuR*sigma1^2-omegas^2)/(nuR-1))))
	rho2 = min(1,max(-1,mean((nuR*sigma1*rho1-omegas*rhos)/(nuR-1))/sigma2))
	#lambda=0
	n <- length(T)
	price <- matrix(0,n,5)
	vol <- matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			price[i,j] <- heston2(S=S,K=K[i,j],kappa1=kappa1,theta1=theta1,sigma1=sigma1,rho1=rho1,lambda1=0,nu1=nu0*nu0R,kappa2=kappa2,theta2=theta2,sigma2=sigma2,rho2=rho2,lambda2=0,nu2=nu0*(1-nu0R),rd=rd[i],rf=rf[i],t=0,T=T[i])
			if (costfun=="vol")
				vol[i,j] <- impVol(data=price[i,j],t0=0,t=T[i],S=S,K=K[i,j],rd=rd[i],rf=rf[i])
			else
				price[i,j] <- price[i,j] / vega[i,j]
		}
	}
	if (2*kappa1*theta1>sigma1^2 & 2*kappa2*theta2>sigma2^2)
		if (costfun=="MSE") sum( (price - data)^2 )
		else {if (costfun=="MAE") sum( abs(price - data) )
		else {if (costfun=="MSPE") sum( (price/data-1)^2 )
		else {if (costfun=="MAPE") sum( abs(price/data-1) )
		else {if (costfun=="vol") sum( (vol - data)^2 ) }}}}
	else
		NA
	}
}
propozycje<-function(varvar,covar,varint,covarint,errType="vega") {
	omegaSqr = varvar/varint
	omega = sqrt(rowMeans(omegaSqr))
	rho = pmax(rep(-1,n),rowMeans(covar/covarint/omega))
	tester(datatable,pair,tau,hestonTSpar[,1],hestonTSpar[,2],hestonTSpar[,3],rho,omega,errType)
}
kappaDurr <- function(theta,nu0,omega,rho,M) {
	1/(2*(theta-nu0))*(8*M*sqrt(nu0)+omega^2/6/nu0*(2-rho^2/2)+omega*rho*nu0)
}
vec2mat <- function(vec,m) {
	n = length(vec)/m
	mat = matrix(0,n,m)
	for (i in 1:n) {
		for (j in 1:m) {
			mat[i,j] = vec[(i-1)*m+j];
		}
	}
	mat
}

cumax <- function(x){
	n=length(x)
	y=numeric(n)
	y[1]=x[1]
	for(i in 2:n){
		y[i]=x[i]-max(x[1:(i-1)])
	}
	y
}

varTSfitAltC <- function(termstruct,modType,tau,a=1,b=6,kappa=2,err=FALSE,max_it=8000,barrier=1,levelup=1.4,leveldn=0.6,omega=0) {
	#m=dim(termstruct)[2]
	#vec=nlHestonTSv(as.numeric(t(params)),as.numeric(t(termstruct)),tau,err)
	#vec2mat(vec,m)
	n=dim(termstruct)[1]
	if (err) k=1
	else k=0
	out = matrix(0,n,3+k)
	if (modType==1)
		omega1=rep(0,n)
	else 
		omega1=omega
	for(i in 1:n){
		params=c(termstruct[i,1]^(1/modType),termstruct[i,6]^(1/modType),kappa)
		out[i,] =nlVarTS(params,termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
		#out[which(out[,1]==0),1] = exp(-200)
		if (out[i,1]==0) {out[i,1]=exp(-200)}
		if (out[i,2]==0) {out[i,2]=exp(-200)}
		if (out[i,1]==Inf) {out[i,1]=exp(700)}
		if (out[i,2]==Inf) {out[i,2]=exp(700)}
		if (i>1) {
			#if(i>5) last=mean(log(out[(i-5):(i-1),1]))
			#else last=log(out[(i-1),1])
			last=max(log(out[1:(i-1),1]))
			if (log(out[i,1])-last > barrier)
				out[i,] =nlVarTS(params*c(levelup,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			#else if (log(out[i,1])-last < -barrier)
			#	out[i,] =nlVarTS(params*c(leveldn,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			#if(i>5) last=mean(log(out[(i-5):(i-1),2]))
			#else last=log(out[(i-1),2])
			last=max(log(out[1:(i-1),2]))
			if (log(out[i,2])-last > barrier)
				out[i,] =nlVarTS(params*c(1,levelup,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			#else if (log(out[i,2])-last < -barrier)
			#	out[i,] =nlVarTS(params*c(1,leveldn,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])				
		}
		if (out[i,1]==0) {out[i,1]=exp(-200)}
		if (out[i,2]==0) {out[i,2]=exp(-200)}
		#errTS[i] = hestonTS(log(out[i,]),termstruct[i,2:5],tau[a:b])
	}
	out
}

varTSfitAlt2C <- function(termstruct,modType,tau,a=1,b=6,kappa=2,err=FALSE,max_it=8000,barrier=1,levelup=1.5,leveldn=0.5,omega=0) {
	#m=dim(termstruct)[2]
	#vec=nlHestonTSv(as.numeric(t(params)),as.numeric(t(termstruct)),tau,err)
	#vec2mat(vec,m)
	n=dim(termstruct)[1]
	if (err) k=1
	else k=0
	out = matrix(0,n,3+k)
	if (modType==1)
		omega1=rep(0,n)
	else 
		omega1=omega
	for(i in 1:n){
		params=c(termstruct[i,1],termstruct[i,6],kappa)
		out[i,] =nlVarTS(params,termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
		if (i>5) {
		  if ((out[i,1])/mean((out[(i-5):(i-1),1]))> exp(barrier))
				out[i,] =nlVarTS(params*c(levelup,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			else if ((out[i,1])/mean((out[(i-5):(i-1),1]))< exp(-barrier))
				out[i,] =nlVarTS(params*c(leveldn,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			if ((out[i,2])/mean((out[(i-5):(i-1),2])) > exp(barrier))
				out[i,] =nlVarTS(params*c(1,levelup,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
			else if ((out[i,2])/mean((out[(i-5):(i-1),2]))< exp(-barrier))
				out[i,] =nlVarTS(params*c(1,leveldn,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])				
		}
		#errTS[i] = hestonTS(log(out[i,]),termstruct[i,2:5],tau[a:b])
	}
	out
}

varTSfitAlt3C <- function(termstruct,modType,tau,a=1,b=6,kappa=2,err=FALSE,max_it=8000,barrier=1,levelup=1.5,leveldn=0.5,omega=0) {
  #m=dim(termstruct)[2]
  #vec=nlHestonTSv(as.numeric(t(params)),as.numeric(t(termstruct)),tau,err)
  #vec2mat(vec,m)
  n=dim(termstruct)[1]
  if (err) k=1
  else k=0
  out = matrix(0,n,3+k)
  if (modType==1)
    omega1=rep(0,n)
  else 
    omega1=omega
  for(i in 1:n){
    params=c(termstruct[i,1],termstruct[i,6],kappa)
    out[i,] =nlVarTS(params,termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
    if (i>5) {
      last=max(out[1:(i-1),1])
      if (max(out[i,1]/last,(out[i,1])/mean((out[(i-5):(i-1),1]))) > exp(barrier))
        out[i,] =nlVarTS(params*c(levelup,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
      else if (min(out[i,1]/last,(out[i,1])/mean((out[(i-5):(i-1),1]))) < exp(-barrier))
        out[i,] =nlVarTS(params*c(leveldn,1,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
      if (max(out[i,2]/last,(out[i,2])/mean((out[(i-5):(i-1),2]))) > exp(barrier))
        out[i,] =nlVarTS(params*c(1,levelup,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
      else if (min(out[i,2]/last,(out[i,2])/mean((out[(i-5):(i-1),2]))) < exp(-barrier))
        out[i,] =nlVarTS(params*c(1,leveldn,1),termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])				
    }
    #errTS[i] = hestonTS(log(out[i,]),termstruct[i,2:5],tau[a:b])
  }
  out
}

varTSfitC <- function(termstruct,modType,tau,a=1,b=6,kappa=2,err=FALSE,max_it=8000,omega=0) {
	#m=dim(termstruct)[2]
	#vec=nlHestonTSv(as.numeric(t(params)),as.numeric(t(termstruct)),tau,err)
	#vec2mat(vec,m)
	n=dim(termstruct)[1]
	if (err) k=1
	else k=0
	out = matrix(0,n,3+k)
	if (modType==1)
		omega1=rep(0,n)
	else 
		omega1=omega
	for(i in 1:n){
		params=c(termstruct[i,1],termstruct[i,6],kappa)
		out[i,] = nlVarTS(params,termstruct[i,a:b],tau[a:b],modType,err,max_it,omega1[i])
		#errTS[i] = hestonTS(log(out[i,]),termstruct[i,2:5],tau[a:b])
	}
	out
}

lagger<-function(x){
	c(NA,x[-length(x)])
}
get_rates<-function(datatable,pair){
	if (pair=="EURUSD") {
		fwd=datatable[,33:38]/10^4
		ois=datatable[,39:44]/100
		fois=datatable[,45:50]/100
	}
	if (pair=="GBPUSD") {
		fwd=(datatable[,33:38]/10^4)
		ois=(datatable[,39:44]/100)
		fois=(datatable[,45:50]/100)
	}
	if (pair=="USDJPY") {
		fwd=(datatable[,33:38]/10^2)
		fois=(datatable[,39:44]/100)
		ois=(datatable[,45:50]/100)
	}
	if (pair=="USDPLN") {
		#fwd=as.numeric(datatable[,33:38]/10^4)
		fwd=(datatable[,33:38]/10^4)
		fois=(datatable[,39:44]/100)
		#ois=as.numeric(datatable[,45:50]/100)
	}
	spot=as.numeric(datatable[,32])
	rf=-log(1/(1+fois[,1]*tau[1]))/tau[1]
	for (i in 2:6) rf=cbind(rf,-log(1/(1+fois[,i]*tau[i]))/tau[i])
	rd=log((1+fwd[,1]/spot)*exp(rf[,1]*tau[1]))/tau[1]
	for (i in 2:6) rd=cbind(rd,log((1+fwd[,i]/spot)*exp(rf[,i]*tau[i]))/tau[i])
	as.matrix(cbind(rd,rf,fwd))
	#fwdr=datatable[,33:38]/10^4/spot
	#rdiff=log(1+fwdr[,1])/tau[1]
	#for (i in 2:6) rdiff=cbind(rdiff,log(1+fwdr[,i])/tau[i])
}
get_K<-function(F,vol,rd,rf,tau){
	K = matrix(NA,dim(F)[1],30)
	K[,1]=F[,1]*exp(vol[,1:5]^2*tau[,1]/2)*exp(vol[,1:5]*sqrt(tau[,1])*qnorm(0.10/exp(-rf[,1]*tau[,1])))
	K[,2] = F[,2]*exp(vol[,6:10]^2*tau[,2]/2)*exp(vol[,2]*sqrt(tau[,2])*qnorm(0.25/exp(-rf[,2]*tau[,2])))
	K[,3] = F[,3]*exp(vol[,11:15]^2*tau[,3]/2)
	K[,4] = F[,4]*exp(vol[,16:20]^2*tau[,4]/2)*exp(-vol[,4]*sqrt(tau[,4])*qnorm(0.25/exp(-rf[,4]*tau[,4])))
	K[,5] = F[,5]*exp(vol[,21:2]^2*tau[,5]/2)*exp(-vol[,5]*sqrt(tau[,5])*qnorm(0.10/exp(-rf[,5]*tau[,5])))
	
	for(i in 2:6) K10p=F*exp(vol[,(i-1)*5+(1:5)]^2*tau/2)*exp(vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
	
	cbind(K10p,K25p,Katmf,K25c,K10c)
}
get_options<-function(vol,F,rd,rf,tau){
	m = length(S)
	n = length(tau)
	out = matrix(NA,m,90)
	for (k in 1:m) {
		volmat=t(matrix(as.numeric(vol[k,]),nrow = 5,ncol = 6))
		K10p=F[k,]*exp(volmat[,1]^2*tau/2)*exp(volmat[,1]*sqrt(tau)*qnorm(0.10/exp(-as.numeric(rf[k,])*tau)))
		K25p=F[k,]*exp(volmat[,2]^2*tau/2)*exp(volmat[,2]*sqrt(tau)*qnorm(0.25/exp(-as.numeric(rf[k,])*tau)))
		Katmf=F[k,]*exp(volmat[,3]^2*tau/2)
		K25c=F[k,]*exp(volmat[,4]^2*tau/2)*exp(-volmat[,4]*sqrt(tau)*qnorm(0.25/exp(-as.numeric(rf[k,])*tau)))
		K10c=F[k,]*exp(volmat[,5]^2*tau/2)*exp(-volmat[,5]*sqrt(tau)*qnorm(0.10/exp(-as.numeric(rf[k,])*tau)))
		K=t(rbind(K10p,K25p,Katmf,K25c,K10c))
		mktPrice=matrix(0,n,5)
		vega=matrix(0,n,5)
		for(i in 1:n){
			for(j in 1:5){
				d1 <- (log(F[k,i]/K[i,j]) + (volmat[i,j]^2/2)*tau[i])/(volmat[i,j]*sqrt(tau[i]))
				d2 <- d1 - volmat[i,j]*sqrt(tau[i])
				vega[i,j] <- K[i,j]*exp(-rd[k,i]*tau[i])*sqrt(tau[i])*dnorm(d2)
				typeo = 1
				mktPrice[i,j] = (typeo * exp(-rd[k,i]*tau[i])*(F[k,i]*pnorm(typeo*d1)-K[i,j]*pnorm(typeo*d2)))/vega[i,j]
			}
		}
		out[k,] = c(as.numeric(t(K)),as.numeric(t(mktPrice)),as.numeric(t(vega)))
	}
	out
}

opterC<-function(params,a,b,i,pair,opType=5,errType=1,modType=1,err=TRUE){
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
		ois=as.numeric(datatable[i,45:50]/100)
	}
	#rd=-log(1/(1+ois*tau))/tau
	#rf=-log((1+fwd/S)*exp(-rd*tau))/tau
	rf=-log(1/(1+fois*tau))/tau
	rd=log((1+fwd/S)*exp(rf*tau))/tau
	F=fwd+S
	K10p=F*exp(vol[,1]^2*tau/2)*exp(vol[,1]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
	K25p=F*exp(vol[,2]^2*tau/2)*exp(vol[,2]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
	Katmf=F*exp(vol[,3]^2*tau/2)
	K25c=F*exp(vol[,4]^2*tau/2)*exp(-vol[,4]*sqrt(tau)*qnorm(0.25/exp(-rf*tau)))
	K10c=F*exp(vol[,5]^2*tau/2)*exp(-vol[,5]*sqrt(tau)*qnorm(0.10/exp(-rf*tau)))
	K=cbind(K10p,K25p,Katmf,K25c,K10c)
	#K <- impStrike(data=vol,t0=0,t=tau,S=S,rd=rd,rf=rf)
	n=length(tau)
	mktPrice=matrix(0,n,5)
	vega=matrix(0,n,5)
	for(i in 1:n){
		for(j in 1:5){
			d1 <- (log(S/K[i,j]) + (rd[i] - rf[i] + vol[i,j]^2/2)*tau[i])/(vol[i,j]*sqrt(tau[i]))
			d2 <- d1 - vol[i,j]*sqrt(tau[i])
			vega[i,j] <- K[i,j]*exp(-rd[i]*tau[i])*sqrt(tau[i])*dnorm(d2)
			typeo = 1
			mktPrice[i,j] = (typeo * (exp(-rf[i]*tau[i])*S*pnorm(typeo*d1)-exp(-rd[i]*tau[i])*K[i,j]*pnorm(typeo*d2)))/vega[i,j]
		}
	}
	nlOpter(params,as.numeric(t(mktPrice)),as.numeric(t(vega)),S,as.numeric(t(K)),rd,rf,tau,opType,errType,modType,err)
}

opterCvParted<-function(i,params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=0,logfile=""){
#,n0,n1) {
	n = dim(params)[1]
	if (sum(as.numeric(is.na(params[,2])))==0) start = 1
	else start = last(which(is.na(params[,2])))+1

	if(i<5) {
	  if (i==1) {
	    start=1; n = 333
	  }
	  else if(i==2){
	    start=334; n = 667
	  }
	  else if(i==3) {
	    start=668; n = 1000
	  }
	  else if(i==4) {
	    start=1001; n = 1333
	  }
	}
	calib = matrix(NA,n+1-start,dim(params)[2]+1)
	for (i in start:n)
	#for (i in n0:n1)
		#if (i==start | mod(i,6)==0) 
		calib[1+i-start,] = nlOpter(params[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust,logfile)
		#else calib[i,] = calib[i-1,]
	calib		
}
opterCvv<-function(vec,params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=0,logfile=""){
  #,n0,n1) {
  n = dim(params)[1]
  if (sum(as.numeric(is.na(params[,2])))==0) start = 1
  else start = last(which(is.na(params[,2])))+1
  calib = matrix(NA,n,dim(params)[2]+1)
  for (i in 1:length(vec))
    #for (i in n0:n1)
    #if (i==start | mod(i,6)==0) 
    calib[vec[i],] = nlOpter(params[vec[i],],mktPrice[vec[i],],vega[vec[i],],F[vec[i],],K[vec[i],],rd[vec[i],],tau,opType,errType,modType,err,max_it,fellerMust,logfile)
  #else calib[i,] = calib[i-1,]
  calib		
}

opterCv<-function(params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=0,logfile=""){
  #,n0,n1) {
  n = dim(params)[1]
  if (sum(as.numeric(is.na(params[,2])))==0) start = 1
  else start = last(which(is.na(params[,2])))+1
  calib = matrix(NA,n,dim(params)[2]+1)
  for (i in start:n)
    #for (i in n0:n1)
    #if (i==start | mod(i,6)==0) 
    calib[i,] = nlOpter(params[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust,logfile)
  #else calib[i,] = calib[i-1,]
  calib		
}

opterCvTrack<-function(params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,step=50) {
	n = dim(params)[1]
	iter = max_it/step
	errors = matrix(0,n,iters)
	calib = cbind(params,rep(0,n))
	for(i in 1:iter) {
		calib = opterCv(calib[-11],mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=step)
		errors[,i] = calib[,11]
	}
	list(calib,errors)
}

opterCv26<-function(params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600) {
	n = dim(params)[1]
	calib = matrix(NA,n,12)
	if (sum(as.numeric(is.na(params[,2])))==0) start = 1
	else start = last(which(is.na(params[,2])))+1
	for (i in start:n) {
		opt1 = nlOpter(params[i,],mktPrice[i,1:5],vega[i,1:5],F[i,1],K[i,1:5],rd[i,1],tau[1],opType,errType,modType,err,max_it)
		opt2 = nlOpter(params[i,],mktPrice[i,6:10],vega[i,6:10],F[i,2],K[i,6:10],rd[i,2],tau[2],opType,errType,modType,err,max_it)
		opt3 = nlOpter(params[i,],mktPrice[i,11:15],vega[i,11:15],F[i,3],K[i,11:15],rd[i,3],tau[3],opType,errType,modType,err,max_it)
		opt4 = nlOpter(params[i,],mktPrice[i,16:20],vega[i,16:20],F[i,4],K[i,16:20],rd[i,4],tau[4],opType,errType,modType,err,max_it)
		opt5 = nlOpter(params[i,],mktPrice[i,21:25],vega[i,21:25],F[i,5],K[i,21:25],rd[i,5],tau[5],opType,errType,modType,err,max_it)
		opt6 = nlOpter(params[i,],mktPrice[i,26:30],vega[i,26:30],F[i,6],K[i,26:30],rd[i,6],tau[6],opType,errType,modType,err,max_it)
		calib[i,] = c(opt1[2],opt2[2],opt3[2],opt4[2],opt5[2],opt6[2],opt1[3],opt2[3],opt3[3],opt4[3],opt5[3],opt6[3])
	}
	calib		
}

refill_table<-function(table){
	n = dim(table)[1]
	m = dim(table)[2]
	out = table
	for (i in 1:n) {
		for (j in 2:m) {
			if(is.na(out[i,j])) out[i,j] = out[i,j-1]
		}
	}
	out		
}
addatend<-function(table,x) {
	n = dim(table)[1]
	m = dim(table)[2]
	out = table
	for (i in 1:n) {
		for (j in 2:m) {
			if(is.na(out[i,j])) out[i,j] = x[i]
		}
	}
	out		
}

calibClean<-function(calibOld,params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=0,barrier=1,levelup=1.5,leveldn=0.5) {
	calib = calibOld
	n = dim(params)[1]
	if (sum(as.numeric(is.na(params[,2])))==0) start = 1
	else start = last(which(is.na(params[,2])))+1
	for (i in start:n) {
		if (i>start) {
			#if(i>5) last=mean(log(calib[(i-5):(i-1),2]))
			#else last=log(calib[(i-1),2])
			last=max(log(calib[1:(i-1),2]))
			if (log(calib[i,2])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(1,levelup,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,2])-mean(log(calib[(i-5):(i-1),2]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(1,leveldn,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)		
			#if(i>5) last=mean(log(calib[(i-5):(i-1),1]))
			#else last=log(calib[(i-1),1])			
			last=max(log(calib[1:(i-1),1]))
			if (log(calib[i,1])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(levelup,1,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,1])-mean(log(calib[(i-5):(i-1),1]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(leveldn,1,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#if(i>5) last=mean(log(calib[(i-5):(i-1),5]))
			#else last=log(calib[(i-1),5])			
			last=max(log(calib[1:(i-1),5]))
			if (log(calib[i,5])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(1,1,1,1,levelup),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,5])-mean(log(calib[(i-5):(i-1),5]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(1,1,1,1,leveldn),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
		}
	}
	calib		
}

calibCleanF<-function(calibOld,params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=0,barrier=1,levelup=1.5,leveldn=0.5,j) {
	calib = calibOld
	n = dim(params)[1]
	if (sum(as.numeric(is.na(params[,2])))==0) start = 1
	else start = last(which(is.na(params[,2])))+1
	for (i in start:j) {
		if (i>start) {
			#if(i>5) last=mean(log(calib[(i-5):(i-1),2]))
			#else last=log(calib[(i-1),2])
			last=max(log(calib[1:(i-1),2]))
			if (log(calib[i,2])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(1,levelup,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,2])-mean(log(calib[(i-5):(i-1),2]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(1,leveldn,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)		
			#if(i>5) last=mean(log(calib[(i-5):(i-1),1]))
			#else last=log(calib[(i-1),1])			
			last=max(log(calib[1:(i-1),1]))
			if (log(calib[i,1])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(levelup,1,1,leveldn,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,1])-mean(log(calib[(i-5):(i-1),1]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(leveldn,1,1,1,1),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#if(i>5) last=mean(log(calib[(i-5):(i-1),5]))
			#else last=log(calib[(i-1),5])			
			last=max(log(calib[1:(i-1),5]))
			if (log(calib[i,5])-last > barrier)
				calib[i,] = nlOpter(params[i,]*c(1,1,1,1,levelup),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
			#else if (log(calib[i,5])-mean(log(calib[(i-5):(i-1),5]))< -barrier)
			#	calib[i,] = nlOpter(params[i,]*c(1,1,1,1,leveldn),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust)
		}
	}
	calib		
}

calibNA<-function(calibOld,params,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,correct=c(2,1/2,1,2,1)) {
	calib = calibOld
	n = dim(params)[1]
	if (sum(as.numeric(is.na(params[,2])))==0) start = 1
	else start = last(which(is.na(params[,2])))+1
	if (dim(params)[2]>5)
		for (i in start:n) {
				if(calibOld[i,11]>4) calib[i,] = nlOpter(as.numeric(calibOld[i,-11])*c(correct,correct),mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust=1)
		}
	else
		for (i in start:n) {
			if(calibOld[i,6]>4) calib[i,] = nlOpter(as.numeric(calibOld[i,-6])*correct,mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,opType,errType,modType,err,max_it,fellerMust=1)
		}
	calib		
}

factor2from1<-function(params,vec) {
	newparams = sweep(params,MARGIN=2,vec,`*`)
	as.matrix(cbind(newparams,newparams))
}

logit<-function(x){1/(1+exp(x))}
alogit<-function(x){log(1/x-1)}
sqrtm<-function(x){sign(x)*sqrt(abs(x))}

err121 <- function(params,omegas,omrho,a) {
	weights = logit(params)
	#omegas1=exp(params[7])^2
	#omrho1=tanh(params[8])*exp(params[7])
	omegas1=max(omegas)*(1+a)
	omegas2=min(omegas)*(1-a)
	omrho1 = min(omrho)*(1-sign(min(omrho))*a)
	omrho2 = max(omrho)*(1+sign(max(omrho))*a)
	#rho=omrho/sqrt(omegas)
	#rho1 = -max(abs(rho))
	#rho2 = min(abs(rho))
	omegasFit = weights*omegas1 + (1-weights)*omegas2
	omrhoFit = weights*omrho1 + (1-weights)*omrho2
	#omrhoFit = weights*rho1*sqrt(omegas1) + (1-weights)*rho2*sqrt(omegas2)
	#sum(abs(median(omegas)-omegasFit)^2)+sum(abs(median(omrho)-omrhoFit)^2)
	sum(abs(omegasFit/median(omegas)-1))+sum(abs(omrhoFit/median(omrho)-1))
	#w = (omegas-omegas2)/(omegas1-omegas2) = (omrho-omrho2)/(omrho1-omrho2)
}

twelve2one <- function(params,omegas,omrho,a=0.1,err=FALSE) {
	params = alogit(params)
	#opt = optim(params,err121,omegas=omegas,omrho=omrho,a=a)
	opt = optimize(err121,c(-5,5),omegas=omegas,omrho=omrho,a=a)
	#if (err) c(logit(opt$par),opt$value)
	#else logit(opt$par)
	if (err) c(logit(opt$minimum),opt$objective)
	else logit(opt$minimum)
}

err128 <- function(params,omegas,omrho,a) {
	weights = logit(params)
	#omegas1=exp(params[7])^2
	#omrho1=tanh(params[8])*exp(params[7])
	omegas1=max(omegas)*(1+a)
	omegas2=min(omegas)*(1-a)
	omrho1 = min(omrho)*(1-sign(min(omrho))*a)
	omrho2 = max(omrho)*(1+sign(max(omrho))*a)
	#rho=omrho/sqrt(omegas)
	#rho1 = -max(abs(rho))
	#rho2 = min(abs(rho))
	omegasFit = weights*omegas1 + (1-weights)*omegas2
	omrhoFit = weights*omrho1 + (1-weights)*omrho2
	#omrhoFit = weights*rho1*sqrt(omegas1) + (1-weights)*rho2*sqrt(omegas2)
	sum(abs(omegas-omegasFit))+sum(abs(omrho-omrhoFit))
	#w = (omegas-omegas2)/(omegas1-omegas2) = (omrho-omrho2)/(omrho1-omrho2)
}

twelve2six <- function(params,omegas,omrho,a=0.1,err=FALSE) {
	params = alogit(params)
	opt = optim(params,err128,omegas=omegas,omrho=omrho,a=a)
	if (err) c(logit(opt$par),opt$value)
	else logit(opt$par)
}

weights1 <- function(omegas,omrho,a=0){
	n=dim(omegas)[1]
	omega1=sqrt(apply(omegas,1,max)*(1+a))
	omega2=sqrt(apply(omegas,1,min)*(1-a))
	#omrho1 = 0.99*omega1
	#omrho2 = -0.99*omega2
	omrho1 = apply(omrho,1,min)*(1-sign(apply(omrho,1,min))*a)
	omrho2 = apply(omrho,1,max)*(1+sign(apply(omrho,1,max))*a)
	weight = numeric(n)
	for (i in 1:n) {
		y = c(omegaLmed[i]^2 - omega2[i]^2,rhoLmed[i]*omegaLmed[i] - omrho2[i])
		x = c(omega1[i]^2-omega2[i]^2,omrho1[i]-omrho2[i])
		lm1 = lm(y~x-1)
		weight[i] = lm1$coeff
		#w = sum(y*x)/sum(x^2)
		if(weight[i]>1 | weight[i]<0)
			weight[i] = logit(optimize(err121,c(-5,5),omegas=omegas[i,],omrho=omrho[i,],a=a)$minimum)
	}
	weight
}

weights2 <- function(omegas,omrho,a=0){
	n=dim(omegas)[1]
	weights = matrix(0,n,6)
	for (i in 1:n) {
		weights[i,] = twelve2six(rep(0.5,6),omegas[i,],omrho[i,],a)
		#weights[i] = twelve2one(0.5,omegas[i,],omrho[i,],a)
	}
	weights
}


weights3 <- function(omegas,omrho,a=0){
	n=dim(omegas)[1]
	omega1=sqrt(apply(omegas,1,max)*(1+a))
	omega2=sqrt(apply(omegas,1,min)*(1-a))
	w = (omegaLmed^2 - omega2^2)/(omega1^2-omega2^2)
}

signabsmax <- function(x) {
	if (abs(max(x)) > abs(min(x)))
		sign(max(x))
	else
		sign(min(x))
}

absmax <- function(x) {
	if (abs(max(x)) > abs(min(x)))
		max(x)
	else
		min(x)
}

absmin <- function(x) {
	if (abs(min(x)) < abs(max(x)))
		min(x)
	else
		max(x)
}

params2factor <- function(omegas,omrho,weights,termstruct,mod=1,a=0,omega=0){
	omega1=sqrt(apply(omegas,1,max)*(1+a))
	omega2=sqrt(apply(omegas,1,min)*(1-a))
	#omrho1 = 0.99*omega1
	#omrho2 = -0.99*omega2
	omrho1 = apply(omrho,1,min)*(1-sign(apply(omrho,1,min))*a)
	omrho2 = apply(omrho,1,max)*(1+sign(apply(omrho,1,max))*a)
	rho1=apply(omrho/sqrt(omegas),1,absmax)
	rho2=apply(omrho/sqrt(omegas),1,absmin)
	varts1=weights*termstruct
	varts2=(1-weights)*termstruct
	if (mod==1 | mod==2) {
		curve1 = varTSfitAltC(varts1,1,tau,2,6,2,TRUE,barrier=1,levelup=1.4,leveldn=0.6)
		curve2 = varTSfitAltC(varts2,1,tau,2,6,2,TRUE,barrier=1,levelup=1.4,leveldn=0.6)
	}
	if (mod==2) {
		omegaM = matrix(rep(omega,6),n,6)
		kappa1 = matrix(rep(curve1[,3],6),n,6); kt = kappa1*tau1; ekt = exp(-kt)
		nu0 = curve1[,1]; theta = curve1[,2]
		timefactor = 2*(nu0-theta)*(1-2*kt*ekt-ekt^2) + theta*(4*ekt-3+2*kt-ekt^2)
		varofvar = omegaM^2/2/tau1^2/kappa1^3*timefactor
		convadj1 = 1-1/8*varofvar/sqrt(varts1)^4
		kappa1 = matrix(rep(curve2[,3],6),n,6); kt = kappa1*tau1; ekt = exp(-kt)
		nu0 = curve2[,1]; theta = curve2[,2]
		timefactor = 2*(nu0-theta)*(1-2*kt*ekt-ekt^2) + theta*(4*ekt-3+2*kt-ekt^2)
		varofvar = omegaM^2/2/tau1^2/kappa1^3*timefactor
		convadj2 = 1-1/8*varofvar/sqrt(varts2)^4
		#convadj = omega		
		vol1 = sqrt(varts1)*convadj1
		vol2 = sqrt(varts2)*convadj2
		curve1 = varTSfitAltC(vol1,1,tau,2,6,0.95,TRUE,barrier=1,levelup=1.4,leveldn=0.6)
		curve2 = varTSfitAltC(vol2,1,tau,2,6,0.95,TRUE,barrier=1,levelup=1.4,leveldn=0.6)	
	}	
	if (mod==3) {
		wo=(omegaLmed^2-omega2^2)/(omega1^2-omega2^2)
		wro=(rhoLmed*omegaLmed-omega2*rho2)/(omega1*rho1-omega2*rho2)
		b=wo/wro
		omega1b = omega1*sqrt(b)
		omega2b = omega2*sqrt((1-b*wro)/(1-wro))
		curve1 = varTSfitAltC(varts1,2,tau,2,6,1.1,TRUE,barrier=1,levelup=1.4,leveldn=0.6,omega=omega1/2/sqrt(2))
		curve2 = varTSfitAltC(varts2,2,tau,2,6,1.1,TRUE,barrier=1,levelup=1.4,leveldn=0.6,omega=omega2/2/sqrt(2))
	}
	cbind(curve1[,2],omega1,omrho1/omega1,curve1[,3],curve1[,1],curve2[,2],omega2,omrho2/omega2,curve2[,3],curve2[,1])
}

calibOutliers<-function(calib,barrier=1) {
	nu0=0;theta=0;omega=0
	nu0d=0;thetad=0;omegad=0
	n = dim(calib)[1]
	if (sum(as.numeric(is.na(calib[,2])))==0) start = 1
	else start = last(which(is.na(calib[,2])))+1
	for (i in start:n) {
		if (i>4+start) {
			if (log(calib[i,2])-mean(log(calib[(i-5):(i-1),2]))> barrier)
				omega=omega+1
			else if (log(calib[i,2])-mean(log(calib[(i-5):(i-1),2]))< -barrier)
				omegad=omegad+1
			if (log(calib[i,1])-mean(log(calib[(i-5):(i-1),1]))> barrier)
				theta=theta+1
			else if (log(calib[i,1])-mean(log(calib[(i-5):(i-1),1]))< -barrier)
				thetad=thetad+1
			if (log(calib[i,5])-mean(log(calib[(i-5):(i-1),5]))> barrier)
				nu0=nu0+1
			else if (log(calib[i,5])-mean(log(calib[(i-5):(i-1),5]))< -barrier)
				nu0d=nu0d+1
		}		
	}
	matrix(c(nu0,nu0d,theta,thetad,omega,omegad),2,3)
}

calibOutliers2<-function(calib,barrier=1) {
	nu0=0;theta=0;omega=0
	n = dim(calib)[1]
	if (sum(as.numeric(is.na(calib[,2])))==0) start = 1
	else start = last(which(is.na(calib[,2])))+1
	for (i in start:n) {
		if (i>start) {
			if (log(calib[i,2])-max(log(calib[1:(i-1),2]))> barrier)
				omega=omega+1
			if (log(calib[i,1])-max(log(calib[1:(i-1),1]))> barrier)
				theta=theta+1
			if (log(calib[i,5])-max(log(calib[1:(i-1),5]))> barrier)
				nu0=nu0+1
		}
	}
	c(nu0,theta,omega)
}

hestonggVol <- function(S,K,rd,rf,t,T,u,I,dtype=1,type=1){
	# type = 1 => Call
	# type = -1 => Put
	i <- complex(real=0,imaginary=1)
	k <- log(K)
	P <- 1/2 + 1/pi * trapz(u,Re(exp(-i*u*k)*I))
	p <- (1-type)/2 + type*P
	d=norm_s_inv(p)
	#d1 = (m + rdiff + v^2/2)/v
	#d2 = (m + rdiff - v^2/2)/v
	#0 = v^2/2 - d1v + m + rdifft
	#0 =-v^2/2 - d2v + m + rdifft
	#v = d1+-sqrt(d^2-2*(m+rdifft))
	#v =-[d2+-sqrt(d^2+2*(m+rdifft))]
	mon = 1#type*sign(log(K/S))
	dtype*(d+mon*sqrt(d^2-dtype*2*(log(S/K)+(rd-rf)*(T-t))))/sqrt(T-t)
}

norm_s_inv <- function(p) {
#pade approximation
	#pi=3.14159
	(sqrt(2*pi)*(p-1/2)-157/231*sqrt(2*pi^3)*(p-1/2)^3)/(1-78/77*pi*(p-1/2)^2+241/2310*pi^2*(p-1/2)^4)
}

impvol_approx <- function(S,K,rd,rf,tau,C,P){
	#pi=3.14159
	beta = sqrt(8/pi)
	d = K*exp(-rd*tau)/S/exp(-rf*tau)
	b = 2/beta * (C+P)/(S*exp(-rf*tau)+K*exp(-rd*tau))
	m = ((1-d)/(1+d))^2
	(b/(1-m/4) + sqrt((b/(1-m/4))^2 - 2*m/(1-m/4)))/sqrt(tau)
}

price_Hardy <- function(S,K,rd,rf,tau,typeo=1){
	#pi=3.14159
	F = S*exp((rd-rf)*tau)
	intrinsic = max(0,typeo*(F-K))
	atmPrice = iv*sqrt(tau)/sqrt(2*pi)
	#intrinsic + atmPrice * (1-0.41*abs(d1))*exp(-abs(d1))
}

opter2phase<-function(params,mktPrice,vega,F,K,rd,tau,errType=1,modType=3,err=TRUE,max_it1=1000,max_it2=500,max_phase=8) {
	calib=c(params,NA)
	i=0
	while (i<max_phase) {
		calib = nlOpter(calib[-11],mktPrice,vega,F,K,rd,tau,8,errType,modType,TRUE,max_it1)
		calib = nlOpter(calib[-11],mktPrice,vega,F,K,rd,tau,2,errType,modType,TRUE,max_it2)
		i=i+1;
		#ratio = calib[11]/err-1
	}
	calib		
}

write.textable<-function(table,filename) {
	write.table(print(table,floating=FALSE,sanitize.colnames.function = identity,sanitize.rownames.function = identity),filename,quote = FALSE)
	con <- file(filename, 'r')
	input = readLines(con)
	close(con)
	file.remove(filename)
	con <- file(filename, 'w') 
	k = length(input)
	#for (i in seq_along(input)) {
	#	input[i] <- gsub('"', '', input[i]) 
	#}
	input = c(input[-c(1,2)])
	writeLines(input, con=con)
	close(con)
	#file.remove(paste0(filename,"_"))
}
heston2d.pricerC <- function(params,m,F,K,rd,vega,mktPrice){
	#n=length(K)
	price=matrix(NA,6,5)
	spots=vec2mat(heston2d_pather(as.numeric(params),m),2*m)
	for (i in 1:6){#tenor loop
		S=exp(spots[i,])*F[i]
		for (j in 1:5){#strike loop
			ij=(i-1)*5+j
			all=mean(pmax(S-K[ij],0))
			#orgprice=batesAttari(as.numeric(params), F[i], K[ij], rd[i], tau[i])
			price[i,j]=(all*exp(-rd[i]*tau[i]))
		}
	}
	price
}
ouou.pricerC <- function(params,m,F,K,rd,vega,mktPrice){
	#n=length(K)
	price=matrix(NA,6,5)
	spots=vec2mat(ouou_pather(as.numeric(params),m),2*m)
	for (i in 1:6){#tenor loop
		S=exp(spots[i,])*F[i]
		for (j in 1:5){#strike loop
			ij=(i-1)*5+j
			all=mean(pmax(S-K[ij],0))
			#orgprice=batesAttari(as.numeric(params), F[i], K[ij], rd[i], tau[i])
			price[i,j]=(all*exp(-rd[i]*tau[i]))
		}
	}
	price
}
ouou2.pricerC <- function(params,m,F,K,rd,vega,mktPrice){
	#n=length(K)
	price=matrix(NA,6,5)
	spots=vec2mat(ouou2_pather(as.numeric(params),m),2*m)
	for (i in 1:6){#tenor loop
		S=exp(spots[i,])*F[i]
		for (j in 1:5){#strike loop
			ij=(i-1)*5+j
			all=mean(pmax(S-K[ij],0))
			#orgprice=batesAttari(as.numeric(params), F[i], K[ij], rd[i], tau[i])
			price[i,j]=(all*exp(-rd[i]*tau[i]))
		}
	}
	price
}

mod<-function(a,b) {
	a%%b
}

test_methods <- function(S,K,vol,rd,rf,pricepervega,vega,out=1){
  params=smileRegVec(S,K,vol,tau)
  SS = as.numeric(params[,1])
  CC = as.numeric(params[,2]*2)
  dK=cbind(K[,2]-K[,1],(K[,3]-K[,1])/2,(K[,4]-K[,2])/2,(K[,5]-K[,3])/2,K[,5]-K[,4])
  GKp=BSMFX(S,K[,1:3],vol[,1:3],rd,rf,0,tau,type=-1)
  GKc=BSMFX(S,K[,3:5],vol[,3:5],rd,rf,0,tau,type=1)
  GK=cbind(GKp[,-3],(GKp[,3]+GKc[,1])/2,GKc[,-1])
  
  F = S*exp((rd-rf)*tau)
  VIX=sqrt(2/tau*exp((rd-rf)*tau)*rowSums(dK*GK/K^2)-1/tau*(F/K[,3]-1)^2)
  P2=2*exp((rd)*tau)*rowSums((1-log(K/F))*dK*GK/K^2)
  P3=3*exp((rd)*tau)*rowSums((2*log(K/F)-log(K/F)^2)*dK*GK/K^2)
  P4=4*exp((rd)*tau)*rowSums((3*log(K/F)^2-log(K/F)^3)*dK*GK/K^2) 
  P1 = exp((rd-rf)*tau) -1 -1/2*P2-1/6*P3-1/24*P4
  
  hestonTSpar3m=varTSfitAltC(as.matrix(t(VIX^2)),1,tau,2,6,2,TRUE,barrier=0.4,levelup=2,leveldn=0.6)
  
  RC2 = (P2-P1^2)
  RC3 = (P3-3*P1*P2+2*P1^3)
  RC4 = (P4-4*P1*P3+6*P1^2*P2-3*P1^4)
  A = RC2
  B = sqrt(1+A)
  exye = (2*B^8-4*B^7+2*B^6+ (2*B^3-2*B^5)*RC2+(2*B^2-B^3)*RC3+1/2*RC4)/B^6
  ey2e = (4*B^8-8*B^7+4*B^6+ (4*B^6-8*B^5+4*B^3)*RC2+(4*B^2-4*B^3)*RC3+RC4)/B^6
  ey2intL = RC2*tau^2/3
  exyintL = RC2*tau/2
  omegaLmed = sqrt(median(ey2e/ey2intL))
  rhoLmed = median(exye/exyintL)/omegaLmed
  
  omegaDA = sqrt(smile2om2(VIX[1],SS[1],CC[1]))
  rhoDA = smile2rho(VIX[1],SS[1],CC[1])

  paramsDurr = c(hestonTSpar3m[2],omegaDA,rhoDA,hestonTSpar3m[3],hestonTSpar3m[1])
  paramsLmed = c(hestonTSpar3m[2],omegaLmed,rhoLmed,hestonTSpar3m[3],hestonTSpar3m[1])
  
  durr=tester1d(paramsDurr,F,as.numeric(t(K)),rep(rd,6),as.numeric(t(pricepervega)),as.numeric(t(vega)),tau,model=1)
  lab=tester1d(paramsLmed,F,as.numeric(t(K)),rep(rd,6),as.numeric(t(pricepervega)),as.numeric(t(vega)),tau,model=1)
  if (out==1) c(paramsDurr,paramsLmed)
  else c(durr,lab)
}
