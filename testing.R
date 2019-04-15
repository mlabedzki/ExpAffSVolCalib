#############################
## Intro                   ##
#############################
options("scipen"=100, "digits"=7)#to force switch of scientific notation
library(Rcpp); Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("nl_expaffsvol.cpp")
source("aff_pricing.R")
library(tseries);
library(xtable)#latex export

tau=c(1/12,2/12,3/12,6/12,1,2)

datatable=read.table("EURUSDall.csv", header = TRUE, sep="\t");pair="EURUSD"
datatable$Date = as.Date(datatable$Date)
n=dim(datatable)[1]
tau1 = sweep(matrix(1,n,6),MARGIN=2,tau,`*`)
P1=matrix(NA,n,6);VIX=P1;P2=P1;P3=P1;P4=P1;P5=P1;P6=P1;SS=P1;CC=P1;mu=P1;Q1=P1;Q2=P1;Q3=P1;Q4=P1;l1=P1;l2=P1;skew=P1;kurt=P1;level=P1;S_Leif=P1[,1];C_Leif=P1[,1]
data.load(datatable,pair)

S = as.numeric(datatable[,32])
vol = datatable[,2:31]/100
rates = get_rates(datatable,pair)
rd = rates[,1:6]
rf = rates[,7:12]
fwd = rates[,13:18]
F = S+fwd
options = get_options(vol,F,rd,rf,tau)
K = options[,1:30]
mktPrice = options[,31:60]
vega = options[,61:90]
remove(rates);remove(options)

#############################
## Variance Term Structure ##
#############################
##Heston case
hrho1m=hcorr(log(S),VIX[,1]^2,21)
hrho3m=hcorr(log(S),VIX[,1]^2,63)

homega1m=volofvol(VIX[,1],21)
homega3m=volofvol(VIX[,1],63)
  
omega = c(VIX[1:62,6],homega3m[-(1:62)])
rho = c(rep(-0.1,62),hrho3m[-(1:62)])
omega1 = cbind(omega,omega,omega,omega,omega,omega)
rho1 = cbind(rho,rho,rho,rho,rho,rho)
VIXadded3m = VIX*sqrt(1 - omega1*rho1*tau1/2 + omega1^2*tau1^2/3/4)

hestonTSpar3m=varTSfitAltC(VIXadded3m^2,1,tau,2,6,2,TRUE,barrier=0.4,levelup=2,leveldn=0.6);

#summary(sqrt(hestonTSpar3m))
#summary(hestonTSpar3m2)
#plot(hestonTSpar3m[,1])
#plot(hestonTSpar3m[,2])

nu0=hestonTSpar3m[,1]
theta=hestonTSpar3m[,2]
kappa=hestonTSpar3m[,3]
kappa1=cbind(kappa,kappa,kappa,kappa,kappa,kappa)
kt=sweep(kappa1,MARGIN=2,tau,`*`)

#library(emdbook) #lambert W function
#a = (VIX[,3]^2-VIX[,6]^2)/(VIX[,1]^2-VIX[,6]^2)
#k = (a*lambertW_base(-exp(-1/a)/a)+1)/a/0.25
#plot(ts(EMA(pmin(k,rep(10,n))),10))
#lines(kappa,col=2)
#summary(pmin(k,rep(5,n)))

##Schoebel-Zhu case
hrho3msz=hcorr(log(S),VIX[,1],63)
homega3msz=sqrt(EMA(c(0,VIX[-1,1]-VIX[-n,1])^2,63)*252)
omegasz = c(VIX[1:62,6]/2,homega3msz[-(1:62)])
rhosz = c(rep(-0.1,62),hrho3msz[-(1:62)])
omegasz1 = cbind(omegasz,omegasz,omegasz,omegasz,omegasz,omegasz)
rhosz1 = cbind(rhosz,rhosz,rhosz,rhosz,rhosz,rhosz)

szTSpar = varTSfitAltC(VIXadded3m^2,2,tau,2,6,1.1,TRUE,max_it=8000,barrier=0.4,levelup=2,leveldn=0.6,omega=omegasz)

#summary(szTSpar)
#summary(sqrt(hestonTSpar3m))
#summary(shzTSpar)
#plot(szTSpar[,2])

## from I. Clark
## E[vol] = sqrt(E[var])*(1-1/8*varofvar/E[var]^2)
## W = vol^2*T
## from gatheral
## var[W_T] = theta*T * omega^2/kappa^2 + O(T)

#timefactor = 2*(nu0-theta)*(1-2*kt*ekt-ekt^2) + theta*(4*ekt-3+2*kt-ekt^2)
#varofvar = omega1^2/2/tau1^2/kappa1^3*timefactor
#convadj = (1-1/8*varofvar/VIXadded3m^4)
#convadj[,1]=pmax(convadj[,1],rep(0.00001,n))
#VolIndex = VIXadded3m*convadj

#########################
## Perfect calibration ##
#########################
#test.feller=feller(theta,homega3m,hrho3m,kappa,nu0);summary(test.feller>0)
#w=which(test.feller<0)
#homega3mf = homega3m
#homega3mf[w] = sqrt(0.9*2*kappa[w]*theta[w])

myparams = cbind(hestonTSpar3m[,2],omega,rho,hestonTSpar3m[,3],hestonTSpar3m[,1])
myparamsf = myparams
myparamsf[,2] = pmin(omega,sqrt(1.99*myparamsf[,1]*myparamsf[,4]))
myparamsz=cbind(szTSpar[,2],omegasz,rhosz,szTSpar[,3],szTSpar[,1])

t0=Sys.time(); perfectCal3 = opterCv(myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600); Sys.time()-t0;
t0=Sys.time(); perfectCal3f = opterCv(myparamsf,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=1); Sys.time()-t0;
t0=Sys.time(); perfectCal3sz = opterCv(myparamsz,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=2,err=TRUE,max_it=1600); Sys.time()-t0;
t0=Sys.time(); perfectCalOUOUsym = opterCv(as.matrix(myparamsz),mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=5,err=TRUE,max_it=1600); Sys.time()-t0;

summary(perfectCalOUOUsym)
perfectCalOUOUsym[170,]
which(perfectCalOUOUsym[,6]>1)
summary(myparamsz)
i=170
nlOpter(myparamsz[i,],mktPrice[i,],vega[i,],F[i,],K[i,],rd[i,],tau,5,1,5,TRUE,1600)
opterC(as.numeric(myparamsz[170,]),1,6,170,pair,opType=5,errType=1,modType=5,max_it=1600,err=TRUE)

perfectCal3alt = calibClean(perfectCal3,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,barrier=0.4,levelup=2)
perfectCal3falt = calibCleanF(perfectCal3f,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=1,barrier=0.4,levelup=2,leveldn=100,n)
perfectCal3falt[,2]=perfectCal3falt[,2]*0.999
test.fellerP=feller(perfectCal3falt);summary(test.fellerP>0)

plot(log(perfectCal3falt[,1]))
plot(log(perfectCal3falt[,2]))
plot(log(perfectCal3falt[,5]))
which((perfectCal3falt[,1])>1)

omegaP = perfectCal3alt[,2]
rhoP = perfectCal3alt[,3]

#perfectCal3faltf = calibNA(perfectCal3falt,myparams,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,c(2,1/2,1,2,1))

calibOutliers2(perfectCal3falt,0.4)

#Schoebel-Zhu
#perfect cal from historic estimates

perfectCal3szalt = calibClean(perfectCal3sz,myparamsz,mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=2,err=TRUE,max_it=1600,barrier=0.4,levelup=2)
calibOutliers2(perfectCal3sz,0.4)
plot((perfectCal3szalt[,2]))
plot((perfectCal3sz[,2]))

plot((perfectCalOUOUsym[,5]))

which((perfectCal3sz[,2])>0.7)

summary(perfectCal3sz)
szTSperf = szTSvar(cbind(perfectCal3sz[,5],perfectCal3sz[,1],perfectCal3sz[,4],perfectCal3sz[,2]))

summary(szTSperf)

#############################
## Specific point
#############################
plot(ts(omegaLmed))
m=350; a=1; b=6
omegaP[m]*rhoP[m]
omegaL0[m,];rhoL0[m,]
omegaL[m,]*rhoL[m,]
median(omegaLO[m,]);median(rhoLO[m,])
median(omegaL[m,]);median(rhoL[m,])
omegaD[m,];rhoD[m,];
opt=opter(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),"Heston","E5",a,b,m,pair="EURUSD");opt
opt2=opter(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),"Heston","E2",a,b,m,pair="EURUSD");opt2
opt2=opter(c(theta[m],homega3mmed[m],hrho3mmed[m],kappa[m],nu0[m]),"Heston","E2",a,b,m,pair="EURUSD");opt2
data.load1(datatable,m)

j=3
plotLines(y=vol,x=K,j)
yhat1=vol[j,3]+(K[j,]/S-1)*SS[m,j]+(K[j,]/S-1)^2*CC[m,j]/2
points(y=yhat1,x=K[j,],col=3);lines(y=yhat1,x=K[j,],col=3);

yy=fitsmile(c(theta[m],omegaD[m,1],rhoD[m,1],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
points(y=yy[j,],x=K[j,],col=2);lines(y=yy[j,],x=K[j,],col=2);

yyL=fitsmile(c(theta[m],omegaL[m,3],rhoL[m,3],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
yyL=fitsmile(c(theta[m],omegaLmed2[m],rhoLmed2[m],kappa[m],nu0[m]),S,vol,rd,rf,0,tau)
points(y=yyL[j,],x=K[j,],col=4);lines(y=yyL[j,],x=K[j,],col=4);

#############################
## 2-factor models tester  ##
#############################
## 1. OUOU ##
ouousymTS = varTSfitAltC(VIXadded3m^2/2,2,tau,2,6,1.1,TRUE,max_it=8000,barrier=0.4,levelup=2,leveldn=0.6,omega=omegaLmed/2/sqrt(2))
myparams1o = factor2from1(perfectCalOUOUsym[,-6],c(1,1,1,1,1))
lab2do1=tester(myparams1o,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2do1[[1]]))


## 2. Bates no Feller ##
myparams1 = factor2from1(perfectCal3alt[,-6],c(0.5,1,1,1,0.5))
lab2db1=tester(myparams1,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2db1[[1]]))

## 3. Bates Feller ##
myparams1f = factor2from1(perfectCal3falt[,-6],c(0.5,1,1,1,0.5))
myparams1f[,2] = pmin(myparams1f[,2],sqrt(1.99*myparams1f[,1]*myparams1f[,4]))
myparams1f[,7] = myparams1f[,2] 
lab2df1=tester(myparams1f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df1[[1]]))

## 3b. Bates Feller - other ##
myparams4f = factor2from1(perfectCal3alt[,-6],c(0.5,1,1,1,0.5))
myparams4f[,2] = pmin(myparams4f[,2],sqrt(1.99*myparams4f[,1]*myparams4f[,4]))
myparams4f[,7] = myparams4f[,2] 

myparamsff = myparams
myparamsff[,2] = pmin(omega,sqrt(0.99*myparamsff[,1]*myparamsff[,4]))
t0=Sys.time(); perfectCalFellerx2 = opterCv(as.matrix(myparamsff),mktPrice,vega,F,K,rd,tau,opType=5,errType=1,modType=1,err=TRUE,max_it=1600,fellerMust=2); Sys.time()-t0;
myparams5f = factor2from1(perfectCalFellerx2[,-6],c(1,1,1,1,1))

myparams6f = factor2from1(myparams,c(0.5,1,1,1,0.5))
myparams6f[,2] = pmin(myparams6f[,2],sqrt(1.99*myparams6f[,1]*myparams6f[,4]))
myparams6f[,7] = myparams6f[,2] 

lab2df4=tester(myparams4f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df4[[1]]))
lab2df5=tester(myparams5f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df5[[1]]))
lab2df6=tester(myparams6f,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2df6[[1]]))

## final testing
t0=Sys.time(); perfectCalOUOU1 = opterCv(myparams1o,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=4,err=TRUE,max_it=1600,logfile="nl_out1o.txt"); Sys.time()-t0;
t0=Sys.time(); perfectCalBates1 = opterCv(myparams1,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,logfile="nl_out1b.txt"); Sys.time()-t0;
t0=Sys.time(); perfectCalBates1f = opterCv(myparams1f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile="nl_out1f.txt"); Sys.time()-t0;

## other testing
t0=Sys.time(); perfectCalBates4f = opterCv(myparams4f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile="nl_out4f.txt"); Sys.time()-t0;
t0=Sys.time(); perfectCalBates5f = opterCv(myparams5f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile="nl_out5f.txt"); Sys.time()-t0;
t0=Sys.time(); perfectCalBates6f = opterCv(myparams6f,mktPrice,vega,F,K,rd,tau,opType=10,errType=1,modType=3,err=TRUE,max_it=1600,fellerMust=1,logfile="nl_out6f.txt"); Sys.time()-t0;

## 4. Simple starting points

simple = cbind(VIXadded3m[,6]^2/2,sqrt(2)*VIXadded3m[,6],0.99,2.1,VIXadded3m[,1]^2/2)
simplepars1 = cbind(simple,sweep(simple,2,c(1,1,-1,1,1),`*`))
lab2ds1=tester(simplepars1,F,K,rd,mktPrice,vega,tau,model=3);summary(ts(lab2ds1[[1]]))

simple = cbind(VIXadded3m[,6]/sqrt(2),sqrt(2)*VIXadded3m[,6]/2,0.99,1.1,VIXadded3m[,1]/sqrt(2))
simplepars2 = cbind(simple,sweep(simple,2,c(1,1,-1,1,1),`*`))
lab2ds2=tester(simplepars2,F,K,rd,mktPrice,vega,tau,model=4);summary(ts(lab2ds2[[1]]))

#############################
## Model error regression  ##
#############################
datatable2=datatable/100
level=datatable2$EURATM1M
slope=datatable2$EURATM2Y-datatable2$EURATM1M
curv=datatable2$EURATM2Y+datatable2$EURATM1M-2*datatable2$EURATM6M
RR1M=datatable2$EUR10P1M-datatable2$EUR10C1M
RR3M=datatable2$EUR10P3M-datatable2$EUR10C3M
RR1Y=datatable2$EUR10P1Y-datatable2$EUR10C1Y
RR2Y=datatable2$EUR10P2Y-datatable2$EUR10C2Y
BF1M=datatable2$EUR10P1M+datatable2$EUR10C1M-2*datatable2$EURATM1M
BF3M=datatable2$EUR10P3M+datatable2$EUR10C3M-2*datatable2$EURATM3M
BF1Y=datatable2$EUR10P1Y+datatable2$EUR10C1Y-2*datatable2$EURATM1Y
BF2Y=datatable2$EUR10P2Y+datatable2$EUR10C2Y-2*datatable2$EURATM2Y

lm1 = lm(perfectCal3alt[,6]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm1 = lm(perfectCal3sz[,6]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2a = lm(perfectCalBates1[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2c = lm(perfectCalBates1f[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)
lm2e = lm(perfectCalOUOU1[,11]/level ~ slope + curv + RR1M + RR2Y + BF1M + BF2Y)

library(texreg)
screenreg(list(lm1a, lm1b, lm2a, lm2b, lm2c, lm2d, lm2e, lm2f))
tab10 = texreg(list(m1, m2, m3, m4, m5, m6a, m6b), digits = 1, dcolumn = FALSE, booktabs = FALSE,use.packages = FALSE, label = "tab:10", caption = "Models explaining the probability of increases and decreases of implied volatility for the EURUSD 1M option in 1 quarter horizon",float.pos = "h!")
tab10 = gsub("hestonTSpar\\[, 1]", "$\\\\nu_{0,H}$", tab10)
tab10 = gsub("shzTSpar2\\[, 1]", "$\\\\nu_{0,S}^2$", tab10)
tab10 = gsub("hestonTSpar\\[, 2]", "$\\\\theta_H$", tab10)
tab10 = gsub("shzTSpar2\\[, 2]", "$\\\\theta_S^2$", tab10)
tab10 = gsub("kappa2shz", "$\\\\xi_S$", tab10)
tab10 = gsub("kappa2", "$\\\\xi_H$", tab10)
write.table(c(tab10,AUROC_stats),"C:/Users/mikolaj.labedzki/Downloads/dane/FX/r3.t10.txt")

#############################
## Tables Printing         ##
#############################
tab0 = rbind(summarypy(hestonTSpar3m[,1]),summarypy(hestonTSpar3m[,2]),summarypy(hestonTSpar3m[,3]),summarypy(sqrt(hestonTSpar3m[,4]/5)))
rownames(tab0)=c("$\\nu_0$","$\\theta$","$\\kappa$","MAE")
cfa0.table <- xtable(tab0)
digits(cfa0.table) <- 4
write.textable(cfa0.table,"r5.t1.txt")  

dat1=as.data.frame(cbind(perfectCalOUOU1[,11],perfectCalBates1[,11]))
dat1=as.data.frame(cbind(perfectCalBates1f[,11],perfectCalOUOU1[,11]))
dat1=as.data.frame(cbind(lab2df1[[1]],lab2do1[[1]]))
dat1=as.data.frame(cbind(par5d,par5l))
ndat1=reshape(dat1, direction="long", varying=list(names(dat1)), v.names="Value")
oneway.test(Value~time, ndat1, var.equal = FALSE)
t.test(Value~time, data = ndat1, paired = TRUE)

tab6=rbind(summarypy(lab2db1[[1]]),summarypy(lab2do1[[1]]),summarypy(lab2df1[[1]]))
modelsmethods=c("Bates/2-stage","OUOU/2-stage","BatesFeller/2-stage")
rownames(tab6)=modelsmethods
cfa6.table <- xtable(tab6)
digits(cfa6.table) <- 4
write.textable(cfa6.table,"r5.t6.txt")

tab7=rbind(summarypy(perfectCalBates1[,11]),summarypy(perfectCalOUOU1[,11]),summarypy(perfectCalBates1f[,11]))
rownames(tab7)=modelsmethods
cfa7.table <- xtable(tab7)
digits(cfa7.table) <- 5
write.textable(cfa7.table,"r5.t7.txt")

#############################
## Plot Printing           ##
#############################

myaxis = seq(as.Date("2010-07-22"), as.Date("2015-08-31"), by = 91)
datatable$Date = as.Date(datatable$Date)
mywidth = 15
myheight = 10

pdf("../images/rhos.pdf", width=mywidth, height=myheight)
plot(rhoP~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", type = "b", lty=1, pch=1, ylab="rho", xlab="Date", ylim=c(-0.75,0.2), cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,hrho3m,col=2, pch = 2, cex=0.5);lines(datatable$Date,hrho3m, col=2, lty = 2, cex=0.5);
abline(h=0,lty=2)
legend("topleft", c("Calib.","Hist. 3m"), lwd=1.4, lty=1:2, pch=c(1,2), col=c(1,2), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/omegas.pdf", width=mywidth, height=myheight)
plot(omegaP~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", type = "b", lty=1, pch=1, col=1, ylab="omega", ylim=c(0.0,0.85), xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,homega3m,col=2, pch = 2, cex=0.5);lines(datatable$Date,homega3m,col=2, lty = 2, cex=0.5)
legend("topleft", c("Calib.","Hist. 3m"), lwd=1.4, lty=1:2, pch=c(1,2), col=c(1,2), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

pdf("../images/durrvslab.pdf", width=mywidth, height=myheight)
plot(hist3m[[1]]perfectCal3alt[,6]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", pch=1, lty=1, t="b", ylab=c("RMSE"), xlab="Date", ylim=c(0.0005,0.028), cex=0.5, cex.lab=1.4, cex.sub=1.4)
points(datatable$Date,perfectCal3alt[,6], col=2, pch=2, cex=0.5);lines(datatable$Date,hist3m[[1]], col=2, lty=2, cex=0.5)
points(datatable$Date,perfectCal3falt[,6], col=3, pch=3, cex=0.5);lines(datatable$Date,perfectCal3falt[,6], col=3, lty=3, cex=0.5)
points(datatable$Date,perfectCal3szalt[,6], col=4, pch=4, cex=0.5);lines(datatable$Date,perfectCal3szalt[,6], col=4, lty=4, cex=0.5)
legend("topleft", c("Hist.-Heston","Calib.-Heston", "Calib.-Heston-Feller", "Calib.-Schöbel-Zhu"), lwd=1.4, pch=1:4, lty=1:4, col=c(1,2,3,4), bty="n", cex=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

library(latticeExtra)
library(reshape)
wide = as.data.frame(cbind(c(" 1m"," 2m"," 3m"," 6m","12m","24m"),hist3m[[2]]))
colnames(wide) = c("tenor","-10p","-25p","50c","25c","10c")
wide.out = melt(wide, measure.vars = c("-10p","-25p","50c","25c","10c"), variable.name = "delta", value.name = "error")
wide.out$value=as.numeric(as.character(wide.out$value))

mat=matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4)
mat=matrix(c(1,0.2,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4)
pdf("../images/better1.pdf", width=7, height=7)
print(cloud(value~tenor+variable, wide.out, main="Errors in Hist3m", panel.3d.cloud=panel.3dbars, distance = 0.3, R.mat = mat, xbase=0.5, ybase=0.5, col.facet='green', scales=list(arrows=FALSE, col=1), par.settings = list(axis.line = list(col = "transparent")), zlim=c(0,0.006), aspect = c(1,0.3,1), par.box = list(col="grey"), drape=FALSE, shade = FALSE, ylab=c("delta"), zlab=c("RMSE")), split = c(1,1,2,1), more = TRUE)
dev.off()

pdf("../images/rmse2d_start.pdf", width=mywidth, height=myheight)
plot(lab2db1[[1]]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=1, lty=1, ylab=c("RMSE"), xlab="Date", ylim=c(0.001,0.031), cex=0.5, cex.lab=1.4, cex.sub=1.4, col=1)
points(datatable$Date,lab2do1[[1]], col=2, pch=2, cex=0.5);lines(datatable$Date,lab2do1[[1]], col=2, lty=2, cex=0.5)
points(datatable$Date,lab2df1[[1]], col=3, pch=3, cex=0.5);lines(datatable$Date,lab2df1[[1]], col=3, lty=3, cex=0.5)
legend("topright", modelsmethods, pch=(1:3), lty=(1:3), col=1:3, bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

# Part 2: Comparison of better starting point and full calibration vs 2phase calibration simple starting point vs Bates old full calib

pdf("../images/rmse2d_end.pdf", width=mywidth, height=myheight)
plot(perfectCalBates1f[,11]~datatable$Date, xlim=c(as.Date("2010-09-10"),as.Date("2015-07-10")), xaxt="n", t="b", pch=1, lty=1, ylab=c("RMSE"), xlab="Date", cex=0.5, cex.lab=1.4, cex.sub=1.4, col=1, ylim=c(0.0003,0.0142))
points(datatable$Date,perfectCalOUOU1[,11], col=2, pch=2, cex=0.5);lines(datatable$Date,perfectCalOUOU1[,11], col=2, lty=2, cex=0.5)
points(datatable$Date,perfectCalBates1[,11], col=3, pch=3, cex=0.5);lines(datatable$Date,perfectCalBates1[,11], col=3, lty=3, cex=0.5)
legend("topright", modelsmethods, pch=(1:3), lty=(1:3), col=1:3, bty="n", cex=1.4, lwd=1.4)
axis(1, myaxis, format(myaxis, "%b %y"), cex.axis = 1)
grid()
dev.off()

convergence1=read.table("proper/nl_out1b.txt",sep=";",fill=TRUE)
convergence1o=read.table("proper/nl_out1o.txt",sep=";",fill=TRUE)
convergence1f=read.table("proper/nl_out1f.txt",sep=";",fill=TRUE)

convergence1o=t(matrix(rep(as.numeric(unlist(t(convergence1o))),6),161,n+5))[-((n+1):(n+5)),]
convergence1f=t(matrix(rep(as.numeric(unlist(t(convergence1f))),6),161,n+5))[-((n+1):(n+5)),]
convergence1=t(matrix(rep(as.numeric(unlist(t(convergence1))),6),161,n+5))[-((n+1):(n+5)),]

#fix(convergence1)
convergence1[,161] = as.numeric(rep(NA,n))
convergence1o[,161] = as.numeric(rep(NA,n))
convergence1f[,161] = as.numeric(rep(NA,n))
convergence1=addatend(convergence1,perfectCalBates1[,11])
convergence1o=addatend(convergence1o,perfectCalOUOU1[,11])
convergence1f=addatend(convergence1f,perfectCalBates1f[,11])
convergence1=cbind(lab2db1[[1]],convergence1)
convergence1o=cbind(lab2do1[[1]],convergence1o)
convergence1f=cbind(lab2df1[[1]],convergence1f)

pdf("../images/rmse2d_conv_log.pdf", width=mywidth, height=myheight)
x=(0:161)*10
x1=1:162
plot(x,log(colMeans(convergence1[,x1])), t="b", lty=1, pch=1,col=1, ylab=c("Natural logarithm of mean RMSE"), xlab="Number of iterations", cex=0.5, cex.lab=1.4, cex.sub=1.4, ylim=c(-6.7,-4.5)) # ylim=c(0.0015,0.0111))ylim=c(-6.7,-4.5)
points(x,log(colMeans(convergence1o[,x1])), col=2, pch=2, cex=0.5);lines(x,log(colMeans(convergence1o[,x1])),col=2,lty=2)
points(x,log(colMeans(convergence1f[,x1])), col=3, pch=3, cex=0.5);lines(x,log(colMeans(convergence1f[,x1])),col=3,lty=3)
legend("topright", modelsmethods, pch=(1:3), lty=(1:3), col=1:3, bty="n", cex=1.4, lwd=1.4)
grid()
dev.off()
