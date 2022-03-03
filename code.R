########## Jointly modeling rating responses and times with fuzzy numbers: An application to psychometric data ###############



########## Data and functions import ###########
rm(list=ls())
source("C:/Users/NC/Downloads/MDPI/utils.R")
load("C:/Users/NC/Downloads/MDPI/data.R")



########## Data analysis ###########

#IRTree model
M=4;N=M-1; I=NROW(D); J=NCOL(D) 
Tm = matrix(c(0,NA,NA,
              1,0,NA,
              1,1,0,
              1,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = c("D/A","Am/Ae","Am"); print(Tm)

Ydata = irtrees::dendrify(as.matrix(D),Tm)
mod_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Ydata) 

#Model parameters estimation
library(lme4)
Alpha_est = matrix(summary(mod_glm)$coefficients$cond[,1],J,N)
Eta_est = ranef(mod_glm)$cond$person


#fIRTree data of depression scores
fuzzyx = getFuzzyNumbers(Eta_est,Alpha_est,Tm, RTS)
C_est = matrix(fuzzyx$Yfuzzy$C,I,J,byrow = TRUE)
L_est = matrix(fuzzyx$Yfuzzy$L,I,J,byrow = TRUE)
R_est = matrix(fuzzyx$Yfuzzy$R,I,J,byrow = TRUE)
W_est = matrix(fuzzyx$Yfuzzy$W,I,J,byrow = TRUE)

#Rescaling of fuzzy data
C_est=mapply(function(j){mapply(function(i) ((4-1)/(max(C_est[,j])))*(C_est[i,j]-max(C_est[,j]))+4,1:NROW(C_est))},1:J)
L_est=mapply(function(j){mapply(function(i) ((4-1)/(max(L_est[,j])))*(L_est[i,j]-max(L_est[,j]))+4,1:NROW(L_est))},1:J)
R_est=mapply(function(j){mapply(function(i) ((4-1)/(max(R_est[,j])))*(R_est[i,j]-max(R_est[,j]))+4,1:NROW(R_est))},1:J)

#Composite indicator for depression: mean for each fuzzy number parameter
meanL=apply(L_est, 1,mean)
meanR=apply(R_est, 1,mean)
meanC=apply(C_est, 1,mean)
meanW=apply(W_est, 1,mean)


#Plots of Depression indicator fuzzy parameters
library(tikzDevice)
#Plot of centers
tikz('C:/Users/NC/Downloads/MDPI/fig3a.tex',width=8.5,height=10.5)
par(mai=c(1,2,0.5,2),mgp=c(3,1,-0.70))
hist(meanC,xlab = "",ylab = "",main = "", col = "#414a6d", breaks = 11,cex.axis=2.75)
title(xlab = "Centers",cex.lab=3,line = 4)
title(ylab = "Frequency",cex.lab=3,line = 4)
dev.off()
#Concentration of values close to one
around_one=meanW[meanW>=0.95&meanW<=1.05]
perc_around_one=length(around_one)/length(meanW)*100

#Plot of right spreads
tikz('C:/Users/NC/Downloads/MDPI/fig3b.tex',width=8.5,height=10.5)
par(mai=c(1,2,0.5,2),mgp=c(3,1,-0.70))
hist(meanR-meanC,xlab = "",ylab = "",main = "", col = "#414a6d", breaks = 11,cex.axis=2.75)
title(xlab = "Left spreads",cex.lab=3,line = 4)
title(ylab = "Frequency",cex.lab=3,line = 4)
dev.off()

#Plot of left spreads
tikz('C:/Users/NC/Downloads/MDPI/fig3c.tex',width=8.5,height=10.5)
par(mai=c(1,2,0.5,2),mgp=c(3,1,-0.70))
hist(meanC-meanL,xlab = "",ylab = "",main = "", col = "#414a6d", breaks = 11,cex.axis=2.75)
title(xlab = "Right spreads",cex.lab=3,line = 4)
title(ylab = "Frequency",cex.lab=3,line = 4)
dev.off()

#Plot of intensification parameter
tikz('C:/Users/NC/Downloads/MDPI/fig3d.tex',width=8.5,height=10.5)
par(mai=c(1,2,0.5,2),mgp=c(3,1,-0.70))
hist(meanW,xlab = "",ylab = "",main = "", col = "#414a6d", breaks = 11,cex.axis=2.75)
title(xlab = "Intensification parameters",cex.lab=3,line = 4)
title(ylab = "Frequency",cex.lab=3,line = 4)
dev.off()




######## Traditional data analysis #########


##### Normal Linear Model

meanDisc=apply(D,1,mean) #mean of crisp ratings for each subject over items
fit=lm(meanDisc~factor(data$religion)+idx_es+factor(data$education))
null=lm(meanDisc~1) #null model
summary(fit)
confint(fit)

#pseudo-R^2 calculation:
n=NROW(meanDisc)
omeg=2*(as.numeric(logLik(fit))-as.numeric(logLik(null)))
lamb=(1/n)*as.numeric(logLik(null))
pseudoR2_normal=-(omeg*(1-lamb))/((omeg+n)*lamb)


##### Log-Normale Linear Model

meanRTS=apply(RTS,1,mean) #mean of RTs for each subject over items
modRTS_full=lm(I(log(meanRTS))~factor(data$religion)+idx_es+factor(data$education))
modRTS_null=lm(I(log(meanRTS))~1) #null model
summary(modRTS_full)
confint(modRTS_full)

#pseudo-R^2 calculation:
omeg=2*(as.numeric(logLik(modRTS_full))-as.numeric(logLik(modRTS_null)))
lamb=(1/n)*as.numeric(logLik(modRTS_null))
pseudoR2_log=-(omeg*(1-lamb))/((omeg+n)*lamb)




######## Fuzzy Normal Linear Model ###########

n = NROW(D)
X=model.matrix(~factor(data$religion)+idx_es+factor(data$education))
J=NCOL(X)-1 #n of predictors

X_null=matrix(1,nrow=n,ncol=1) #null model
J_null=NCOL(X_null)-1 

#Response variable: depression composite indicator
outcome=matrix(c(meanL,meanC,meanR,meanW),nrow=n)

#Model fitting
resIna = optim(par = rep(1,J+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J+1),1e-10), upper = rep(Inf,J+2), y = outcome, X = X, control = list(trace = 3, maxit = 1000), hessian = TRUE)
resIna_null = optim(par = rep(1,J_null+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J_null+1),1e-10), upper = rep(Inf,J_null+2), y = outcome, X = X_null, control = list(trace = 3, maxit = 1000), hessian = TRUE)

b_est = resIna$par[-(J+1+1)] #coefficients
sigma_beta = sqrt(diag(solve(resIna$hessian)))[-(J+1+1)] #standard errors
sigmay_est = resIna$par[(J+1+1)] 

CI = cbind(b_est-qt(1-0.05/2,n-J-1)*sigma_beta,
           b_est+qt(1-0.05/2,n-J-1)*sigma_beta)

estimates=cbind(b_est,sigma_beta,CI)

#pseudo-R^2 calculation:
llik_null=likelihood_dombi_fn(pars=resIna_null$par,X=X_null,y=outcome)
llik_full=likelihood_dombi_fn(pars=resIna$par,X=X,y=outcome)
omeg=2*(llik_full-llik_null)
lamb=(1/n)*llik_null
pseudoR2_fuzzy=-(omeg*(1-lamb))/((omeg+n)*lamb)


#y_star = defuzzyfication of depression
y_hat = X%*%b_est 
y_star=mapply(function(i) integrate(function(x){x*dombi_fn(x = x,l = meanL[i] ,m = meanC[i],r = meanR[i],w = meanW[i])$mu*dnorm(x,X[i,]%*%b_est,sd = sigmay_est)},lower = meanL[i]-1e-02,upper = meanR[i]+1e-02)$value/integrate(function(x){dombi_fn(x = x,l = meanL[i],m = meanC[i],r = meanR[i],w = meanW[i])$mu*dnorm(x,mean = X[i,]%*%b_est,sd = sigmay_est)},lower = meanL[i]-1e-02,upper = meanR[i]+1e-02)$value,1:n)
y_star=matrix(y_star,160,1)

#Residuals
r=y_star-y_hat
summary(r)

#Plot of fitted Fuzzy Normal Linear Model
ips=cbind(meanL,meanC,meanR)
ics=data.frame("religiousness"=factor(data$religion),"emotional_stability"=idx_es,"university"=factor(data$education))
jcs=ics$emotional_stability
len=501
db = matrix(NA,160,len); dbnmu = matrix(NA,160,len)
for(i in 1:160){
  db[i,] = dombi_fn(NULL,meanL[i],meanC[i],meanR[i],meanW[i],len)$mu
  dbnmu[i,] = scales::rescale(db[i,],from = range(db[i,]),to=c(jcs[i],jcs[i]+1))
}
dbnx = matrix(NA,160,len)
for(i in 1:160){dbnx[i,] = dombi_fn(NULL,meanL[i],meanC[i],meanR[i],meanW[i],len)$x}

tikzDevice::tikz(file='C:/Users/NC/Downloads/MDPI/fig5.tex',width=6,height=5)
par(mai=c(1.15, 0.85, 0.65, 0.15))
cols=c("palegreen4","red4","blue4","goldenrod4")
plot(X[1:NROW(meanC),3],y=meanC[1:NROW(meanC)],col="white", lwd=0,xlim=c(min(jcs),max(jcs+1)),ylim = c(min(meanC)-0.1,max(meanC)), bty="n", main = "",xlab = "", ylab = "",cex.axis=1.1)
title(xlab = "emotional$\\_$stability",cex.lab=1.25,line = 2.5)
title(ylab = "Fuzzy ratings",cex.lab=1.25,line = 3)
set.seed(8542)
for(i in sample(1:160,25)){
  lines(dbnmu[which(ics$religiousness=="1"&ics$university=="1")[i],],dbnx[which(ics$religiousness=="1"&ics$university=="1")[i],], col = cols[1],lty=5,lwd = 2)
  points(y=meanC[which(ics$religiousness=="1"&ics$university=="1")[i]],x=max(dbnmu[which(ics$religiousness=="1"&ics$university=="1")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="1"&ics$university=="1")[i],1],x0=jcs[which(ics$religiousness=="1"&ics$university=="1")[i]],ips[which(ics$religiousness=="1"&ics$university=="1")[i],3],x1=jcs[which(ics$religiousness=="1"&ics$university=="1")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="2"&ics$university=="1")[i],],dbnx[which(ics$religiousness=="2"&ics$university=="1")[i],], col = cols[2],lty=5,lwd = 2)
  points(y=meanC[which(ics$religiousness=="2"&ics$university=="1")[i]],x=max(dbnmu[which(ics$religiousness=="2"&ics$university=="1")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="2"&ics$university=="1")[i],1],x0=jcs[which(ics$religiousness=="2"&ics$university=="1")[i]],ips[which(ics$religiousness=="2"&ics$university=="1")[i],3],x1=jcs[which(ics$religiousness=="2"&ics$university=="1")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="1"&ics$university=="3")[i],],dbnx[which(ics$religiousness=="1"&ics$university=="3")[i],], col = cols[3],lty=5,lwd = 2)
  points(y=meanC[which(ics$religiousness=="1"&ics$university=="3")[i]],x=max(dbnmu[which(ics$religiousness=="1"&ics$university=="3")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="1"&ics$university=="3")[i],1],x0=jcs[which(ics$religiousness=="1"&ics$university=="3")[i]],ips[which(ics$religiousness=="1"&ics$university=="3")[i],3],x1=jcs[which(ics$religiousness=="1"&ics$university=="3")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="2"&ics$university=="3")[i],],dbnx[which(ics$religiousness=="2"&ics$university=="3")[i],], col = cols[4],lty=5,lwd = 2)
  points(y=meanC[which(ics$religiousness=="2"&ics$university=="3")[i]],x=max(dbnmu[which(ics$religiousness=="2"&ics$university=="3")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="2"&ics$university=="3")[i],1],x0=jcs[which(ics$religiousness=="2"&ics$university=="3")[i]],ips[which(ics$religiousness=="2"&ics$university=="3")[i],3],x1=jcs[which(ics$religiousness=="2"&ics$university=="3")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
}
abline(a=b_est[1],b=b_est[3],col = cols[1], lwd = 2)
abline(a=b_est[1]+b_est[2],b=b_est[3],col = cols[2], lwd = 2)
abline(a=b_est[1]+b_est[4],b=b_est[3],col = cols[3], lwd = 2)
abline(a=b_est[1]+b_est[2]+b_est[4],b=b_est[3],col = cols[4], lwd = 2)

add_legend("bottom",fill = cols,legend = c("religiousness:no,university:no","religiousness:yes,university:no","religiousness:no,university:yes","religiousness:yes,university:yes"),border = FALSE,bty = "n",ncol = 2,cex=1.2)
dev.off()





##### #Fuzzy Normal linear model with W={1,...,1} #####
#Response variable
meanW_one=meanW/meanW
outcome2=matrix(c(meanL,meanC,meanR,meanW_one),nrow=n)

#Model fitting
resIna2 = optim(par = rep(1,J+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J+1),1e-10), upper = rep(Inf,J+2), y = outcome2, X = X, control = list(trace = 3, maxit = 1000), hessian = TRUE)
resIna_null2 = optim(par = rep(1,J_null+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J_null+1),1e-10), upper = rep(Inf,J_null+2), y = outcome2, X = X_null, control = list(trace = 3, maxit = 1000), hessian = TRUE)

#pseudo-R^2
llik_null2=likelihood_dombi_fn(pars=resIna_null2$par,X=X_null,y=outcome2)
llik_full2=likelihood_dombi_fn(pars=resIna2$par,X=X,y=outcome2)

omeg2=2*(llik_full2-llik_null2)
lamb2=(1/n)*llik_null2
pseudoR2_fuzzy2=-(omeg2*(1-lamb2))/((omeg2+n)*lamb2)
