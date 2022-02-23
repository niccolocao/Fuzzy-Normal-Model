########## DASS ###############
rm(list=ls())
source("C:/Users/NC/Downloads/MDPI_template/MDPI_template/utils.R")
data <- read.delim("C:/Users/NC/Downloads/DASS_data_21.02.19/DASS_data_21.02.19/data.csv")
data$ID=1:NROW(data)

#Item selection of Depression subscale
D=data[,c(7,13,28,37,46,49,61,70,76,91,100,109,112,124,NCOL(data))]

#Sample estraction
set.seed(3400404)
sub=sample(D$ID,250,replace = FALSE)
D=D[sub,]
data=data[sub,]


########### Preprocessing #########
#Recodification of Religiousness and University predictors
table(data$religion)
table(data$education)

data$religion[data$religion==2]=1
data$religion[data$religion==3]=2
data$religion[data$religion==4]=2
data$religion[data$religion==5]=2
data$religion[data$religion==6]=2
data$religion[data$religion==7]=2
data$religion[data$religion==8]=2
data$religion[data$religion==9]=2
data$religion[data$religion==10]=2
data$religion[data$religion==11]=2
data$religion[data$religion==12]=2

D=D[-c(which(data$religion==0)),]
data=data[-c(which(data$religion==0)),]

D=D[-c(which(data$gender==3)),] #just considering binary gender
data=data[-c(which(data$gender==3)),]

D=D[-c(which(data$age==13)),] #administration for over 14 years old
data=data[-c(which(data$age==13)),]

data$education[data$education==4]=3
data$education[data$education==2]=1
table(data$education)




#Response times for each item
RTS = data[,c(7,13,28,37,46,49,61,70,76,91,100,109,112,124)+2]
RTS=RTS/1000


#Outliers removing in RTs: over 2 sd from the mean of RTs
outliers=mapply(function(j) abs(RTS[,j]-mean(RTS[,j]))>2*sd(RTS[,j]),1:NCOL(RTS))
RTS=RTS[-which(outliers[,1]==TRUE),]
RTS=RTS[-which(outliers[,2]==TRUE),]
RTS=RTS[-which(outliers[,3]==TRUE),]
RTS=RTS[-which(outliers[,4]==TRUE),]
RTS=RTS[-which(outliers[,5]==TRUE),]
RTS=RTS[-which(outliers[,6]==TRUE),]
RTS=RTS[-which(outliers[,7]==TRUE),]
RTS=RTS[-which(outliers[,8]==TRUE),]
RTS=RTS[-which(outliers[,9]==TRUE),]
RTS=RTS[-which(outliers[,10]==TRUE),]
RTS=RTS[-which(outliers[,11]==TRUE),]
RTS=RTS[-which(outliers[,12]==TRUE),]
RTS=RTS[-which(outliers[,13]==TRUE),]
RTS=RTS[-which(outliers[,14]==TRUE),]



D=D[-which(outliers[,1]==TRUE),];data=data[-which(outliers[,1]==TRUE),]
D=D[-which(outliers[,2]==TRUE),];data=data[-which(outliers[,2]==TRUE),]
D=D[-which(outliers[,3]==TRUE),];data=data[-which(outliers[,3]==TRUE),]
D=D[-which(outliers[,4]==TRUE),];data=data[-which(outliers[,4]==TRUE),]
D=D[-which(outliers[,5]==TRUE),];data=data[-which(outliers[,5]==TRUE),]
D=D[-which(outliers[,6]==TRUE),];data=data[-which(outliers[,6]==TRUE),]
D=D[-which(outliers[,7]==TRUE),];data=data[-which(outliers[,7]==TRUE),]
D=D[-which(outliers[,8]==TRUE),];data=data[-which(outliers[,8]==TRUE),]
D=D[-which(outliers[,9]==TRUE),];data=data[-which(outliers[,9]==TRUE),]
D=D[-which(outliers[,10]==TRUE),];data=data[-which(outliers[,10]==TRUE),]
D=D[-which(outliers[,11]==TRUE),];data=data[-which(outliers[,11]==TRUE),]
D=D[-which(outliers[,12]==TRUE),];data=data[-which(outliers[,12]==TRUE),]
D=D[-which(outliers[,13]==TRUE),];data=data[-which(outliers[,13]==TRUE),]
D=D[-which(outliers[,14]==TRUE),];data=data[-which(outliers[,14]==TRUE),]





#Ten Item Presonality Inventory cleaning
TIPI=data[,c(132:141,NCOL(data))] 
TIPI[TIPI==0]=NA
TIPI$ID=1:NROW(TIPI)
TIPI=na.omit(TIPI)
D=D[c(TIPI$ID),]
data=data[c(TIPI$ID),]
RTS=RTS[c(TIPI$ID),]
TIPI=TIPI[,-NCOL(TIPI)]
D=D[,-NCOL(D)]

#Compund indicator of emotional_stability
TIPI[,c(6,2,8,4,10)]=abs(TIPI[,c(6,2,8,4,10)]-8)#item reversed
es=TIPI[,c(4,9)]
alpha_es=psych::alpha(es)
alpha_es$total[1]
idx_es = unlist(mapply(function(i) alpha_es$total[1]*sum(es[i,]) + (1-alpha_es$total[1])*mean(na.omit(rowSums(es))),1:NROW(es)))


########## Data analysis ###########
#IRTree models
M=4;N=M-1; I=NROW(D); J=NCOL(D) 
Tm = matrix(c(0,NA,NA,
              1,0,NA,
              1,1,0,
              1,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = c("D/A","Am/Ae","Am"); print(Tm)

Ydata = irtrees::dendrify(as.matrix(D),Tm)
mod1_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Ydata) 
mod2_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+1|person), family = binomial, data = Ydata)

AIC(mod1_glm,mod2_glm)

require(lme4)
Alpha_est = matrix(summary(mod1_glm)$coefficients$cond[,1],J,N)
Eta_est = ranef(mod1_glm)$cond$person



#fIRTree data
fuzzyx = getFuzzyNumbers(Eta_est,Alpha_est,Tm, RTS)
C_est = matrix(fuzzyx$Yfuzzy$C,I,J,byrow = TRUE)
L_est = matrix(fuzzyx$Yfuzzy$L,I,J,byrow = TRUE)
R_est = matrix(fuzzyx$Yfuzzy$R,I,J,byrow = TRUE)
W_est = matrix(fuzzyx$Yfuzzy$W,I,J,byrow = TRUE)


#Rescaling of fuzzy data
C_est=mapply(function(j){mapply(function(i) ((4-1)/(max(C_est[,j])))*(C_est[i,j]-max(C_est[,j]))+4,1:NROW(C_est))},1:J)
L_est=mapply(function(j){mapply(function(i) ((4-1)/(max(L_est[,j])))*(L_est[i,j]-max(L_est[,j]))+4,1:NROW(L_est))},1:J)
R_est=mapply(function(j){mapply(function(i) ((4-1)/(max(R_est[,j])))*(R_est[i,j]-max(R_est[,j]))+4,1:NROW(R_est))},1:J)

#Calculation of the mean for each fuzzy number parameter
menL=apply(L_est, 1,mean)
menR=apply(R_est, 1,mean)
menC=apply(C_est, 1,mean)
menW=apply(W_est, 1,mean)


#Plot of intensification parameter
library(tikzDevice)
tikz('C:/Users/NC/Downloads/MDPI_template/MDPI_template/fig3a.tex',width=8.5,height=10.5)
par(mfrow=c(1,2))
hist(menW,xlab= "Intensification parameter",main = "",cex.lab=1.5)
dev.off()

#Plot of the difference between lower-bounds and centroids
tikz('C:/Users/NC/Downloads/MDPI_template/MDPI_template/fig3b.tex',width=8.5,height=10.5)
par(mfrow=c(1,2))
hist(menR-menC,xlab= "Lower Bounds-Centroids differences",main = "",cex.lab=1.5)
dev.off()

#Plot of the difference between centroids and upper-bounds
tikz('C:/Users/NC/Downloads/MDPI_template/MDPI_template/fig3c.tex',width=8.5,height=10.5)
par(mfrow=c(1,2))
hist(menC-menL,xlab= "Centroids-Upper Bounds differences",main = "",cex.lab=1.5)
dev.off()


######## Traditional data analysis procedure #########


#Normal Linear Model on the mean of crisp ratings for each subject over items
menDisc=apply(D,1,mean)
fit=lm(menDisc~factor(data$religion)+idx_es+factor(data$education))
null=lm(menDisc~1) #null model
summary(fit)
confint(fit)
#pseudo-R^2
omeg=2*(as.numeric(logLik(fit))-as.numeric(logLik(null)))
lamb=(1/n)*as.numeric(logLik(null))
pseudoR2_normal=-(omeg*(1-lamb))/((omeg+n)*lamb)

#Log-Normale Linear Model on the mean of RTs for each subject over items
menRTS=apply(RTS,1,mean)
modRTS_full=lm(I(log(menRTS))~factor(data$religion)+idx_es+factor(data$education))
modRTS_null=lm(I(log(menRTS))~1) #null model
summary(modRTS_full)
confint(modRTS_full)
#pseudo-R^2
omeg=2*(as.numeric(logLik(modRTS_full))-as.numeric(logLik(modRTS_null)))
lamb=(1/n)*as.numeric(logLik(modRTS_null))
pseudoR2_log=-(omeg*(1-lamb))/((omeg+n)*lamb)


######## Fuzzy Normal Linear Model ###########
n = NROW(D)
X=model.matrix(~factor(data$religion)+idx_es+factor(data$education))

X_null=matrix(1,nrow=n,ncol=1) #null model
J_null=NCOL(X_null)-1 
J=NCOL(X)-1 #predictors

#Response variable
outcome=matrix(c(menL,menC,menR,menW),nrow=n)

#Model fitting
resIna = optim(par = rep(1,J+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J+1),1e-10), upper = rep(Inf,J+2), y = outcome, X = X, control = list(trace = 3, maxit = 1000), hessian = TRUE)
resIna_null = optim(par = rep(1,J_null+2), fn = inaccuracy_dombi_fn, method = "L-BFGS-B", lower = c(rep(-Inf,J_null+1),1e-10), upper = rep(Inf,J_null+2), y = outcome, X = X_null, control = list(trace = 3, maxit = 1000), hessian = TRUE)

b_est = resIna$par[-(J+1+1)] #coefficients
sigma_beta = sqrt(diag(solve(resIna$hessian)))[-(J+1+1)] #standard errors


CI = cbind(b_est-qt(1-0.05/2,n-J-1)*sigma_beta,
           b_est+qt(1-0.05/2,n-J-1)*sigma_beta)

estimates=cbind(b_est,sigma_beta,CI)

#pseudo-R^2
llik_null=likelihood_dombi_fn(pars=resIna_null$par,X=X_null,y=outcome)
llik_full=likelihood_dombi_fn(pars=resIna$par,X=X,y=outcome)

omeg=2*(llik_full-llik_null)
lamb=(1/n)*llik_null
pseudoR2_fuzzy=-(omeg*(1-lamb))/((omeg+n)*lamb)


#y_star = defuzzyfication of the outcome variable
y_hat = X%*%b_est 
y_star=mapply(function(i) integrate(function(x){x*dombi_fn(x = x,l = menL[i] ,m = menC[i],r = menR[i],w = menW[i])$mu*dnorm(x,X[i,]%*%b_est,sd = sigmay_est)},lower = menL[i]-1e-02,upper = menR[i]+1e-02)$value/integrate(function(x){dombi_fn(x = x,l = menL[i],m = menC[i],r = menR[i],w = menW[i])$mu*dnorm(x,mean = X[i,]%*%b_est,sd = sigmay_est)},lower = menL[i]-1e-02,upper = menR[i]+1e-02)$value,1:n)
y_star=matrix(y_star,160,1)

#Residuals
H = X%*%solve(t(X)%*%X)%*%t(X)
r=y_star-y_hat
stan = mapply(function(i) r[i]/sqrt(1-diag(H)[i]), 1:NROW(r))
stud = mapply(function(i) r[i]/(((sigmay_est)^2)*sqrt(1-diag(H)[i])), 1:NROW(r))
hist(r)


#Plot of fitted Fuzzy Normal Linear Model
ips=cbind(menL,menC,menR)
ics=data.frame("religiousness"=factor(data$religion),"emotional_stability"=idx_es,"university"=factor(data$education))
jcs=ics$emotional_stability
len=501
db = matrix(NA,160,len); dbnmu = matrix(NA,160,len)
for(i in 1:160){
  db[i,] = dombi_fn(NULL,menL[i],menC[i],menR[i],menW[i],len)$mu
  dbnmu[i,] = scales::rescale(db[i,],from = range(db[i,]),to=c(jcs[i],jcs[i]+1))
}
dbnx = matrix(NA,160,len)
for(i in 1:160){dbnx[i,] = dombi_fn(NULL,menL[i],menC[i],menR[i],menW[i],len)$x}

tikzDevice::tikz(file='C:/Users/NC/Downloads/MDPI_template/MDPI_template/fig5.tex',width=6,height=5)

par(mai=c(1.15, 0.85, 0.45, 0.15))
cols=c("palegreen4","red4","blue4","goldenrod4")
plot(X[1:NROW(menC),3],y=menC[1:NROW(menC)],col="white", lwd=0,xlim=c(min(jcs),max(jcs+1)),ylim = c(min(menC)-0.1,max(menC)), bty="n", main = "",xlab = "emotional$\\_$stability", ylab = "Fuzzy ratings")
set.seed(8542)
for(i in sample(1:160,25)){
  lines(dbnmu[which(ics$religiousness=="1"&ics$university=="1")[i],],dbnx[which(ics$religiousness=="1"&ics$university=="1")[i],], col = cols[1],lty=5,lwd = 2)
  points(y=menC[which(ics$religiousness=="1"&ics$university=="1")[i]],x=max(dbnmu[which(ics$religiousness=="1"&ics$university=="1")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="1"&ics$university=="1")[i],1],x0=jcs[which(ics$religiousness=="1"&ics$university=="1")[i]],ips[which(ics$religiousness=="1"&ics$university=="1")[i],3],x1=jcs[which(ics$religiousness=="1"&ics$university=="1")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="2"&ics$university=="1")[i],],dbnx[which(ics$religiousness=="2"&ics$university=="1")[i],], col = cols[2],lty=5,lwd = 2)
  points(y=menC[which(ics$religiousness=="2"&ics$university=="1")[i]],x=max(dbnmu[which(ics$religiousness=="2"&ics$university=="1")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="2"&ics$university=="1")[i],1],x0=jcs[which(ics$religiousness=="2"&ics$university=="1")[i]],ips[which(ics$religiousness=="2"&ics$university=="1")[i],3],x1=jcs[which(ics$religiousness=="2"&ics$university=="1")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="1"&ics$university=="3")[i],],dbnx[which(ics$religiousness=="1"&ics$university=="3")[i],], col = cols[3],lty=5,lwd = 2)
  points(y=menC[which(ics$religiousness=="1"&ics$university=="3")[i]],x=max(dbnmu[which(ics$religiousness=="1"&ics$university=="3")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="1"&ics$university=="3")[i],1],x0=jcs[which(ics$religiousness=="1"&ics$university=="3")[i]],ips[which(ics$religiousness=="1"&ics$university=="3")[i],3],x1=jcs[which(ics$religiousness=="1"&ics$university=="3")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
  
  lines(dbnmu[which(ics$religiousness=="2"&ics$university=="3")[i],],dbnx[which(ics$religiousness=="2"&ics$university=="3")[i],], col = cols[4],lty=5,lwd = 2)
  points(y=menC[which(ics$religiousness=="2"&ics$university=="3")[i]],x=max(dbnmu[which(ics$religiousness=="2"&ics$university=="3")[i],]),pch=23,bg="antiquewhite4",cex=0.7);
  segments(y0=ips[which(ics$religiousness=="2"&ics$university=="3")[i],1],x0=jcs[which(ics$religiousness=="2"&ics$university=="3")[i]],ips[which(ics$religiousness=="2"&ics$university=="3")[i],3],x1=jcs[which(ics$religiousness=="2"&ics$university=="3")[i]],col = "antiquewhite4",lty = 3, lwd = 2);
}
abline(a=b_est[1],b=b_est[3],col = cols[1], lwd = 2)
abline(a=b_est[1]+b_est[2],b=b_est[3],col = cols[2], lwd = 2)
abline(a=b_est[1]+b_est[4],b=b_est[3],col = cols[3], lwd = 2)
abline(a=b_est[1]+b_est[2]+b_est[4],b=b_est[3],col = cols[4], lwd = 2)

add_legend("bottom",fill = cols,legend = c("religiousness:no,university:no","religiousness:yes,university:no","religiousness:no,university:yes","religiousness:yes,university:yes"),border = FALSE,bty = "n",ncol = 2)

dev.off()
