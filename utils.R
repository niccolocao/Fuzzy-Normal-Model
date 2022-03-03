
dombi_fn = function(x=NULL,l,m,r,w=1,len=101){
  if(is.null(x)){x=seq(from=l,to=r,length.out=len)}
  mux=rep(1e-99,length(x))
  mux[x<=m&x>l]= 1/(1+((m-x[x<=m&x>l])/(x[x<=m&x>l]-l))^w)
  mux[x>m&x<r]= 1/(1+((r-x[x>m&x<r])/(x[x>m&x<r]-m))^-w)
  return(list(x=x,mu=mux))
}


inaccuracy_dombi_fn = function(pars,y,X){
  beta=pars[1:NCOL(X)]
  sigma= pars[NCOL(X)+1]
  lb=y[,1]
  m=y[,2]
  ub=y[,3]
  w=y[,4]
  n=NROW(X)
  mu = X%*%beta
  #print(c(beta,sigma))
  fy = -sum(mapply(function(i){integrate(function(x){dombi_fn(x = x,m = m[i],l = lb[i], r = ub[i], w = w[i])$mu*dnorm(x,mu[i],sd = sigma,log=TRUE)},lb[i],ub[i])$value},1:n))
  return(fy)
}

likelihood_dombi_fn = function(pars,y,X){
  beta=pars[1:NCOL(X)]
  sigma= pars[NCOL(X)+1]
  lb=y[,1]
  m=y[,2]
  ub=y[,3]
  w=y[,4]
  n=NROW(X)
  mu = X%*%beta
  fy = sum(mapply(function(i){log(integrate(function(x){dombi_fn(x = x,m = m[i],l = lb[i], r = ub[i], w = w[i])$mu*dnorm(x,mu[i],sd = sigma)},lb[i]-1e-02,ub[i]+1e-02)$value)},1:n))
  return(fy)
}


getFuzzyNumbers =function(Eta=NULL,Alpha=NULL,Tm=NULL, Time=NULL){
  Alpha=as.matrix(Alpha); Eta=as.matrix(Eta); Time=as.matrix(Time)
  J=NROW(Alpha);M=NROW(Tm);I=NROW(Eta);N=NCOL(Eta)
  Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}
  RTS_median = mapply(function(i) median(RTS[,i]),1:NCOL(RTS)) #Mediana
  ECDF = mapply(function(i) ecdf(RTS[,i]), 1:J) #Distribuzione cumulata empirica; stimatore della FDR
  PY = matrix(NA,I*J,M) #JIxM (J nested within I)
  ETA = kronecker(Eta,matrix(1,J,1)) #JIxN (rep each Eta J times)
  ALPHA = kronecker(matrix(1,I,1),Alpha) #JIxN (rep each Alpha I times)
  for(m in 1:M){
    TM = matrix(1,I*J,1)%*%Tm[m,] 
    DM = matrix(1,I*J,1)%*%Dm[m,]
    PY[,m] = apply( (exp((ETA+ALPHA)*TM)/(matrix(1,J*I,N)+exp((ETA+ALPHA))))^DM, 1, prod)
  }
  xsup=seq(0,1,length.out = M)
  YIJ_stats=matrix(NA,I*J,6); colnames(YIJ_stats) = c("C","S","L","R","W","Kauff")
  YIJ_stats[,1] = mapply(function(i)sum(PY[i,]*xsup),1:(I*J))
  YIJ_stats[,2] = mapply(function(i)1/sum(PY[i,]*(xsup-YIJ_stats[i,1])^2),1:(I*J))
  YIJ_stats[,3:4] = t(mapply(function(i){triangular_fn(mi = YIJ_stats[i,1],si = YIJ_stats[i,2])$par},1:(I*J)))[,1:2]
  YIJ_stats[,5] = as.vector(t(mapply(function(j){mapply(function(i){ECDF[[j]](RTS_median[j])-ECDF[[j]](RTS[i,j])+1},1:I)},1:J))) #TEMPO
  YIJ_stats[,6] = mapply(function(i){kaufmann_index(dombi_fn(l = YIJ_stats[i,3],r = YIJ_stats[i,4],m = YIJ_stats[i,1],w = YIJ_stats[i,5])$mu)},1:(I*J))
  YIJ_stats = data.frame(rep(1:I,each=J),rep(1:J,I),YIJ_stats); colnames(YIJ_stats)[1:2]=c("sbj","itm")
  PY = data.frame(rep(1:I,each=J),rep(1:J,I),PY); colnames(PY)=c("sbj","itm",paste0("m",1:M))
  return(list(PY=PY,Yfuzzy=YIJ_stats))
}


kaufmann_index = function(mux){
  return(2*sum(abs(mux-(mux>=0.5)))/length(mux))
}



triangular_fn = function(mi,si,plotx=FALSE){
  #mi: mode 
  #si: precision
  #Source: Williams, T. M. (1992). Practical use of distributions in network analysis. Journal of the Operational Research Society, 43(3), 265-270.
  lb=0;ub=1 #standard beta function
  mux = (1+mi*si)/((1+mi*si)+(1+si*(1-mi)))
  vx = 1/si; w=3.5
  cx = sqrt(w*vx - 3*(mi-mux))
  bx =  (cx + 3*mi - 3*mux)/2
  ax = mi-bx; 
  ai=ifelse(ax<lb,lb,ax)
  bi=ax+cx; bi=ifelse(bi>ub,ub,bi)
  mi=ax+bx
  xsup=seq(lb,ub,length.out=100)
  fy=extraDistr::dtriang(xsup,ai,bi,mi); fy=fy/max(fy); 
  if(plotx){plot(xsup,fy,bty="n",type="l")}
  return(list(xsup=xsup,mux=fy,par=c(ai,bi,mi)))
}


add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
