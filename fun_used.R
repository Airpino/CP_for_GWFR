library(readr)
library(fda)
#library(cluster)
#library(vows)
#library(clusterSim)
#library(fossil)
#library(geofd)
#library(geoR)
#library(MASS)
library(Rlinsolve)
#library(rworldmap)
#library(mapIT) 
#library(emulator)

index_band3=function(stb, rk, Yb, cp, c)
{
  ST=(colMeans(stb))
  Rk=matrix(0, l, 1)
  #browser()
  for(i in 1:l){
    #  Rk[i]=quantile(rk[cp[,i]==1,i],1-alpha)
    Rk[i]=quantile(rk[,i],1-alpha)
  }
  Cp=rowMeans(cp)
  RK=sqrt(quantile(Rk,1-alpha))
  RK
  YM=colMeans(Yb)
  YM
  
  grafico=T
  if(grafico){
    plf=function(Yf, bse, c, Ylim)
    {
      n=dim(Yf)[1]
      if(is.null(n))
        n=1
      Y=Yf%*%t(bse)
      xylim=list(xlim=c(min(x), max(x)), ylim=Ylim)
      plot(NULL, NULL, type = "l", xlim = xylim$xlim, ylim = xylim$ylim, xlab=c[1], ylab=c[2])
      #browser()
      for(i in 1:n){lines(x, Y[i,])}
      title(c[3])
    }
    
    v=c(min((YM%*%t(bse)-RK*ST), Y[,]%*%t(bse)), max((YM%*%t(bse)+RK*ST), Y[,]%*%t(bse)))
    #browser()
    plf(Y[,], bse, c, v)
    lines(x, YM%*%t(bse), col="blue", lwd = 3)
    lines(x, (YM%*%t(bse)+RK*ST), col="red", lwd = 3)  
    lines(x, (YM%*%t(bse)-RK*ST), col="red", lwd = 3)  
  }
  lu=(YM%*%t(bse)+RK*ST)
  m=YM%*%t(bse)
  lw=(YM%*%t(bse)-RK*ST)
  
  luc=Data2fd(x, t(lu), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(luc-YM)
    return(w.ff)
  }
  b=integrate(fint, lower=0, upper=1, stop.on.error=F)[[1]]
  b
  
  lwc=Data2fd(x, t(lw), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(YM - lwc)
    return(w.ff)
  }
  a=integrate(fint, lower=0, upper=1, stop.on.error=F)[[1]]
  a
  Width=b+a
  Width
  l=nrow(Y)
  B = eval.basis(x, bs)
  inout=matrix(0, l, 1)
  for(i in 1:(l)){
    a=lw[1,]-B%*%(Y[i,])
    b=B%*%(Y[i,])-lu[1,]
    inout[i]=1-((length(a[a>0]!=0))||(length(b[b>0]!=0)))
  }
  Cova=mean(inout)*100
  Cova
  inout=matrix(0, l, 1)
  for(i in 1:(l)){
    fint=function(x){
      #browser()
      B = eval.basis(x, bs)
      a=B%*%(luc-lwc)
      b=B%*%(lwc-Y[i,])
      b[b<0]=0
      c=B%*%(Y[i,]-luc)
      c[c<0]=0
      return(a+(2/alpha)*b+(2/alpha)*c)
    }
    inout[i]=integrate(fint, lower=0, upper=1, stop.on.error=F)[[1]]
  }
  IntScore=mean(inout[])
  IntScore
  #browser()
  ### covalpha funzionale
  inout=matrix(0, l, length(x))
  for(i in 1:(l)){
    a=lw[1,]-B%*%(Y[i,])
    b=B%*%(Y[i,])-lu[1,]
    inout[i, ((1:length(x))[(a<0)&(b<0)]) ]=1
  }
  Covat=colMeans(inout)*100
  #####   indici che mi servono
  return(list(c(Cova, Width, IntScore), Covat, YM%*%t(bse), RK, ST))
  #
}

index_band=function(stb, rk, Yb, cp, c)
{
  ST=colMeans(stb)
  Rk=matrix(0, l, 1)
  for(i in 1:l){
    #  Rk[i]=quantile(rk[cp[,i]==1,i],1-alpha)
    Rk[i]=quantile(rk[,i],1-alpha)
  }
  Cp=rowMeans(cp)
  RK=sqrt(quantile(Rk,1-alpha))
  RK
  YM=colMeans(Yb)
  YM
  
  grafico=T
  if(grafico){
    plf=function(Ys, bse, c, Ylim)
    {
      n=dim(Ys)[1]
      if(is.null(n))
        n=1
      Y=Ys%*%t(bse)
      xylim=list(xlim=c(min(x), max(x)), ylim=Ylim)
      plot(NULL, NULL, type = "l", xlim = xylim$xlim, ylim = xylim$ylim, xlab=c[1], ylab=c[2])
      for(i in 1:n){lines(x, Y[i,])}
      title(c[3])
    }
    
    v=c(min((YM%*%t(bse)-RK*ST), Y[,]), max((YM%*%t(bse)+RK*ST), Y[,]))
    
    plf(Y[,], bse, c, v)
    lines(x, YM%*%t(bse), col="blue", lwd = 3)
    lines(x, (YM%*%t(bse)+RK*ST), col="red", lwd = 3)  
    lines(x, (YM%*%t(bse)-RK*ST), col="red", lwd = 3)  
  }
  
  lu=(YM%*%t(bse)+RK*ST)
  m=YM%*%t(bse)
  lw=(YM%*%t(bse)-RK*ST)
  
  luc=Data2fd(x, t(lu), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(luc-YM)
    return(w.ff)
  }
  b=integrate(fint, lower=0, upper=1)[[1]]
  b
  
  lwc=Data2fd(x, t(lw), bs)$coef
  fint=function(x){
    B = eval.basis(x, bs)
    w.ff=B%*%(YM - lwc)
    return(w.ff)
  }
  a=integrate(fint, lower=0, upper=1)[[1]]
  a
  Width=b+a
  Width
  
  B = eval.basis(x, bs)
  inout=matrix(0, 2*l, 1)
  for(i in trial){#1:(2*l)){
    a=lw[1,]-B%*%(Y[i,])
    b=B%*%(Y[i,])-lu[1,]
    inout[i]=1-((length(a[a>0]!=0))||(length(b[b>0]!=0)))
  }
  Cova=mean(inout[trial])*100
  Cova
  
  inout=matrix(0, 2*l, 1)
  for(i in trial){#1:(2*l)){
    fint=function(x){
      B = eval.basis(x, bs)
      a=B%*%(luc-lwc)
      b=B%*%(lwc-Y[i,])
      b[b<0]=0
      c=B%*%(Y[i,]-luc)
      c[c<0]=0
      return(a+(2/alpha)*b+(2/alpha)*c)
    }
    inout[i]=integrate(fint, lower=0, upper=1)[[1]]
  }
  IntScore=mean(inout[trial])
  IntScore
  
  #####   indici che mi servono
  return(c(Cova, Width, IntScore))
  #
}
##compute the numerator and denominator matrices B and A respectively
s.computAB<-function(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=rep(1, ngrid), delta=delta, s=s.dip){
  
  nbasis=ncol(B.grid.weight)
  n.curve=ncol(coef)
  A=B1=B2=matrix(0, nbasis, nbasis)
  
  
  
  for(i in 1:(n.curve-1)){
    for(j in (i+1):n.curve){
      V=arrayVp[,,i]+arrayVp[,,j]
      varf=diag(B.grid%*%V%*%t(B.grid))
      #### here, the use of spatially dipendens
      ebeta=s[i,j]*coef[,i]-s[j,i]*coef[,j]
      fhat=as.vector(B.grid%*%ebeta)
      
      ### calculate the squared mean
      temp=B.grid.weight*((fhat^2+varf)*sqrt(weight)*delta)
      A1=apply(temp, 2, sum)
      A=A+(A1%*%t(A1))
      
      ### calculate the variance
      M=t(B.grid)%*%(B.grid*(weight*delta))
      temp2=diag(B.grid%*%V%*%M%*%V%*%t(B.grid))
      B1=B1+(t(B.grid.weight)%*%(B.grid.weight*temp2*delta))
      
      tempv=t(B.grid.weight)%*%(B.grid*(fhat*sqrt(weight)*delta))
      B2=B2+(tempv%*%V%*%t(tempv))
    }
  }
  
  
  Bv=2*B1+4*B2
  
  return(list(A=A, Bv=Bv))
}


########################################################################################
#### coef: an nbasisxn spline coefficient matrix 
#### arrayVp: nbasisXnbasisXn array,  Bayesian variance of the spline coefficient
#### B.grid: ngridXnbasis spline design matrix for the smoothed curves
#### B.grid.weight: spline design matrix for the weight function
#### t.grid: the equally spaced grid points 
########################################################################################

s.weight.minCV=function(coef=coef, arrayVp=arrayVp, 
                        B.grid=B.grid, B.grid.weight=B.grid.weight, t.grid=t.grid,
                        niter=100, tol=1e-4, s.dip=s.dip){
  
  ngrid=nrow(B.grid.weight) ###number of grid points
  delta=(max(t.grid)-min(t.grid))/ngrid
  fit=s.computAB(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=rep(1, ngrid),  delta=delta, s=s.dip)
  svdB=svd(fit$Bv)
  Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v) ##B^{-1/2}
  svdBA=svd(Bh%*%fit$A%*%Bh)
  qnew=Bh%*%svdBA$u[,1]
  CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
  weight=as.vector((B.grid.weight%*%qnew)^2)
  weight=weight/sum(weight*delta)
  
  #    par(mfrow=c(3,3))
  diff=iter=1
  while(diff>tol & iter<niter){
    cat("iteration=", iter,
        "difference=", diff,
        "CVar^2=", CVar2, "\n")
    
    qold=qnew
    fit=s.computAB(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=weight, delta=delta, s=s.dip)
    svdB=svd(fit$Bv)
    Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v)
    svdBA=svd(Bh%*%fit$A%*%Bh)
    qnew=Bh%*%svdBA$u[,1]; 
    qnew=(qnew+qold)/2
    weight=as.vector((B.grid.weight%*%qnew)^2);  
    weight=weight/sum(weight*delta)
    diff=mean(abs(B.grid.weight%*%(qold-qnew)))
    iter=iter+1
    CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
    #        plot(t.grid, weight, pch=20, col='red', ylab='Weight', xlab='t')
  }        
  cat("iteration=", iter,
      "difference=", diff,
      "CVar^2=", CVar2, "\n")
  
  return(list(weight=weight, q=qnew, B.grid.weight=B.grid.weight, CVar2=CVar2, iter=iter))
}

########################################################################################
####  
#### 
#### 
#### 
#### 
########################################################################################
Gh=function(g, coord, Eu.d)
{
  emp.trace.vari <- trace.variog(coords=coord,
                                 L2norm=g, bin=F,
                                 uvec="default",breaks="default")
  nugget.fix=NULL
  max.dist.variogram=NULL
  sigma2.0<-quantile(emp.trace.vari$v,0.75, na.rm=T)
  phi.0<-quantile(emp.trace.vari$u,0.75, na.rm=T)
  if(is.null(nugget.fix)){
    fix.nugget<-FALSE
    nugget<-0
  }else{
    fix.nugget=TRUE
    nugget<-nugget.fix
  }
  if (is.null(max.dist.variogram))
    max.dist.variogram<-max(emp.trace.vari$u, na.rm=T)   
  models=fit.tracevariog(emp.trace.vari, models=c("exponential","spherical", "matern", "circular", "cubic", "powered.exponential", "pure.nugget","wave","gaussian", "gneiting", "cubic", "cauchy"), ##cubic is ok ##circular is better
                         sigma2.0,phi.0,
                         max.dist=max.dist.variogram,fix.nugget=fix.nugget,
                         nugget=nugget)
  
  sigm2 <- models$best$cov.pars[1]
  Sig <-sigm2 - cov.spatial(Eu.d,cov.model=models$best$cov.model,
                            cov.pars=models$best$cov.pars,
                            kappa=models$best$kappa)
  return(Sig)
}



Beta_plot=function(L, x1, x2, bse1, bse2, phi=30, teta=45, flag=F)
{
  l=dim(bse1)[1]
  Z=list()
  h=length(L)
  for(i in 1:h)
  {
    Z[[i]]=bse1%*%L[[i]]%*%t(bse2)
    persp(x1, x2, Z[[i]],
          #theta = 45, phi = 30, 
          phi=phi,theta=teta,d=10,#col="blue",ticktype="detaile", 
          scale = flag #sltheta = -120, shade = 0.75, border = NA, box = FALSE)
          #, xlab = "t", ylab = "tau", zlab=expression(beta(t,tau))
          , xlab = " ", ylab = " ", zlab=" "
    )
  }
  return("Beta_inter"=Z)
}



##################
###   cross validation per il valore della h
##################
cvh=function(ud, Y, X,
             #Wr=diag(dim(Y)[1]),
             Wr=rep(1,nrow(Y)), #ora Wr è un vettore
             J_f1, J_f2=J_f1, flag)
{
  h=1:20
  cv=h
  m=nrow(Y)
  ######## ----------------------------
  H=length(X)
  M=ncol(X[[1]])
  Xbar=matrix(0, nrow = nrow(Y), ncol =M*H )
  
  for(i in 1:H){  
    Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f1)
  }
  d<-crossprod(Xbar, Y)
  ############## ----------------------
  for(j in h)
  {
    
    W=exp(-ud/h[j])*Wr
    #browser()
    # res=lmr(Y, X, W,#ora W è un vettore 
    #         J_f1, flag)
    res=lmr_3(Y, X, Xbar, W,#ora W è un vettore 
              J_f1, flag)
    #residual=matrix(0, m, 1)
    yb=Y-res$HatMatrix%*%Y
    # cv[j]=sum(diag(yb%*%J_f2%*%t(yb)))
    cv[j]=sum((yb%*%J_f2)*yb)
  }
  
  W<-exp(-ud/which.min(cv))
  return(W)
}

cvh_2=function(ud, Y, X, Xbar,d,
               #Wr=diag(dim(Y)[1]),
               Wr=rep(1,nrow(Y)), #ora Wr è un vettore
               J_f1, J_f2=J_f1, flag)
{
  h=1:20
  cv=h
  m=nrow(Y)
  
  for(j in h)
  {
    
    W=exp(-ud/h[j])*Wr
    
    # res=lmr_2(Y, X, Xbar, d, W,#ora W è un vettore 
    #           J_f1, flag)
    res=lmr_3(Y, X, Xbar,  W,#ora W è un vettore 
              flag)
    #residual=matrix(0, m, 1)
    yb=Y-res$HatMatrix%*%Y
    cv[j]=sum((yb%*%J_f2)*yb)
  }
  
  W<-exp(-ud/which.min(cv))
  return(W)
}


##################
###   cross validation per i pesi rispetto alla posizione s
##################
cvW=function(Ds, Y, X, 
             #Wr=diag(dim(Y)[1]),
             Wr=rep(1,nrow(Y)), #ora Wr è un vettore
             J_f1, J_f2=J_f1, flag)
{
  
  m=dim(Ds)[1]
  cv=1:m
  mi=Inf
  
  ######## ----------------------------
  H=length(X)
  M=ncol(X[[1]])
  Xbar=matrix(0, nrow = nrow(Y), ncol =M*H )
  
  for(i in 1:H){  
    Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f1)
  }
  d<-crossprod(Xbar, Y)
  ############## ----------------------
  
  
  for(j in 1:m)
  {
    ud=Ds[j,]
    #W=cvh(ud, Y, X, Wr, J_f1, J_f2, flag)
    W=cvh_2(ud, Y, X, Xbar,d,
            #Wr=diag(dim(Y)[1]),
            Wr, #ora Wr è un vettore
            J_f1, J_f2=J_f1, flag)
    
    #browser()
    Wbf=W*Wr
    #residual=matrix(0, m, 1)
    #res=lmr(Y, X, Wbf, J_f1, flag)
    
    # res=lmr_2(Y, X, Xbar, d, Wbf,#ora W è un vettore 
    #           J_f1, flag)
    res=lmr_3(Y, X, Xbar,  Wbf,#ora W è un vettore 
              flag)
    yb=Y-res$HatMatrix%*%Y
    # cv[j]=sum(diag(yb%*%J_f2%*%t(yb)))#vedi sopra
    cv[j]<-sum((yb%*%J_f2)*yb)
    if(cv[j]<mi) {
      Wb=W
      mi=cv[j]
    }
  }
  return(Wb)
}


lmr=function(Y, X, #W=diag(dim(Y)[1]), 
             W=rep(1,nrow(Y)),#ora W è un vettore
             J_f, flag=0)
{
  # see 
  #https://stackoverflow.com/questions/20562177/get-hat-matrix-from-qr-decomposition-for-weighted-least-square-regression
  Y=as.matrix(Y)
  #coefficienti della Y: ci sono 100 unità, funzioni approssimanti 4
  N=nrow(Y)
  K=ncol(Y)
  #matrici dei coefficienti delle X_g, con g=1,...,H
  H=length(X)
  M=ncol(X[[1]])
  # vbar=list()
  # for(i in 1:H){vbar[[i]]=X[[i]]}
  mh=M*H
  q=0
  Xbar=matrix(0, nrow = N, ncol = mh)
  
  for(i in 1:H){  
    Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f)
  }
  if(flag){
    unos=matrix(1, nrow=N, ncol=1)
    Xbar=cbind(unos, Xbar)
    q=1
  }
  #A=t(Xbar)%*%W%*%Xbar
  #A=(mmult(t(Xbar),W))%*%Xbar
  
  A<-tryCatch(
    {
      solve((mmult(t(Xbar),W))%*%Xbar)
    },
    error=function(cond) {
      ginv((mmult(t(Xbar),W))%*%Xbar)
    })
  
  #d=t(Xbar)%*%Y
  d<-crossprod(Xbar, Y)
  B=list()
  unos=matrix(0, nrow=K, ncol=1)
  b=matrix(0, ncol=K, nrow = (M+q))
  if(flag){    
    for(j in 1:H){
      B[[j]]=b[-q, ]
    }
  }else{
    for(j in 1:H){
      B[[j]]=b
    }
  }
  #  for(i in 1:dim(d)[2]){b[,i]=lsolve.cgs(A, d[,i], xinit = NA, reltol = 1e-15, maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), adjsym = TRUE, verbose = TRUE)$x  }
  
  # b=lsolve.cgs(A, d, xinit = NA, reltol = 1e-15, 
  #              maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), 
  #              adjsym = TRUE, verbose = TRUE)$x
  #### ---------------------------
  # browser()
  #b=solve(A,d) #PERCHE' non VA bene?
  b=A%*%d
  ### ---------------------------
  unos=b[q, ]
  for(i in 1:H){
    B[[i]]=b[(1+q+(M*(i-1))):(M*i+q),]
  }
  
  # A<-tryCatch(
  #  {
  #    solve(A)
  #  },
  #  error=function(cond) {
  #    ginv(A)
  #  })
  H=mmult(Xbar%*%A%*%t(Xbar),W)
  return(list("beta"=B, "HatMatrix"=H, "beta_0"=unos))
}

lmr_1.1=function(Y, X, #W=diag(dim(Y)[1]), 
                 W=rep(1,nrow(Y)),#ora W è un vettore
                 J_f, flag=0){
  # see 
  #https://stackoverflow.com/questions/20562177/get-hat-matrix-from-qr-decomposition-for-weighted-least-square-regression
  Y=as.matrix(Y)
  #coefficienti della Y: ci sono 100 unità, funzioni approssimanti 4
  N=nrow(Y)
  K=ncol(Y)
  #matrici dei coefficienti delle X_g, con g=1,...,H
  H=length(X)
  M=ncol(X[[1]])
  # vbar=list()
  # for(i in 1:H){vbar[[i]]=X[[i]]}
  mh=M*H
  q=0
  Xbar=matrix(0, nrow = N, ncol = mh)
  
  for(i in 1:H){  
    Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f)
  }
  if(flag){
    unos=matrix(1, nrow=N, ncol=1)
    Xbar=cbind(unos, Xbar)
    q=1
  }
  #A=t(Xbar)%*%W%*%Xbar
  #A=(mmult(t(Xbar),W))%*%Xbar
  
  A<-tryCatch(
    {
      solve((mmult(t(Xbar),W))%*%Xbar)
    },
    error=function(cond) {
      ginv((mmult(t(Xbar),W))%*%Xbar)
    })
  
  #d=t(Xbar)%*%Y
  # d<-crossprod(Xbar, Y)
  # B=list()
  # unos=matrix(0, nrow=K, ncol=1)
  # b=matrix(0, ncol=K, nrow = (M+q))
  # if(flag){    
  #   for(j in 1:H){
  #     B[[j]]=b[-q, ]
  #   }
  # }else{
  #   for(j in 1:H){
  #     B[[j]]=b
  #   }
  # }
  #  for(i in 1:dim(d)[2]){b[,i]=lsolve.cgs(A, d[,i], xinit = NA, reltol = 1e-15, maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), adjsym = TRUE, verbose = TRUE)$x  }
  
  # b=lsolve.cgs(A, d, xinit = NA, reltol = 1e-15, 
  #              maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), 
  #              adjsym = TRUE, verbose = TRUE)$x
  #### ---------------------------
  # browser()
  #b=solve(A,d) #PERCHE' non VA bene?
  # b=A%*%d
  ### ---------------------------
  # unos=b[q, ]
  # for(i in 1:H){
  #   B[[i]]=b[(1+q+(M*(i-1))):(M*i+q),]
  # }
  
  # A<-tryCatch(
  #  {
  #    solve(A)
  #  },
  #  error=function(cond) {
  #    ginv(A)
  #  })
  H=mmult(Xbar%*%A%*%t(Xbar),W)
  return(list("HatMatrix"=H))
}

lmr_2=function(Y, X, Xbar, d, #W=diag(dim(Y)[1]), 
               W=rep(1,nrow(Y)),#ora W è un vettore
               J_f, flag=0)
{
  Y=as.matrix(Y)
  #coefficienti della Y: ci sono 100 unità, funzioni approssimanti 4
  N=nrow(Y)
  K=ncol(Y)
  #matrici dei coefficienti delle X_g, con g=1,...,H
  H=length(X)
  M=ncol(X[[1]])
  
  mh=M*H
  q=0
  
  if(flag){
    unos=matrix(1, nrow=N, ncol=1)
    
    Xbar_OLD=Xbar
    Xbar=cbind(unos, Xbar)
    d<-rbind(crossprod(unos,Y),d)
    q=1
  }
  
  # X_bb<-diag(sqrt(W))%*%Xbar
  # #mmult(Xbar,sqrt(W))
  #   Y_bb<-diag(sqrt(W))%*%Y
  #   #mmult(Y,sqrt(W))
  #   browser()
  #A=t(Xbar)%*%W%*%Xbar
  #A=(mmult(t(Xbar),W))%*%Xbar
  
  A<-tryCatch(
    {
      solve((mmult(t(Xbar),W))%*%Xbar)
    },
    error=function(cond) {
      ginv((mmult(t(Xbar),W))%*%Xbar)
    })
  
  #d=t(Xbar)%*%Y
  d<-crossprod(Xbar, Y)
  
  B=list()
  unos=matrix(0, nrow=K, ncol=1)
  b=matrix(0, ncol=K, nrow = (M+q))
  if(flag){    
    for(j in 1:H){
      B[[j]]=b[-q, ]
    }
  }else{
    for(j in 1:H){
      B[[j]]=b
    }
  }
  
  b=A%*%d
  ### ---------------------------
  unos=b[q, ]
  for(i in 1:H){
    B[[i]]=b[(1+q+(M*(i-1))):(M*i+q),]
  }
  
  # A<-tryCatch(
  #  {
  #    solve(A)
  #  },
  #  error=function(cond) {
  #    ginv(A)
  #  })
  H=mmult(Xbar%*%A%*%t(Xbar),W)
  # H1<-mmult(Xbar%*%solve((mmult(t(Xbar),W))%*%Xbar)%*%t(Xbar),W)
  # H2<-diag(1/sqrt(W))%*%X_bb%*%solve(crossprod(X_bb))%*%t(X_bb)%*%diag(sqrt(W))
  # H3<-Xbar%*%solve((mmult(t(Xbar),W))%*%Xbar)%*%t(Xbar)%*%diag(W)
  # H4<-diag(1/sqrt(W))%*%tcrossprod(qr.Q(qr(X_bb)))%*%diag(sqrt(W))
  # browser()
  return(list("beta"=B, "HatMatrix"=H, "beta_0"=unos))
}

lmr_3=function(Y, X, Xbar, #W=diag(dim(Y)[1]), 
               W=rep(1,nrow(Y)),#ora W è un vettore
               flag=0)
{
  Y=as.matrix(Y)
  #coefficienti della Y: ci sono 100 unità, funzioni approssimanti 4
  N=nrow(Y)
  K=ncol(Y)
  #matrici dei coefficienti delle X_g, con g=1,...,H
  H=length(X)
  M=ncol(X[[1]])
  
  mh=M*H
  q=0
  
  if(flag){
    unos=matrix(1, nrow=N, ncol=1)
    
    Xbar_OLD=Xbar
    Xbar=cbind(unos, Xbar)
    q=1
  }
  
  H<-tryCatch(
    {
      mmult(Xbar%*%solve((mmult(t(Xbar),W))%*%Xbar)%*%t(Xbar),W)
      # solve((mmult(t(Xbar),W))%*%Xbar)
    },
    error=function(cond) {
      mmult(Xbar%*%ginv((mmult(t(Xbar),W))%*%Xbar)%*%t(Xbar),W)
      #ginv((mmult(t(Xbar),W))%*%Xbar)
    })
  #H=mmult(Xbar%*%A%*%t(Xbar),W)
  
  return(list("HatMatrix"=H))
}
## ----CPP funs ----------------
library( Rcpp )

#  Source code for our function
cppFunction(
  'NumericMatrix mmult( NumericMatrix m , NumericVector v , bool byrow = true ){
  if( byrow );
    if( ! m.nrow() == v.size() ) stop("Non-conformable arrays") ;
  if( ! byrow );
    if( ! m.ncol() == v.size() ) stop("Non-conformable arrays") ;

  NumericMatrix out(m) ;

  if( byrow ){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if( ! byrow ){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
}'
)

cppFunction(
  'NumericMatrix vmmult( NumericVector v, NumericMatrix m  ){
    if( ! m.ncol() == v.size() ) stop("Non-conformable arrays") ;

  NumericMatrix out(m) ;

    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  
  return out ;
}'
)




#  Source code for our function
cppFunction(
  'double myf_c( NumericMatrix m , NumericMatrix n ){

  double out=0;

  for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out = out + m(i,j) * n(j,i);
      }
    }
  
  return out ;
}'
)
#### ----------------------------

############### old functions
##compute the numerator and denomator matrices B and A respectively
s.computAB_OLD<-function(coef=coef, arrayVp=arrayVp,
                         B.grid=B.grid, B.grid.weight=B.grid.weight, 
                         weight=rep(1, ngrid), delta=delta, s=s.dip){
  
  nbasis=ncol(B.grid.weight)
  n.curve=ncol(coef)
  A=B1=B2=matrix(0, nbasis, nbasis)
  
  
  
  for(i in 1:(n.curve-1)){
    for(j in (i+1):n.curve){
      V=arrayVp[,,i]+arrayVp[,,j]
      varf=diag(B.grid%*%V%*%t(B.grid))
      #### here, the use of spatially dipendens
      ebeta=s[i,j]*coef[,i]-s[j,i]*coef[,j]
      fhat=as.vector(B.grid%*%ebeta)
      
      ### calculate the squared mean
      temp=B.grid.weight*((fhat^2+varf)*sqrt(weight)*delta)
      A1=apply(temp, 2, sum)
      A=A+(A1%*%t(A1))
      
      ### calculate the variance
      M=t(B.grid)%*%(B.grid*(weight*delta))
      temp2=diag(B.grid%*%V%*%M%*%V%*%t(B.grid))
      B1=B1+(t(B.grid.weight)%*%(B.grid.weight*temp2*delta))
      
      tempv=t(B.grid.weight)%*%(B.grid*(fhat*sqrt(weight)*delta))
      B2=B2+(tempv%*%V%*%t(tempv))
    }
  }
  
  
  Bv=2*B1+4*B2
  
  return(list(A=A, Bv=Bv))
}


########################################################################################
#### coef: an nbasisxn spline coefficient matrix 
#### arrayVp: nbasisXnbasisXn array,  Bayesian variance of the spline coefficient
#### B.grid: ngridXnbasis spline design matrix for the smoothed curves
#### B.grid.weight: spline design matrix for the weight function
#### t.grid: the equally spaced grid points 
########################################################################################

s.weight.minCV_OLD=function(coef=coef, arrayVp=arrayVp,
                            B.grid=B.grid, B.grid.weight=B.grid.weight, 
                            t.grid=t.grid, niter=100, tol=1e-4, s.dip=s.dip){
  
  ngrid=nrow(B.grid.weight) ###number of grid points
  delta=(max(t.grid)-min(t.grid))/ngrid
  fit=s.computAB_OLD(coef=coef, arrayVp=arrayVp, 
                     B.grid=B.grid, B.grid.weight=B.grid.weight, 
                     weight=rep(1, ngrid),  delta=delta, s=s.dip)
  svdB=svd(fit$Bv)
  Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v) ##B^{-1/2}
  svdBA=svd(Bh%*%fit$A%*%Bh)
  qnew=Bh%*%svdBA$u[,1]
  CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
  weight=as.vector((B.grid.weight%*%qnew)^2)
  weight=weight/sum(weight*delta)
  
  #    par(mfrow=c(3,3))
  diff=iter=1
  while(diff>tol & iter<niter){
    cat("iteration=", iter,
        "difference=", diff,
        "CVar^2=", CVar2, "\n")
    
    qold=qnew
    fit=s.computAB_OLD(coef=coef, arrayVp=arrayVp, B.grid=B.grid, B.grid.weight=B.grid.weight, weight=weight, delta=delta, s=s.dip)
    svdB=svd(fit$Bv)
    Bh=svdB$u%*%diag(svdB$d^(-1/2))%*%t(svdB$v)
    svdBA=svd(Bh%*%fit$A%*%Bh)
    qnew=Bh%*%svdBA$u[,1]; 
    qnew=(qnew+qold)/2
    weight=as.vector((B.grid.weight%*%qnew)^2);  
    weight=weight/sum(weight*delta)
    diff=mean(abs(B.grid.weight%*%(qold-qnew)))
    iter=iter+1
    CVar2=t(qnew)%*%fit$Bv%*%qnew/(t(qnew)%*%fit$A%*%qnew)
    #        plot(t.grid, weight, pch=20, col='red', ylab='Weight', xlab='t')
  }        
  cat("iteration=", iter,
      "difference=", diff,
      "CVar^2=", CVar2, "\n")
  
  return(list(weight=weight, q=qnew, B.grid.weight=B.grid.weight, CVar2=CVar2, iter=iter))
}

########################################################################################
####  
#### 
#### 
#### 
#### 
########################################################################################
Gh_OLD=function(g, coord, Eu.d)
{
  emp.trace.vari <- trace.variog(coords=coord,
                                 L2norm=g, bin=F,
                                 uvec="default",breaks="default")
  nugget.fix=NULL
  max.dist.variogram=NULL
  sigma2.0<-quantile(emp.trace.vari$v,0.75, na.rm=T)
  phi.0<-quantile(emp.trace.vari$u,0.75, na.rm=T)
  if(is.null(nugget.fix)){
    fix.nugget<-FALSE
    nugget<-0
  }else{
    fix.nugget=TRUE
    nugget<-nugget.fix
  }
  if (is.null(max.dist.variogram))
    max.dist.variogram<-max(emp.trace.vari$u, na.rm=T)   
  models=fit.tracevariog(emp.trace.vari, models=c("exponential","spherical", "matern", "circular", "cubic", "powered.exponential", "pure.nugget","wave","gaussian", "gneiting", "cubic", "cauchy"), ##cubic is ok ##circular is better
                         sigma2.0,phi.0,
                         max.dist=max.dist.variogram,fix.nugget=fix.nugget,
                         nugget=nugget)
  
  sigm2 <- models$best$cov.pars[1]
  Sig <-sigm2 - cov.spatial(Eu.d,cov.model=models$best$cov.model,
                            cov.pars=models$best$cov.pars,
                            kappa=models$best$kappa)
  return(Sig)
}



Beta_plot_OLD=function(L, x1, x2, bse1, bse2, phi=30, teta=45, flag=F)
{
  l=dim(bse1)[1]
  Z=list()
  h=length(L)
  for(i in 1:h)
  {
    Z[[i]]=bse1%*%L[[i]]%*%t(bse2)
    persp(x1, x2, Z[[i]],
          #theta = 45, phi = 30, 
          phi=phi,theta=teta,d=10,#col="blue",ticktype="detaile", 
          scale = flag #sltheta = -120, shade = 0.75, border = NA, box = FALSE)
          #, xlab = "t", ylab = "tau", zlab=expression(beta(t,tau))
          , xlab = " ", ylab = " ", zlab=" "
    )
  }
  return("Beta_inter"=Z)
}



##################
###   cross validation per il valore della h
##################
cvh_OLD=function(ud, Y, X, Wr=diag(dim(Y)[1]), J_f1, J_f2=J_f1, flag)
{
  h=1:20
  cv=h
  m=dim(Y)[1]
  for(j in h)
  {
    W=diag(exp(-ud/h[j]))%*%Wr
    res=lmr_OLD(Y, X, W, J_f1, flag)
    
    residual=matrix(0, m, 1)
    yb=Y-res$HatMatrix%*%Y
    for(i in 1:m){
      residual[i]=t(yb[i,])%*%J_f2%*%(yb[i,])
    }
    cv[j]=sum(residual)
  }
  W=diag(exp(-ud/which.min(cv)))
  return(W)
}

##################
###   cross validation per i pesi rispetto alla posizione s
##################
cvW_OLD=function(Ds, Y, X, Wr=diag(dim(Y)[1]), J_f1, J_f2=J_f1, flag)
{
  m=dim(Ds)[1]
  cv=1:m
  ud=Ds[1,]
  W=cvh_OLD(ud, Y, X, Wr, J_f1, J_f2, flag)
  Wbf=W%*%Wr
  residual=matrix(0, m, 1)
  res=lmr_OLD(Y, X, Wbf, J_f1, flag)
  yb=Y-res$HatMatrix%*%Y
  for(i in 1:m){
    residual[i]=t(yb[i,])%*%J_f2%*%(yb[i,])
  }
  cv[1]=sum(residual)
  mi=cv[1]
  Wb=W
  for(j in 2:m)
  {
    ud=Ds[j,]
    W=cvh_OLD(ud, Y, X, Wr, J_f1, J_f2, flag)
    Wbf=W%*%Wr
    residual=matrix(0, m, 1)
    res=lmr_OLD(Y, X, Wbf, J_f1, flag)
    yb=Y-res$HatMatrix%*%Y
    for(i in 1:m){
      residual[i]=t(yb[i,])%*%J_f2%*%(yb[i,])
    }
    cv[j]=sum(residual)
    if(cv[j]<mi) {
      Wb=W
      mi=cv[j]
    }
  }
  return(Wb)
}

lmr_OLD=function(Y, X, W=diag(dim(Y)[1]), J_f, flag=0)
{
  Y=as.matrix(Y)
  #coefficienti della Y: ci sono 100 unità, funzioni approssimanti 4
  N=dim(Y)[1]
  K=dim(Y)[2]
  #matrici dei coefficienti delle X_g, con g=1,...,H
  H=length(X)
  M=dim(X[[1]])[2]
  vbar=list()
  for(i in 1:H){vbar[[i]]=X[[i]]}
  mh=M*H
  q=0
  Xbar=matrix(0, nrow = N, ncol = mh)
  for(i in 1:H){  Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f)}
  if(flag){
    unos=matrix(1, nrow=N, ncol=1)
    Xbar=cbind(unos, Xbar)
    q=1
  }
  A=t(Xbar)%*%W%*%Xbar
  d=t(Xbar)%*%Y
  B=list()
  unos=matrix(0, nrow=K, ncol=1)
  b=matrix(0, ncol=K, nrow = (M+q))
  if(flag){    for(j in 1:H){B[[j]]=b[-q, ]}
  }else{
    for(j in 1:H){B[[j]]=b}
  }
  #  for(i in 1:dim(d)[2]){b[,i]=lsolve.cgs(A, d[,i], xinit = NA, reltol = 1e-15, maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), adjsym = TRUE, verbose = TRUE)$x  }
  b=lsolve.cgs(A, d, xinit = NA, reltol = 1e-15, maxiter = dim(A)[1],preconditioner = diag(dim(A)[1]), adjsym = TRUE, verbose = TRUE)$x
  
  #b=solve(A,d)
  unos=b[q, ]
  for(i in 1:H){B[[i]]=b[(1+q+(M*(i-1))):(M*i+q),]}
  A=ginv(A, tol = 1e-330)
  H=Xbar%*%A%*%t(Xbar)%*%W
  return(list("beta"=B, "HatMatrix"=H, "beta_0"=unos))
}

##### NEW
##################
###   cross validation per i pesi rispetto alla posizione s
##################

cvW_new2<-function(Ds, #matrice distanze geo
                   Y, 
                   X, 
                   Wr=rep(1,nrow(Y)), #ora Wr è un vettore
                   J_f1, 
                   J_f2=J_f1, 
                   flag)#mi dice se aggiungere o no il termine noto nel modello
{
  
  m=dim(Ds)[1] #righe della matrice di distanze geo
  cv=1:m
  mi=Inf
  
  ######## ----------------------------
  H=length(X)
  M=ncol(X[[1]])
  Xbar=matrix(0, nrow = nrow(Y), ncol =M*H )
  
  for(i in 1:H){  
    Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f1)
  }
  #d<-crossprod(Xbar, Y)
  ############## ------Metto il vettore unitario se flag=1 ----------------
  if(flag){
    unos=matrix(1, nrow=nrow(Y), ncol=1)
    
    Xbar_OLD=Xbar
    Xbar=cbind(unos, Xbar)
  }
  ##############-------
  for(j in 1:m)
  {
    ud=Ds[j,]
    
    {
      J_f2=J_f1
      h2=1:20
      cv2=h2
      
      for(j2 in h2)
      {
        
        W=exp(-ud/h2[j2])*Wr
        
        XW=mmult(t(Xbar),sqrt(W))
        XWX=tcrossprod(XW)
        yb=Y-mmult(emulator::quad.form.inv(XWX,t(Xbar)),W)%*%Y
        #yb=Y-Xbar%*%solve(t(Xbar)%*%W%*$Xbar)%*%Xbar%*%W%*%Y
        
        cv2[j2]=sum((yb%*%J_f2)*yb)
      }
      
      W<-exp(-ud/which.min(cv2))
      #return(W)
    }
    Wbf=W*Wr
    XW=mmult(t(Xbar),sqrt(Wbf))
    XWX=tcrossprod(XW)
    yb=Y-mmult(emulator::quad.form.inv(XWX,t(Xbar)),Wbf)%*%Y    
    
    cv[j]<-sum((yb%*%J_f2)*yb)
    if(cv[j]<mi) {
      Wb=Wbf
      mi=cv[j]
    }
  }
  return(Wb)
}

cppFunction(depends = "RcppArmadillo", ' 
NumericVector c_f3(arma::mat Xbar,
                arma::mat Y,
                Rcpp::NumericVector W) {
  
arma::mat Xptilde;
 arma::mat Yptilde; 
NumericVector yb;
   Xptilde=mmult(Xbar.t(),sqrt(W));
         Yptilde=mmult(Y.t(),sqrt(W));
      arma::mat tmp1;
      arma::mat tmp2;
            tmp1 =Xptilde*Xptilde.t();
      tmp2 = (Xptilde*Yptilde.t());
      yb = Y-Xbar*arma::solve(tmp1, tmp2);
      
  return yb;
}
', includes ={'arma::mat mmult( arma::mat m , NumericVector v , bool byrow = true ){
  
  arma::mat out(m) ;
  int i,j;
  if( byrow ){
    for (j = 0; j < m.n_cols; j++) {
      for (i = 0; i < m.n_rows; i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if( ! byrow ){
    for (i = 0; i < m.n_rows; i++) {
      for (j = 0; j < m.n_cols; j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
};
'})


cvW_new3<-function(Ds, #distanze geografiche
                   Xbar,
                   Y, #Ys
                   Wr,
                   J_f2=J_f1){
  m=dim(Ds)[1] #righe della matrice di distanze geo
  cv=1:m
  mi=Inf 
  ##############-------
  for(j in 1:m)
  {
    ud=Ds[j,]
    
    {
      
      h2=1:20
      cv2=h2
      
      for(j2 in h2)
      {
        
        W=exp(-ud/h2[j2])*Wr
        Xptilde<-mmult(t(Xbar),sqrt(W))
        Yptilde<-mmult(t(Y),sqrt(W))
        yb=Y-Xbar%*%solve(tcrossprod(Xptilde),Xptilde%*%t(Yptilde))
        cv2[j2]=sum((yb%*%J_f2)*yb)
      }
      
      W<-exp(-ud/which.min(cv2))
      
    }
    Wbf=W*Wr
    
    Xptilde<-mmult(t(Xbar),sqrt(Wbf))
    Yptilde<-mmult(t(Y),sqrt(Wbf))
    yb=Y-Xbar%*%solve(tcrossprod(Xptilde),Xptilde%*%t(Yptilde))
    cv[j]<-sum((yb%*%J_f2)*yb)
    if(cv[j]<mi) {
      Wb=Wbf
      mi=cv[j]
    }
  }
  return(Wb)  
}


cvW_new4<-function(Ds, #distanze geografiche
                   Xbar,
                   Y, #Ys
                   Wr,
                   J_f2=J_f1){
  m=dim(Ds)[1] #righe della matrice di distanze geo
  cv=1:m
  mi=Inf 
  ##############-------
  for(j in 1:m)
  {
    ud=Ds[j,]
    
    {
      
      h2=1:20
      cv2=h2
      
      for(j2 in h2)
      {
        
        W=exp(-ud/h2[j2])*Wr
        
        yb=c_f3(Xbar,Y,W)
        #if(sum(abs(yb-yb2))>1e-11) browser()
        cv2[j2]=sum((yb%*%J_f2)*yb)
      }
      
      W<-exp(-ud/which.min(cv2))
      
    }
    Wbf=W*Wr
    
    yb=c_f3(Xbar,Y,Wbf)
    
    cv[j]<-sum((yb%*%J_f2)*yb)
    if(cv[j]<mi) {
      Wb=Wbf
      mi=cv[j]
    }
  }
  return(Wb)  
}



