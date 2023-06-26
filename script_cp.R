rm(list=ls(all=TRUE))
library(fda)

source("script_funzioni_funzionale.R")
source("script_funzioni_gw.R")
source("script_funzioni_pesi.R")
source("script_funzioni_plot.R")


n=200
X=list()
k=list()
####---- dati generati molro random e a caso ----####
for (i in 1:4) {
  set.seed(i)
  k[[i]]=sample(9:15, 1)
  X[[i]]=mvrnorm(175, mvrnorm(n, 0, 1e+6), 1e+10*diag(n)  )
}

k=list(k1=k[[1]], k2=list(k[[2]], k[[3]], k[[4]]))
coord=mvrnorm(n, rep(0, 2), diag(2)  )

####----- come si usa step1-----#######
buff=step1(X,
           type_of_basis="bspline",
           parameters_basis=k,
           plots=F)

Y=buff$Y
X=buff$X
J=buff$Jx
Jy=buff$Jy
fl=1
DMAT=as.matrix(dist(coord))

bs=create.bspline.basis(c(0,1), k$k1)
bse=eval.basis(seq(0,1, length.out=75) , bs)

#####----- cross validation h -----#####
a=cvh(X, Y, J, Jy, DMAT, kernel="Gaus", fl)

#####----- conformal prediction con una banda -----#####
Z=zeta(X, J, fl)
set.seed(0)
train <- sample(nrow(coord), floor(0.5*nrow(coord)), replace = FALSE)
test <- (1:nrow(coord))[-train]
l=length(test)
R=matrix(0, 1, l-1)

start<-Sys.time()

index_X=c(test[1], train)
index_Y=index_X
Zs=Z[index_X,]
Ys=Y[index_Y,]
Ds=DMAT[index_X,index_X]
y0=wr1(1, Zs, Ys, Ds, kernel = "Gaus", h=list(h=a$the_best_h))$yihat
#da definire
S=sqrt(colMeans(((Ys - matrix(y0, nrow = l+1, ncol = k$k1, byrow = T))%*%t(bse))^2))
for (j in (2:l)) {
  index_Y=c( test[j], train)
  Ys=Y[index_Y,]
  y=wr1(1, Zs, Ys, Ds, kernel = "Gaus", h=list(h=a$the_best_h))$yihat
  R[j-1]=max((abs(y0- y)%*%t(bse))/S, na.rm =T)
  cat(" - ",j)
}

end<-Sys.time()
cat("\n")
print(end-start)
if(1){tempo_trasc=(end-start)}else{tempo_trasc=tempo_trasc+(end-start)}
print(paste0("estimated end at ", Sys.time()+tempo_trasc))

Yf=bse%*%t(Y)
p1=CURVE_PLO(Yf,TIME=seq(0,1,length.ou=75),xlab = TeX("$\\tau$"),ylab=TeX("$Y(\\tau)$"),
               showTIT = T)
show(p1)


alpha=0.05
arg=c("t", "Prediction band", paste(expression(D), "- SUP"))
yf=y0%*%t(bse)
r=sqrt(quantile(R, 1-alpha))
I=seq(0,1,length.ou=75)
p1<-BAND_PLO(Yf[, ],TIME=I,I, yf, r,S,xlab=arg[1],ylab=arg[2],TITLE_lab=arg[3],showTIT=T)
show(p1)
