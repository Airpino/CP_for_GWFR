library(MASS)
library(fda)

source("script_funzioni_funzionale.R")
source("script_funzioni_gw.R")
source("script_funzioni_pesi.R")

n=200
X=list()
k=list()
####---- dati generati molro random e a caso ----####
for (i in 1:4) {
  set.seed(i)
  k[[i]]=sample(9:15, 1)
  X[[i]]=mvrnorm(175, rep(0, n), diag(n)  )
}

k=list(k1=k[[1]], k2=list(k[[2]], k[[3]], k[[4]]))
coord=mvrnorm(n, rep(0, 2), diag(2)  )
####----- come si usa step1-----#######

buff=step1(X,
           type_of_basis="bspline",
           parameters_basis=k,
           plots=T)

Y=buff$Y
X=buff$X
J=buff$Jx
Jy=buff$Jy
DMAT=as.matrix(dist(coord))
fl=1
#####----- cross validation h -----#####
a=cvh(X, Y, J, Jy, DMAT, kernel="Gaus", fl)

#####----- per avere solo una stima -----#####
Z=zeta(X, J, fl)
buff=wr1(2, Z, Y, DMAT, kernel = "Gaus", h=list(h=a$the_best_h))
y2=buff$yihat
betai=buff$betai

#####----- per avere tutto il modello -----#####
buff=wr(X, Y, J, Jy, DMAT, kernel = "Gaus", h=list(h=a$the_best_h), fl)
Yhat=buff$Yh
beta=buff$beta

