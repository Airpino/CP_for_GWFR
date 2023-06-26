##########----- funzioni per il gwr: la zeta --------####################
## X := list of covariate, the generic element is a matrix nXki where n is the number
##      of sites and ki is the number of coefficent of covariate
## J := list Jphi, the generic element is the matrix of inner prodoct of functional basis
##      of dimension kiXki where ki is the number of coefficent of covariate
## fl := a flag, 1 to include b0 and 0 to not include
### return
## Z := return the matrix Z
zeta <- function(X, J, fl){
  n <- nrow(X[[1]])
  l <- length(X)
  if(fl){
    Z <- matrix(1, nrow = n, 1)
    for (i in 1:l) {
      Z <- cbind(Z, X[[i]]%*%J[[i]])
    }
  } else{
    Z <- X[[1]]%*%J[[1]]
    for (i in 2:l) {
      Z <- cbind(Z, X[[i]]%*%J[[i]])
    }
  }
  return(Z)
}

##########----- funzioni per il gwr: i beta --------####################
## Z := the Z matrix of covariate 
## Y := the matrix of coefficente 
## w := the vector of weight
### return 
## beta :=  la matrice dei beta
## hat_riga := il vettore riga della matrice hat
Bs <- function(Z, Y, w){
  A<-tryCatch(
    {
      solve(t(Z)%*%(w*Z))
    },
    error=function(cond) {
      ginv(t(Z)%*%(w*Z))
    })
  C <- (A%*%t(Z)*w)
  B <- C%*%Y
  return(list(beta=B, hat_riga=C))
}
##########----- funzione gwr --------####################
## X := list of covariate, the generic element is a matrix nXki where n is the number
##      of sites and ki is the number of coefficent of covariate
## Y := the matrix of coefficente 
## J := list Jphi, the generic element is the matrix of inner prodoct of functional basis
##      of dimension kiXki where ki is the number of coefficent of covariate
## Jy := inner prodoct of functional basis of covariete 
## DMAT := the matrix of distance
## kernel := type of kernel
## h := the list of paramiter of type of kernel 
## fl := a flag, 1 to include b0 and 0 to not include
### return 
## beta :=  the list of la matrix of beta, each element for each site
## hat_Matrix := the hat matrix
## residuo := il residuo del modello
## Yh := the Y estimeted
wr <- function(X, Y, J, Jy, DMAT, kernel="Gaus", h=list(h=1), fl){
  n <- nrow(Y)
  Z <- zeta(X, J, fl)
  beta <- list()
  Hat_Matrix <- matrix(0, n, n)
  for (i in 1:n){
    wi <- peso(DMAT, i, kernel, h)
    buff <- Bs(Z, Y, wi)
    beta[[i]] <- buff$beta
    Hat_Matrix[i,] <- Z[i,]%*%buff$hat_riga
  }
  Yh <- Hat_Matrix%*%Y
  yb <- Y-Yh 
  r <- sqrt(sum((yb%*%Jy)*yb))
  return(list(beta=beta, Hat_Matrix=Hat_Matrix, residuo=r, Yh=Yh))
}

wr1 <- function(i, Zs, Ys, Ds, kernel = "Gaus", h=list(h=1)){
  wi <- peso(Ds, i, kernel, list(h=a$the_best_h))
  buff <- Bs(Zs, Ys, wi)
  y <- Zs[i,]%*%buff$beta
  return(list(yihat=y, betai=buff$beta))
}