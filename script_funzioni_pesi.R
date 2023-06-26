##########----- funzione vettore peso --------####################
## X := list of coefficient matrices of the covariates, 
##       the generic element is a matrix nXki where n is the number
##       of sites and ki is the number of coefficients of covariate
## DMAT := the matrix of geographic distances
## si := site 
## kernel := type of kernel
## par_kernel := list of parameters for the kernel
### return 
## w := the vector of weigths
peso <- function(DMAT, si, kernel="Gaus", par_kernel=list(h=NULL))
{
  u <- DMAT[si,]
  if(kernel=="Gaus")
  {
    h <- par_kernel$h
    w <- exp(-0.5*u^2/h)
  }
  return(w)
}

##########----- funzione cross validation --------####################
## X := list of covariate, the generic element is a matrix nXki where n is the number
##      of sites and ki is the number of coefficent of covariate
## Y := the matrix of coefficente 
## J := list Jphi, the generic element is the matrix of inner prodoct of functional basis
##      of dimension kiXki where ki is the number of coefficent of covariate
## Jy := inner prodoct of functional basis of covariete 
## DMAT := the matrix of distance
## kernel := type of kernel
## fl := a flag, 1 to include b0 and 0 to not include
## h := the vector of h to be checked
### return 
## the_best_h := the best h 
## cvh := the vextor of cv for each h (in {1,2,....,19,20})
cvh <- function(X, Y, J, Jy, 
                DMAT, 
                kernel="Gaus", 
                fl, 
                h=c(1:20)){
  n <- nrow(Y)
  Z <- zeta(X, J, fl)
  min_crit <- Inf
  cv=1:length(h)
  
  for (j in 1:length(h)) {
    Yh <- matrix(0, n, ncol = ncol(Y))
    for (i in 1:n){
      wi <- peso(DMAT, i, kernel, list(h=h[j]))
      buff <- Bs(Z[-i,], Y[-i,], wi[-i])
      Yh[i,] <- Z[i,]%*%buff$beta
    }
    yb <- Y-Yh 
    r <- (sum((yb%*%Jy)*yb))
    cv[j] <- r
    if(min_crit >r)
      opt_h <- h[j]
  }
  return(list(the_best_h=opt_h, cvh=cv))
}