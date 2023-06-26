##########----- approssimazione funzionle --------####################
## dati := a matrix pXn, weher p are the time and n are the sites
## type_of_basis := "fourier" := type of basis to use
## parameters_basis := list(k=9) := a named list of parameters for the basis
## plots := T if you want to see the plots
### return
## Xc := the matrix of coefficent, nXk 
## J := the matrix of inner product, kXk
fv=function(dati, 
            type_of_basis="fourier", 
            parameters_basis= list(k=NULL), 
            plots=F)
{
  x=seq(0,1, length.out = length(dati[,1]))
  if(plots) matplot(x, dati, type="l")
  if (type_of_basis=="fourier"){
    k=parameters_basis$k
    b = create.fourier.basis(range(x), nbasis=k)
  }
  if (type_of_basis=="bspline"){
    k=parameters_basis$k
    b = create.bspline.basis(range(x), nbasis=k)
  }
  df1=Data2fd( dati, argvals=x, basisobj=b)
  if(plots) plot(df1)
  X1=df1$coefs
  W <- eval.penalty(b, rng=range(x),Lfdobj=0)
  return(list(Xc=t(X1), J=W))
}



##########----- approssimazione funzionle variabili --------####################
## lp := a list of dataframes where the first is the dependent functional variable
## type_of_basis := "fourier" := type of basis to use
## parameters_basis := list(k1=9, k2=list()) := a named list of parameters for the basis
## plots := T if you want to see the plots
### return
## Y := matrix of coefficent of response variable
## X := list of matrix of coefficent of covariate
## Jx := list of matrix of inner product of basis function of covariate
## Jy := matrix of inner product of basis function of response
step1<-function(lp,
            type_of_basis="fourier",
            parameters_basis=list(k1=9, k2=10),
            plots=T #if you want to see the plots
){
  #LOAD AND FUNCTIONALIZE input data
  ### Response variable -----  
  buff<-lp[[1]]
  
  dati<-as.matrix(buff)
  k=parameters_basis$k1
  buff=fv(dati, "bspline", parameters_basis=list(k=k), plots)
  Y=buff$Xc
  Jy=buff$J
  #### HO CREATO La Y ---- 
  
  ###### FUNCTIONAL COVARIATES  -------
  ### 
  if( length(parameters_basis$k2)==1) 
    k=parameters_basis$k2 else
      k=parameters_basis$k2
  X<-list()
  Jx=list()
  #browser()
  for(covariates in 2:length(lp)){
    buff=lp[[covariates]] 
    dati<-as.matrix(buff)
    buff<-fv(dati, "bspline", parameters_basis=list(k=k[[covariates-1]]), plots)
    X[[covariates-1]]=buff$Xc
    Jx[[covariates-1]]=buff$J
  }
    #----------
  return(list(Y=Y,X=X, Jx=Jx, Jy=Jy))
}













































# library(vows)
# ##########----- approssimazione funzionle --------####################
# ## dat := a matrix pXn, weher p are the time and n are the sites
# ## type_of_basis := "fourier" := type of basis to use
# ## parameters_basis := list(k=9, types=1, t) := a named list of parameters for the basis
# ## plots := T if you want to see the plots
# fv2=function(dati, 
#             type_of_basis="fourier", 
#             parameters_basis= list(k=NULL, t=NULL), 
#             plots){
#   x=seq(0,1, length.out = length(dati[,1]))
#   if(plots) matplot(x, dati, type="l")
#   if (type_of_basis=="fourier"){
#     k=parameters_basis$k
#     b = create.fourier.basis(range(x), nbasis=k)
#   }
#   if (type_of_basis=="bspline"){
#     k=parameters_basis$k
#     b = create.bspline.basis(range(x), nbasis=k)
#   }
#   n<-ncol(dati)
#   B = eval.basis(x, b)
#   W <- eval.penalty(b, rng=range(x),Lfdobj=0)
#   arrayVp=array(NA, c(k, k, n))
#   coef=matrix(NA, k, n)
#   for (i in 1:n){
#     fit0=gam(dati[,i]~B-1)
#     arrayVp[,,i]=fit0$Vp
#     coef[,i]=fit0$coefficients
#   }
#   B.grid=eval.basis(t., bsb)
#   Y=B.grid%*%t(coef)
#   if(plots) matplot(Y, type = "l")
#   
#   return(list(Xc=coef, J=W, arrayVp=arrayVp))
# }
