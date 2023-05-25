rm(list=ls(all=TRUE))
source("fun_used.R", encoding = 'UTF-8')

#####   figure degli inquinanti
load("Dati_applicazione_indquinanti.RData")
#####cose che servono
t=1:365
k=12
bs=create.bspline.basis(range(t), k)
bse=eval.basis(t, bs)
#####PM10
PM10=idt12[[1]][[1]]
X=bse%*%t(PM10)
matplot(X, type = "l", xlab = "Day", ylab = "PM10")
#####SO2
SO2=idt12[[1]][[2]]
X=bse%*%t(SO2)
matplot(X, type = "l", xlab = "Day", ylab = expression(SO[2]) )
######PM2.5
PM2.5=idt12[[2]]
X=bse%*%t(PM2.5)
matplot(X, type = "l", xlab = "Day", ylab = "PM2.5")


X=list(idt12[[1]][[1]], idt12[[1]][[2]])
Y=idt12[[2]]
D=idt12[[3]]
Wq=idt12[[4]]
J_f1=idt12[[5]]
J_f2=J_f1
mn=dim(Y)

DMAT=D
#DMAT<-DMAT/(max(DMAT))*15

H=length(X)
M=ncol(X[[1]])
Xbar=matrix(0, nrow = nrow(Y), ncol =M*H )

for(i in 1:H){  
  Xbar[,(1+(M*(i-1))):(i*M)]=(X[[i]]%*%J_f1)
}
flag=1
fl=flag
############## ------Metto il vettore unitario se flag=1 ----------------
if(flag){
  unos=matrix(1, nrow=nrow(Y), ncol=1)
  
  Xbar_OLD=Xbar
  Xbar=cbind(unos, Xbar)
}

lol=1

n=mn[1]
z=0
np=100
x=seq(0,1, length.out = np)
bs=create.bspline.basis(range(x), k)
bse=eval.basis(x, basisobj=bs) 
set.seed(z)
smpl=sample(1:n, 15)
trial=(1:n)[-smpl]
#browser()
###########################################
l=length(trial)
###########################################
alpha=0.05

Yb1=matrix(0, l, k)
cp1=matrix(0, l, l)
rk1=matrix(0, l, l)
stb1=matrix(0, l, np)

Yb2=matrix(0, l, k)
cp2=matrix(0, l, l)
rk2=matrix(0, l, l)
stb2=matrix(0, l, np)

Yb3=matrix(0, l, k)
cp3=matrix(0, l, l)
rk3=matrix(0, l, l)
stb3=matrix(0, l, np)

Yb4=matrix(0, l, k)
cp4=matrix(0, l, l)
rk4=matrix(0, l, l)
stb4=matrix(0, l, np)

Yb5=matrix(0, l, k)
cp5=matrix(0, l, l)
rk5=matrix(0, l, l)
stb5=matrix(0, l, np)
#DMAT<-as.matrix(dist(coord))
ll=length(smpl)+1
######### MEGA FOR ----------
for(t in 1:l){
  index=c(smpl, trial[t])
  Xs=list()
  for (i in 1:H) {  Xs[[i]]=(X[[i]][index,])}
  Ys=(Y[index,])
  
  start<-Sys.time()
  print(t)
  Wr=rep(1,nrow(Ys))
    W = cvW_new4(DMAT[index,index], #distanze geografiche
                 Xbar[index,],
                 Y[index,], #Ys
                 Wr,
                 J_f2=J_f1)
    
  w=W
  res=lmr_1.1(Ys, Xs, w, #vettore
              J_f1, flag=fl)
  
  tmpM<-res$HatMatrix%*%Ys
  ym=colSums(tmpM)/ll
  st=sqrt(colSums(((tmpM-ym)%*%t(bse))^2)/ll)
  ####################  HOMO misura di non-conf legato ai pesi del modello D -----------------------
  stb1[t,]=st
  
  r2=matrix(0, ll, 1)
  yd=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
  st_tmp=sqrt(colSums((((Ys)-ym)%*%t(bse))^2)/ll)
  
  for(i in 1:ll){yd[i,]=yd[i,]/st_tmp}
  yb=t(Data2fd(x, t(yd), bs)$coef)
  r2=rowSums((yb%*%Wq)*yb)
  
  cp1[t, t]=length(r2[r2<r2[ll]])<=alpha*ll
  rk1[t, t]=r2[ll]
  Yb1[t,]=ym
  
  #################### HOMO misura di non-conf legato ai pesi del modello DH -----------------------
  stb2[t,]=st
  
  r2=rowSums((yb%*%J_f2)*yb)/w
  
  cp2[t, t]=length(r2[r2<r2[ll]])<=alpha*ll
  rk2[t, t]=r2[ll]
  Yb2[t,]=ym
  
  ####################  HOMO misura di non-conf (VANTINI) Dsup -----------------------
  yd=(tmpM-ym)%*%t(bse) #residui sui siti scelti (campione)
  st_tmp=matrix(0, np, 1)
  for(i in 1:np){st_tmp[i]=max(abs((yd)[,i]))}
  fint=function(y){
    B = eval.basis(y, bs)
    coeff=Data2fd(x, st_tmp, bs)$coefs
    return(t(B%*%coeff))
  }
  a=integrate(fint, lower=0, upper=1)[[1]]
  st_tmp=st_tmp/a
  stb3[t,]=st_tmp
  
  yd=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
  for(i in 1:np){st_tmp[i]=max(abs((yd)[,i]))}
  fint=function(y){
    B = eval.basis(y, bs)
    coeff=Data2fd(x, st_tmp, bs)$coefs
    return(t(B%*%coeff))
  }
  a=integrate(fint, lower=0, upper=1)[[1]]
  st_tmp=st_tmp/a
  for(i in 1:ll){r2[i]=max(abs(yd[i,]/st_tmp))}
  
  cp3[t, t]=length(r2[r2<r2[ll]])<=alpha*ll
  rk3[t, t]=r2[ll]
  Yb3[t,]=ym
  
  ####################  ETERO Dh-----------------------
  ym=colSums(w*(res$HatMatrix%*%Ys))/sum(w)
  st_E=sqrt(colSums(w*(((res$HatMatrix%*%Ys)-ym)%*%t(bse))^2)/sum(w))
  stb4[t,]=st_E
  
  r2=matrix(0, ll, 1)
  yd_1=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
  st_tmp=sqrt(colSums(w*((Ys-ym)%*%t(bse))^2)/sum(w))
  
  for(i in 1:ll){yd_1[i,]=yd_1[i,]/st_tmp}
  yb=t(Data2fd(x, t(yd_1), bs)$coef) #coefficienti della trasformazione funz degli scarti dalla media standardizzati
  r2=rowSums((yb%*%J_f2)*yb)/w
  
  cp4[t, t]=length(r2[r2<r2[ll]])<=alpha*ll
  rk4[t, t]=r2[ll]
  Yb4[t,]=ym
  
  ####################  ETERO Dsup-----------------------
  yd=(tmpM-ym)%*%t(bse)
  st_tmp=matrix(0, np, 1)
  for(i in 1:np){st_tmp[i]=max(abs((yd)[,i]))}
  fint=function(y){
    B = eval.basis(y, bs)
    coeff=Data2fd(x, st_tmp, bs)$coefs
    return(t(B%*%coeff))
  }
  a=integrate(fint, lower=0, upper=1)[[1]]
  st_tmp=st_tmp/a
  stb5[t,]=st_tmp
  
  r2=matrix(0, ll, 1)
  yd=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
  for(i in 1:np){st[i]=max(abs((yd)[,i]))}
  fint=function(y){
    B = eval.basis(y, bs)
    coeff=Data2fd(x, st_tmp, bs)$coefs
    return(t(B%*%coeff))
  }
  a=integrate(fint, lower=0, upper=1)[[1]]
  st_tmp=st_tmp/a
  for(i in 1:ll){r2[i]=max(abs(yd[i,]/st_tmp))}
  
  cp5[t, t]=length(r2[r2<r2[ll]])<=alpha*ll
  rk5[t, t]=r2[ll]
  Yb5[t,]=ym
  
  ####################  -----------------------
  for(j in (1:l)[-t]){
    
    cat(" - ",j)
    index2=c(smpl, trial[j])
    Ys=(Y[index2,])
    
    W = cvW_new4(DMAT[index,index], #distanze geografiche
                 Xbar[index,],
                 Y[index2,], #Ys
                 Wr,
                 J_f2=J_f1)
    
    w=W
    res=lmr_1.1(Ys, Xs, w, #vettore
                J_f1, flag=fl)
    ym=colSums((res$HatMatrix%*%Ys))/ll
    ####################  HOMO con D-----------------------
    r2=matrix(0, ll, 1)
    yd_1=tcrossprod((Ys-matrix(ym,nrow(Ys),ncol(Ys),byrow = T)),bse) #residui sui siti scelti (campione)
    yd=yd_1
    st=sqrt(colSums((((Ys)-ym)%*%t(bse))^2)/ll)
    for(i in 1:ll){
      yd[i,]=yd[i,]/st}
    
    yb=t(Data2fd(x, t(yd), bs)$coef)
    r2=rowSums((yb%*%Wq)*yb)  
    cp1[j, t]=length(r2[r2<r2[ll]])<=alpha*ll
    rk1[j, t]=r2[ll]
    
    ####################  HOMO con Dh-----------------------
    r2=rowSums((yb%*%J_f2)*yb)/w
    
    cp2[j, t]=length(r2[r2<r2[ll]])<=alpha*ll
    rk2[j, t]=r2[ll]
    
    ####################  HOMO con Dsup -----------------------
    yd=yd_1
    for(i in 1:np){st[i]=max(abs((yd)[,i]))}#ricalcolo st
    fint=function(y){
      B = eval.basis(y, bs)
      coeff=Data2fd(x, st, bs)$coefs
      return(t(B%*%coeff))
    }
    a=integrate(fint, lower=0, upper=1)[[1]]
    st=st/a
    for(i in 1:ll){r2[i]=max(abs(yd[i,]/st))}
    
    cp3[j, t]=length(r2[r2<r2[ll]])<=alpha*ll
    rk3[j, t]=r2[ll]
    
    ####################  ETERO con Dh-----------------------
    ym=colSums(w*(res$HatMatrix%*%Ys))/sum(w)
    
    r2=matrix(0, ll, 1)
    yd=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
    st=sqrt(colSums(w*((Ys-ym)%*%t(bse))^2)/sum(w))
    
    for(i in 1:ll){
      yd[i,]=yd[i,]/st
    }
    yb=t(Data2fd(x, t(yd), bs)$coef)
    
    r2=rowSums((yb%*%J_f2)*yb)/w    
    cp4[j, t]=length(r2[r2<r2[ll]])<=alpha*ll
    rk4[j, t]=r2[ll]
    
    #################### ETERO con Dsup -----------------------
    r2=matrix(0, ll, 1)
    yd=t(t(Ys)-ym)%*%t(bse) #residui sui siti scelti (campione)
    for(i in 1:np){st[i]=max(abs((yd)[,i]))}
    fint=function(y){
      B = eval.basis(y, bs)
      coeff=Data2fd(x, st, bs)$coefs
      return(t(B%*%coeff))
    }
    a=integrate(fint, lower=0, upper=1)[[1]]
    st=st/a
    for(i in 1:ll){r2[i]=max(abs(yd[i,]/st))}
    
    cp5[j, t]=length(r2[r2<r2[ll]])<=alpha*ll
    rk5[j, t]=r2[ll]
    
    ####################  -----------------------
  }
  
  end<-Sys.time()
  cat("\n")
  print(end-start)
  if(t==1){tempo_trasc=(end-start)}else{tempo_trasc=tempo_trasc+(end-start)}
  print(paste0("estimated end at ", Sys.time()+(l-t)*tempo_trasc/t))
}

covt=matrix(0, 5, 100)
I=matrix(0, 1, 4)
IP=matrix(0, 1, 4)
IV=matrix(0, 1, 4)
NIP=matrix(0, 1, 4)
NIV=matrix(0, 1, 4)
bande=list()

buff=index_band3(stb1, rk1, Yb1, cp1, c("t", "PM2.5", paste("Prediction band - ", expression(D), "- HOMO")))
I[1,2:4]=buff[[1]]
covt[1,]=buff[[2]]
bande[[1]]=list(centro=buff[[3]], raggio=buff[[4]], modulazione=buff[[5]])
I[1,1]=mean(covt[1,])
buff=index_band3(stb2, rk2, Yb2, cp2, c("t", "PM2.5", paste("Prediction band - ", expression(D[h]), "- HOMO")))
IP[1,2:4]=buff[[1]]
covt[2,]=buff[[2]]
IP[1,1]=mean(covt[2,])
bande[[2]]=list(centro=buff[[3]], raggio=buff[[4]], modulazione=buff[[5]])
buff=index_band3(stb3, rk3, Yb3, cp3, c("t", "PM2.5", paste("Prediction band - ", expression(D[supp]), "- HOMO")))
IV[1,2:4]=buff[[1]]
covt[3,]=buff[[2]]
IV[1,1]=mean(covt[3,])
bande[[3]]=list(centro=buff[[3]], raggio=buff[[4]], modulazione=buff[[5]])
buff=index_band3(stb4, rk4, Yb4, cp4, c("t", "PM2.5", paste("Prediction band - ", expression(D[h]), "- HETERO")))
NIP[1,2:4]=buff[[1]]
covt[4,]=buff[[2]]
NIP[1,1]=mean(covt[4,])
bande[[4]]=list(centro=buff[[3]], raggio=buff[[4]], modulazione=buff[[5]])
buff=index_band3(stb5, rk5, Yb5, cp5, c("t", "PM2.5", paste("Prediction band - ", expression(D[supp]), "- HETERO")))
NIV[1,2:4]=buff[[1]]
covt[5,]=buff[[2]]
NIV[1,1]=mean(covt[5,])
bande[[5]]=list(centro=buff[[3]], raggio=buff[[4]], modulazione=buff[[5]])
dfres<-as.data.frame(
  rbind(I[1,], IP[1,], IV[1,]))
rownames(dfres)<-c("D","Dh","Dsup")
names(dfres)<-c("CovAlL %", "CovAlG %","Width","InterScore")
dfres2<-as.data.frame(
  rbind(NIP[1,], NIV[1,])
)
rownames(dfres2)<-c("Dh","Dsup")
names(dfres2)<-c("CovAlL %", "CovAlG %","Width","InterScore")

print(dfres)
print(dfres2)

save.image("app_inq_15.RData")







