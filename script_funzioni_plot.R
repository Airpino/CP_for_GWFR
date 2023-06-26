library(tidyverse)
library(latex2exp)
# plf2=function(x, Y, c, Ylim)
# {
#   n=dim(Y)[1]
#   if(is.null(n))
#     n=1
#   xylim=list(xlim=c(min(x), max(x)), ylim=Ylim)
#   plot(NULL, NULL, type = "l", xlim = xylim$xlim, ylim = xylim$ylim, xlab=c[1], ylab=c[2])
#   for(i in 1:n){lines(x, Y[i,])}
#   title(c[3])
# }
# 
# index_band6=function(alpha, ST, rk, YM, Yb, c, x=NULL, grafico=T, indicatori=F)
# {
#   RK=sqrt(quantile(rk,1-alpha))
#   if(is.null(x)){x=1:length(ST)}
#   if(grafico){
#     v=c(min((YM-RK*ST), Yb), max((YM+RK*ST), Yb))
#     plf2(x, Yb, c, v)
#     lines(x, YM, col="blue", lwd = 3)
#     lines(x, (YM+RK*ST), col="red", lwd = 3)  
#     lines(x, (YM-RK*ST), col="red", lwd = 3)  
#   }
#   return(c(RK, ST))
# } 
# 

CURVE_PLO<-function(X,
                    TIME=c(1:365),
                    xlab="DAY",
                    ylab="PM2.5",
                    TITLE_lab="TITLE",showTIT=T,alpha=0.2){
  
  X_df<-as_tibble(X) %>% mutate(Day=TIME) %>% dplyr::select(Day,everything())
  X_df <- X_df %>% pivot_longer(-Day,names_to = "Curve",values_to = "PM2.5")
  p<-ggplot(X_df) + geom_line(aes(x=Day,y=PM2.5,group=as.factor(Curve)),alpha=alpha)+
    theme_bw()+
    xlab(xlab)+ylab(ylab)
  
  if(showTIT){  
    p<-p+ggtitle(TITLE_lab)+
      theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}

BAND_PLO<-function(X,
                   TIME=c(1:365),
                   x,m,r,s,
                   xlab="DAY",
                   ylab="PM2.5",
                   TITLE_lab="TITLE",showTIT=T){
  X_df<-as_tibble(X) %>% mutate(Day=TIME) %>% dplyr::select(Day,everything())
  band2<-data.frame(x=x,
                    ym=as.vector(m)-unname(r)*s,
                    yM=as.vector(m)+unname(r)*s,
                    m=as.vector(m))
  X_df <- X_df %>% pivot_longer(-Day,names_to = "Curve",values_to = "PM2.5")
  p<-ggplot(X_df) + geom_line(aes(x=Day,y=PM2.5,group=as.factor(Curve)),alpha=0.4)+
    theme_bw()+
    geom_ribbon(data=band2, aes(x=x,ymin=ym,ymax=yM),alpha=0.1,color="red",
                linewidth=1.2)+
    geom_line(data=band2, aes(x=x,y=m),color="blue",linewidth=1.2,linetype=1)+
    xlab(xlab)+ylab(ylab)
  
  if(showTIT){  
    p<-p+ggtitle(TITLE_lab)+
      theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}