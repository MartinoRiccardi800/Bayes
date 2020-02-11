rm(list=ls())
library(rstan)
library(mice)
library(readr)
library(coda)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)

setwd("C:/Users/jacop/Documents/Universita/QUINTO ANNO/PRIMO_SEMESTRE/BAYESIAN_STATISTICS/Project/Martino/csv_files/csv_files")


# Datasets import
coo=read.csv("coordinates.csv", header = TRUE, sep = ",", dec = ".")
cp=read.csv("cumProd.csv", header = TRUE, sep = ",", dec = ".")
rp=read.csv("rateProd.csv", header = TRUE, sep = ",", dec = ".")

# New dataframe with ordered oil wells
cpp=data.frame(cp[,1],cp[,12],cp[,23],cp[,34],cp[,45],cp[,47:50],cp[,2:11],cp[,13:22],cp[,24:33],cp[,35:44],cp$PRODU50,cp$Time)
labels(cpp)
# Divide the oil wells according to the time they have been created (3 different temporal instants)
fase=rep(1,1,50)
p=dim(cp)[2]
for(i in 1:p-1){
  if (length(which(cpp[,i]==0))>37) {
    if (length(which(cpp[,i]==0))>74){
      fase[i]=3
    }
    else {
      fase[i]=2
    }
  }
}
n1=length(which(fase==1))
n2=length(which(fase==2))
n3=length(which(fase==3))
n1+n2+n3

# Colors for the plot
colore=rep(1,1,50)
colore[which(fase==1)]=rainbow(3)[1]  #red
colore[which(fase==2)]=rainbow(3)[2]  #green
colore[which(fase==3)]=rainbow(3)[3]  #blue

# Plot 1
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
legend(6000, 2000, legend="Phase 1", col="red", pch=17, cex=0.8)
# Plot 2
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
points(coo$x[which(fase==2)],coo$y[which(fase==2)],col=colore[which(fase==2)],pch=17)
legend(6000, 2000, legend=c("Phase 1", "Phase 2"), col=c("red", "green"), pch=17, cex=0.8)
# Plot 3
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
points(coo$x[which(fase==2)],coo$y[which(fase==2)],col=colore[which(fase==2)],pch=17)
points(coo$x[which(fase==3)],coo$y[which(fase==3)],col=colore[which(fase==3)],pch=17)
legend(6000, 2000, legend=c("Phase 1", "Phase 2", "Phase 3"), col=c("red", "green", "blue"), pch=17, cex=0.8)


# Cumulative oil extraction over time
x11()
matplot(cpp, type = 'l', lty = 1, xlab = "Time", ylab = "Oil quantity",col = colore)

#MIA PARTE
c=seq(2,26,by=1) # vettore che contiene i numeri delle colonne di cpp...
# ...con i pozzi nel primo slot

r=seq(1,38,by=1) # istanti di tempo del primo slot 

size_cpp=dim(cpp)
pozzi=coo[c,1] # pozzi del primo slot
#View(pozzi) # vedi pozzi

m=length(pozzi) # numero pozzi attivi nel primo slot
tempo=cpp[r,51] # vettore con tutte le date fino al 38° istante

##t contiene tutti i numeri di giorni cumulati fino al 38° istante
t=0;
for (i in 1:length(tempo)-1){
  k=difftime(strptime(tempo[i+1], format="%Y-%m-%d"), strptime(tempo[i], format="%Y-%m-%d"), units="days")
  t=rbind(t,t[i]+k)
}
######
### n contiene tutti 25 volte il numero 38 (ogni pozzo ha 38 istanti di tempo)
n=c(rep(length(t), length(pozzi)))
###
cp=cpp[r,c] #prendo da cpp solo le colonne e righe che ci servono 

##Metto tutte le funzioni cumulate in un unico vettore (nell'ordine dei pozzi)
cumul=NULL
for (i in 1:m){
  cumul=c(cumul,cp[,i])
}
####
###Stessa cosa per gli istanti di tempo
tt=c(rep(t,length(pozzi)));

poid=NULL
for(j in 1:m) 
{ 
  po<-c(rep(j,length(t)))
  poid=c(poid,po)
}
m=as.integer(m)
poid=as.integer(poid)

X <- tt
X=as.matrix(X)
k=1
G=matrix(rep(0,950*25),nrow = 950)
for(i in 1:25){
  G[k:(k+37),i]=rep(1,38)
  k=k+38
}
colnames(G) <- c("2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")
Y <- as.vector(cumul)

N <- length(X)
p_fix <- 2
ngr <- 25

##-------------------------------------
## MODELLO 4 Hierarchical Regression con DIPENDENZA SPAZIALE


# Matrice delle distanze tra i pozzi (Distanza Euclidea)
# Costruisco la matrice D(50,50)
nrow=seq(1,50,1)
ncol=nrow
D=matrix(rep(0,2500),nrow = 50)
for(i in ncol){
  for(j in ncol){
    D[i,j]=sqrt((coo$x[i]-coo$x[j])^2+(coo$y[i]-coo$y[j])^2)
  }
  D[i,i]=0
}
# Tolgo il primo pozzo (INJE1)
Dmod=D[-1,-1]
Dmod
# Costruisco la matrice "Near" dei neighboors per ciascun pozzo: ogni riga i rappresenta 
# le distanze tra il pozzo i e tutti gli altri (compreso se stesso, infatti la prima colonna è zero)
Near=matrix(rep(0,625),nrow = 25)
for(i in 1:25){
  Near[i,]=sort(D[i,1:25])
}
# Estraggo da Near la sottomatrice Near10 che rappresenta i 10 nearest neighbor per ogni pozzo
Near10=Near[,c(2:11)]
dim(Near10)
nx1=dim(Near10)[1]
nx2=dim(Near10)[2]
# Calcolo la distanza media dei 10 nearest neighboor per ogni pozzo
avdist=rowMeans(Near10)
avdist


###################################### Priors

# tau2--->cov(rate production)
RP=cbind(rp$PRODU2,rp$PRODU3,rp$PRODU4,rp$PRODU5,rp$PRODU6,rp$PRODU7,rp$PRODU8,rp$PRODU9,rp[,2:11],rp[,13:19])
varRP=rep(0,25)
for(i in 1:25){
  varRP[i]=var(RP[1:38,i])
}
sigma2_betaPR=mean(varRP)
sigma2_betaPR

# theta
meanRP=colMeans(RP[1:38,])
fit=lm(meanRP ~ -1+avdist )
BETA.LS=fit$coef
BETA.LS
S2.LS=summary(fit)$sigma^2
S2.LS
# scelgo i parametri della inv_gamma s.t. la media sia S2.LS (la media di una invgamma(a,b) è mu=b/(a-1),
# quindi scelgo per esempio a=2, b=25000


## MODELLO 4 Hierarchical Regression CON dipendenza spaziale

data_POZZI4 <-list(N = N, 
                   ngr = 25,  
                   Y = Y, 
                   X = as.vector(tt),
                   D = as.vector(avdist),
                   G = as.matrix(G))

inits4 <- function() 
{
  list(beta = rep(1400, ngr),
       sigma2_beta = rep(20000, ngr),
       tau2= 2,
       theta=0.8,
       sigma2_theta=20000)
}

# FIT THE MODEL

Normal_REGR_Spatial <- stan(file = "Normal_REGR_Spatial_ST.stan", 
                            data = data_POZZI4,
                            chains = 2, 
                            iter = 6000, 
                            warmup = 1000, 
                            thin= 10, 
                            seed = 42, 
                            init = inits4,
                            algorithm = 'NUTS')

save(fitNormal_REGR_Spatial, file="Normal_REGR_Spatial")
#load("Normal_REGR_est.dat")

x11()
rstan::traceplot(Normal_REGR_Spatial, pars = "sigma2_beta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial, pars = "beta", inc_warmup = FALSE)
x11()
rstan::traceplot(Normal_REGR_Spatial, pars = "sigma2_theta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial, pars = "theta", inc_warmup = FALSE)

x11()
plot_post4 <- Normal_REGR_Spatial %>% 
  rstan::extract("beta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post4 %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

x11()
plot_post4 <- Normal_REGR_Spatial %>% 
  rstan::extract("theta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post4 %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# Diagnostic ------------------------------------------------------------

coda_chain4 <- As.mcmc.list(Normal_REGR_Spatial, pars = c("beta", "theta"))
summary(coda_chain4)

# Gelman and rubin's convergence diagnostic
gelman.diag(coda_chain4, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain4, frac1=0.1, frac2=0.5)

# autocorrelation and plot
x11()
acfplot(coda_chain4, lag.max = 30)

my_WAIC <- function(fit, param){
  llik   <- rstan::extract(fit, param)[[1]]
  p_WAIC <- sum(apply(llik, 2, var))
  lppd   <- sum(apply(llik, 2, function(x) log(mean(exp(x)))))
  WAIC   <- - 2 * lppd + 2 * p_WAIC
  return(WAIC)
}

m4=my_WAIC(Normal_REGR_Spatial, "log_lik")


# PLOT (PROVA)

vett=rainbow(25)
colvet=NULL
for(k in 1:25){
  colvet=c(colvet,rep(vett[k],38))
}
x11()
plot(tt,Y,pch=20, col=colvet)
for(k in 1:25){
abline(0,mean(jpar4[,k]),col=vett[k])
}



# Provo togliendo sigma2_theta


data_POZZI4bis <-list(N = N, 
                   ngr = 25,  
                   Y = Y, 
                   X = as.vector(tt),
                   D = as.vector(avdist),
                   G = as.matrix(G))

inits4bis <- function() 
{
  list(beta = rep(1400, ngr),
       sigma2_beta = rep(20000, ngr),
       tau2= 20000,
       theta=0.72)
}

# FIT THE MODEL

Normal_REGR_Spatialbis <- stan(file = "Normal_REGR_Spatial_ST2.stan", 
                            data = data_POZZI4bis,
                            chains = 2, 
                            iter = 30000, 
                            warmup = 3000, 
                            thin= 10, 
                            seed = 42, 
                            init = inits4bis,
                            algorithm = 'NUTS')


save(Normal_REGR_Spatialbis, file="fitNormal_spatbia.dat")
load("fitNormal_spatbia.dat")

x11()
rstan::traceplot(Normal_REGR_Spatialbis, pars = "sigma2_beta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatialbis, pars = "beta", inc_warmup = FALSE)
x11()
rstan::traceplot(Normal_REGR_Spatialbis, pars = "theta", inc_warmup = FALSE)


x11()
plot_post4bis <- Normal_REGR_Spatialbis %>% 
  rstan::extract("beta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post4bis %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

x11()
plot_post4bis <- Normal_REGR_Spatialbis %>% 
  rstan::extract("theta") %>% 
  as.data.frame() %>% 
  map_df(as_data_frame, .id = 'param')

plot_post4bis %>% 
  ggplot(aes(value, fill = param)) + 
  geom_density() + 
  facet_wrap(~param, scales = 'free') + 
  scale_fill_locuszoom() + 
  theme_minimal() +
  theme(legend.position="none")

# Diagnostic ------------------------------------------------------------

coda_chain4bis <- As.mcmc.list(Normal_REGR_Spatialbis, pars = c("beta", "theta"))
summary(coda_chain4bis)

# Gelman and rubin's convergence diagnostic
gelman.diag(coda_chain4bis, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain4bis, frac1=0.1, frac2=0.5)

# autocorrelation and plot
x11()
acfplot(coda_chain4bis, lag.max = 30)

my_WAIC(Normal_REGR_Spatialbis, "log_lik")
