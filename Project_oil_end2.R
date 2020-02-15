rm(list=ls())

install.packages("rstan")
install.packages("mice")
install.packages("readr")
install.packages("coda")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("purrr")
install.packages("ggsci")
install.packages("car")
install.packages("BayesFactor")
install.packages("devtools")
library(devtools)
install_github("rasmusab/bayesian_first_aid")

library(rstan)
library(mice)
library(readr)
library(coda)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggsci)
library(car)
library(BayesFactor)
library(BayesianFirstAid)

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


######

p_fix <- 2
ngr <- 25


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



######################################
tempotot=cpp[,51]

##t contiene tutti i numeri di giorni cumulati fino al 38° istante
t=0;
for (i in 1:length(tempotot)-1){
  k=difftime(strptime(tempotot[i+1], format="%Y-%m-%d"), strptime(tempotot[i], format="%Y-%m-%d"), units="days")
  t=rbind(t,t[i]+k)
}

# Extract y for all initial 25 wells
y=NULL
for (i in 2:26){
  y=c(y,cpp[,i])
    
}

# Try to fit a frequentist linear regression to get the linear behaviour or not
reg1=t
reg1=rep(reg1,25)

reg2=rep(0,length(t))
for (i in 39:length(t)){
  reg2[i]=t[i]-t[38]

  
}
reg2=rep(reg2,25)
reg3=rep(0,length(t))
for (i in 76:length(t)){
  reg3[i]=t[i]-t[75]
  
  
}
reg3=rep(reg3,25)
mat=data.frame(cbind(y,reg1,reg2,reg3))
attach(mat)
fit1=lm(y~-1+reg1+reg2+reg3,data=mat)
summary(fit1)
fit2=lm(y~-1+reg1+I(reg2^2)+I(reg3^2),data=mat)
summary(fit2)
detach(mat)

x11()
plot(fit2)

# Computing the dummy matrix useful for the stan regression file
XX <- t
XX=as.matrix(XX)
k=1
GG=matrix(rep(0,2925*75),nrow = 2925)
for(i in 0:24){
  GG[k:(k+116),3*i+1]=rep(1,117)
  k=k+117
}
k=1
for(i in 0:24){
  GG[(k+38):(k+116),3*i+2]=rep(1,79)
  k=k+117
}
k=1
for(i in 0:24){
  GG[(k+75):(k+116),3*i+3]=rep(1,42)
  k=k+117
}
full_dummy=GG
full_dummy=full_dummy[1:117,1:3]


N <- 117*25
ngr <- 75

# Computing the average distance matrix
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
  Near[i,]=sort(Dmod[i,1:25])
}
# Estraggo da Near la sottomatrice Near10 che rappresenta i 10 nearest neighbor per ogni pozzo
Near10=Near[,c(2:11)]
dim(Near10)
nx1=dim(Near10)[1]
nx2=dim(Near10)[2]
# Calcolo la distanza media dei 10 nearest neighboor per ogni pozzo
avdist=rowMeans(Near10)
avdist

Near1=matrix(rep(0,39*39),nrow = 39)
for(i in 1:39){
  Near1[i,]=sort(Dmod[i,1:39])
}
# Estraggo da Near la sottomatrice Near10 che rappresenta i 10 nearest neighbor per ogni pozzo
Near15=Near1[,c(2:11)]
dim(Near15)
nx1=dim(Near15)[1]
nx2=dim(Near15)[2]
# Calcolo la distanza media dei 10 nearest neighboor per ogni pozzo
avdist1=rowMeans(Near15)
avdist1

Near2=matrix(rep(0,49*49),nrow = 49)
for(i in 1:49){
  Near2[i,]=sort(Dmod[i,1:49])
}
# Estraggo da Near la sottomatrice Near10 che rappresenta i 10 nearest neighbor per ogni pozzo
Near20=Near2[,c(2:11)]
dim(Near20)
nx1=dim(Near20)[1]
nx2=dim(Near20)[2]
# Calcolo la distanza media dei 10 nearest neighboor per ogni pozzo
avdist2=rowMeans(Near20)
avdist2
avgdist=rep(0,75)
for ( i in (0:24))
  avgdist[3*i+1]=avdist[i+1]
for ( i in (0:24))
  avgdist[3*i+2]=avdist1[i+1]
for ( i in (0:24))
  avgdist[3*i+3]=avdist2[i+1]

# Try with differences of avgdist
for ( i in (0:24))
  avgdist[3*i+2]=avgdist[3*i+2]-avgdist[3*i+1]
for ( i in (0:24))
  avgdist[3*i+3]=avgdist[3*i+3]-(avgdist[3*i+1]+avgdist[3*i+2])


# New Y vector in order not to change the stan file
yy=y

for(z in 1:25){
    p=1+(z-1)*117
    i=39+(z-1)*117
    j=75+(z-1)*117
    h=117+(z-1)*117
    yy[c(i:j)]=y[c(i:j)]-y[i-1]
    yy[c((j+1):h)]=y[c((j+1):h)]-y[j]
    yy[c(p:(i-1))]=y[c(p:(i-1))]
}

tempo=reg1
for(z in 1:25){
  p=1+(z-1)*117
  i=39+(z-1)*117
  j=75+(z-1)*117
  h=117+(z-1)*117
  tempo[c(i:j)]=reg1[c(i:j)]-reg1[i-1]
  tempo[c(j+1:h)]=reg1[c(j+1:h)]-reg1[j]
  #tempo[c(p:(i-1))]=reg1[c(p:(i-1))]
}
tempo=rep(tempo[1:117],25)

x11()
plot(tempo[1:117],yy[1:117])
points(tempo[1:117],yy[1:117])

data_POZZI_end <-list(N = N, 
                      ngr = 75,  
                      Y = yy, 
                      X = as.vector(tempo),
                      D = as.vector(avgdist),
                      G = as.matrix(full_dummy))

inits4end <- function() 
{
  list(beta = rep(1400, ngr),
       sigma2_beta = rep(20000, ngr),
       tau2= 20000,
       theta=0.72)
}

# FIT THE MODEL

Normal_REGR_Spatial_end <- stan(file = "Normal_REGR_Spatial_ST2_final7.stan", 
                               data = data_POZZI_end,
                               chains = 1, 
                               iter = 3000, 
                               warmup = 1000, 
                               thin= 10, 
                               seed = 42, 
                               init = inits4end,
                               algorithm = 'NUTS')

#save(Normal_REGR_Spatialbis, file="Normal_REGR_Spatialbis_est.dat")

paramsend=rstan::extract(Normal_REGR_Spatial_end, c("beta", "theta"), perm = T)
jparend=cbind(paramsend[[1]], paramsend[[2]])

x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "sigma2_beta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "beta", inc_warmup = FALSE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "sigma2_theta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "theta", inc_warmup = FALSE)

x11()
plot_post4 <- Normal_REGR_Spatial_end %>% 
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

# In order to fit the model also for the second slot of wells I extract the new y and ditance vectors
y2=NULL
for (i in 27:40){
  y2=c(y2,cpp[,i])
  
}

avgdist2=rep(0,14*3)
for ( i in (0:13))
  avgdist2[3*i+1]=0
for ( i in (0:13))
  avgdist2[3*i+2]=avdist1[i+25+1]
for ( i in (0:13))
  avgdist2[3*i+3]=avdist2[i+25+1]

# Try with differences of avgdist
for ( i in (0:13))
  avgdist2[3*i+2]=avgdist2[3*i+2]-avgdist2[3*i+1]
for ( i in (0:13))
  avgdist2[3*i+3]=avgdist2[3*i+3]-(avgdist2[3*i+1]+avgdist2[3*i+2])


yy2=y2

for(z in 1:14){
  p=1+(z-1)*117
  i=39+(z-1)*117
  j=75+(z-1)*117
  h=117+(z-1)*117
  yy2[c(i:j)]=y2[c(i:j)]-y2[i-1]
  yy2[c((j+1):h)]=y2[c((j+1):h)]-y2[j]
  yy2[c(p:(i-1))]=y2[c(p:(i-1))]
}

x11()
plot(tempo[1:117],yy2[1:117])

ngr2=14*3
N2=length(yy2)
tempo2=tempo[1:N2]

data_POZZI_end2 <-list(N = N2, 
                      ngr = ngr2,  
                      Y = yy2, 
                      X = as.vector(tempo2),
                      D = as.vector(avgdist2),
                      G = as.matrix(full_dummy))

inits4end2 <- function() 
{
  list(beta = rep(600, ngr2),
       sigma2_beta = rep(20000, ngr2),
       tau2= 20000,
       theta=0.72)
}

# FIT THE MODEL

Normal_REGR_Spatial_end2 <- stan(file = "Normal_REGR_Spatial_ST2_final7.stan", 
                                data = data_POZZI_end2,
                                chains = 1, 
                                iter = 3000, 
                                warmup = 1000, 
                                thin= 10, 
                                seed = 42, 
                                init = inits4end2,
                                algorithm = 'NUTS')

paramsend2=rstan::extract(Normal_REGR_Spatial_end2, c("beta", "theta"), perm = T)
jparend2=cbind(paramsend2[[1]], paramsend2[[2]])

#Test for equal thetas
theta1=jparend[,dim(jparend)[2]]
theta2=jparend2[,dim(jparend2)[2]]
ttestBF(x=theta1,y=theta2)

betamean=colMeans(jparend)
b2sum=0
for(i in 0:24){
  b2sum=b2sum+betamean[3*i + 2]
}
betamean2=colMeans(jparend2)
b21sum=0
for(i in 0:13){
  b21sum=b21sum+betamean2[3*i + 2]
}




b2vect=NA
for(i in 0:24){
  b2vect[i+1]=betaavg[3*i+2]
}
b2vect2=NA
for(i in 0:13){
  b2vect2[i+1]=betaavg2[3*i+2]
}
bayes.t.test(x=-b2vect, y=14/25*b2vect2, alternative = "two.sided")
b2intercept=NA
for(i in 0:13){
  b2intercept[i+1]=betaavg2[3*i+1]
}
bayes.t.test(x=b2intercept, mu=0.5)


#test the coefficient of increase
bayes.t.test(x=-b2vect, y=0.85*14/25*b2vect2)
bayes.t.test(x=-b2vect, y=1.02*14/25*b2vect2)
bayes.t.test(x=-b2vect, y=1.07*14/25*b2vect2)
bayes.t.test(x=-b2vect, y=1.1*14/25*b2vect2)

