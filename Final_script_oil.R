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
library(car)
library(BayesFactor)
library(gplots)

setwd("C:/Users/jacop/Documents/Universita/QUINTO ANNO/PRIMO_SEMESTRE/BAYESIAN_STATISTICS/Project/Martino/csv_files/csv_files")


# Datasets import
coo=read.csv("coordinates.csv", header = TRUE, sep = ",", dec = ".")
cp=read.csv("cumProd.csv", header = TRUE, sep = ",", dec = ".")
rp=read.csv("rateProd.csv", header = TRUE, sep = ",", dec = ".")

## Exploratory Analysis ------------------

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

# Plot 1 - Location of the wells opened firstly
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
legend(6000, 2000, legend="Phase 1", col="red", pch=17, cex=0.8)
# Plot 2 - - Location of the wells opened secondly
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
points(coo$x[which(fase==2)],coo$y[which(fase==2)],col=colore[which(fase==2)],pch=17)
legend(6000, 2000, legend=c("Phase 1", "Phase 2"), col=c("red", "green"), pch=17, cex=0.8)
# Plot 3 - Location of the wells opened thirdly
x11()
plot(coo$x[which(fase==1)],coo$y[which(fase==1)],col=colore[which(fase==1)],pch=17,xlim = c(0,7500),ylim = c(0,7500), xlab = "X coordinate", ylab="Y coordinate", title("Oil well coordinates"))
points(coo$x[which(fase==2)],coo$y[which(fase==2)],col=colore[which(fase==2)],pch=17)
points(coo$x[which(fase==3)],coo$y[which(fase==3)],col=colore[which(fase==3)],pch=17)
legend(6000, 2000, legend=c("Phase 1", "Phase 2", "Phase 3"), col=c("red", "green", "blue"), pch=17, cex=0.8)


# Plot of the Cumulative oil extraction over time
x11()
matplot(cpp, type = 'l', lty = 1, xlab = "Time", ylab = "Oil quantity",col = colore)


p_fix = 2
ngr = 2

## Preliminary quantities for the STAN model---------------
# t (time)
# y (response: extrated oil)


# Converting the time variable from "date" to "day"
tempotot=cpp[,dim(cpp)[2]]
# At time t1=tempotot[39] and t2=tempotot[76] new wells are opened
t1=39
t2=76
ttot=length(tempotot) # Total number of temporal instants

t=0;
for (i in 1:length(tempotot)-1){
  k=difftime(strptime(tempotot[i+1], format="%Y-%m-%d"), strptime(tempotot[i], format="%Y-%m-%d"), units="days")
  t=rbind(t,t[i]+k)
}


# Extract y for all initial n1 wells (the well number 1 in the dataset is not considered because never working)
y=NULL
for (i in 2:(n1+1)){
  y=c(y,cpp[,i])
}

# Try to fit a frequentist linear regression to get the linear behaviour or not
reg1=t
reg1=rep(reg1,n1)

reg2=rep(0,length(t))
for (i in t1:length(t)){
  reg2[i]=t[i]-t[t1-1]
}

reg2=rep(reg2,n1)
reg3=rep(0,length(t))
for (i in t2:length(t)){
  reg3[i]=t[i]-t[t2-1]
}

reg3=rep(reg3,n1)
mat=data.frame(cbind(y,reg1,reg2,reg3))
attach(mat)
fit1=lm(y~-1+reg1+reg2+reg3,data=mat)
summary(fit1)
fit2=lm(y~-1+reg1+I(reg2^2)+I(reg3^2),data=mat)
summary(fit2)
detach(mat)

x11()
plot(fit2)

# Computing the dummy matrix (useful for the STAN regression file)
# We need 3 parameters for each well (3 because there are three instants of time when new wells are opened)

N <- 117*25  # Number of observations of the response variable
ngr <- 75 # Number of parameters

XX <- t
XX=as.matrix(XX)
k=1
GG=matrix(rep(0,N*ngr),nrow = N)
for(i in 0:(n1-1)){
  GG[k:(k+(ttot-1)),3*i+1]=rep(1,ttot)
  k=k+ttot
}
k=1
for(i in 0:(n1-1)){
  GG[(k+(t1-1)):(k+(ttot-1)),3*i+2]=rep(1,79)
  k=k+ttot
}
k=1
for(i in 0:(n1-1)){
  GG[(k+ngr):(k+(ttot-1)),3*i+3]=rep(1,42)
  k=k+ttot
}
full_dummy2=GG
full_dummy=GG[1:ttot,1:3]


# Computing the average distance matrix: it is filled with the Euclidean distances between all the wells

nrow=seq(1,(dim(cpp)[2])-1,1)
ncol=nrow
D=matrix(rep(0,2500),nrow =50)
for(i in ncol){
  for(j in ncol){
    D[i,j]=sqrt((coo$x[i]-coo$x[j])^2+(coo$y[i]-coo$y[j])^2)
  }
  D[i,i]=0
}
# Delete well 1 (INJE1), the well not working
Dmod=D[-1,-1]
Dmod

# Building the matrix "Near" at each step when new wells are opened
# NB: the distance between a well and itself is later removed

Near=matrix(rep(0,n1*n1),nrow = n1)
for(i in 1:n1){
  Near[i,]=sort(Dmod[i,1:n1])
}
# Extracting the 10 nearest neighbors for each well (for the first 25 wells)
Near10=Near[,c(2:11)]
dim(Near10)
nx1=dim(Near10)[1]
nx2=dim(Near10)[2]
# Computing the average distange (for the first 25 wells)
avdist=rowMeans(Near10)
avdist

Near1=matrix(rep(0,t1*t1),nrow = t1)
for(i in 1:t1){
  Near1[i,]=sort(Dmod[i,1:t1])
}
# Extracting the 10 nearest neighbors for each well for the first 25 wells after opening new 14 wells
Near15=Near1[,c(2:11)]
dim(Near15)
nx1=dim(Near15)[1]
nx2=dim(Near15)[2]
# Computing the average distange (for the first 25 wells) after opening new 14 wells
avdist1=rowMeans(Near15)
avdist1

Near2=matrix(rep(0,49*49),nrow = 49)
for(i in 1:49){
  Near2[i,]=sort(Dmod[i,1:49])
}
# Extracting the 10 nearest neighbors for each well for the first 25 wells after opening the last 10 wells
Near20=Near2[,c(2:11)]
dim(Near20)
nx1=dim(Near20)[1]
nx2=dim(Near20)[2]
# Computing the average distange (for the first 25 wells) after openingthe last 10 wells
avdist2=rowMeans(Near20)
avdist2
avgdist=rep(0,n1*3)
for ( i in (0:(n1-1)))
  avgdist[3*i+1]=avdist[i+1]
for ( i in (0:(n1-1)))
  avgdist[3*i+2]=avdist1[i+1]
for ( i in (0:(n1-1)))
  avgdist[3*i+3]=avdist2[i+1]

# Now we are interested in the difference in the average distance
for ( i in (0:(n1-1)))
  avgdist[3*i+2]=avgdist[3*i+2]-avgdist[3*i+1]
for ( i in (0:(n1-1)))
  avgdist[3*i+3]=avgdist[3*i+3]-(avgdist[3*i+1]+avgdist[3*i+2])


# New Y vector in order to have a less computationally heavy STAN file
yy=y

for(z in 1:n1){
    p=1+(z-1)*ttot
    i=t1+(z-1)*ttot
    j=(t2-1)+(z-1)*ttot
    h=ttot+(z-1)*ttot
    yy[c(i:j)]=y[c(i:j)]-y[i-1]
    yy[c((j+1):h)]=y[c((j+1):h)]-y[j]
    yy[c(p:(i-1))]=y[c(p:(i-1))]
}

tempo=reg1
for(z in 1:n1){
  p=1+(z-1)*ttot
  i=t1+(z-1)*ttot
  j=(t2-1)+(z-1)*ttot
  h=ttot+(z-1)*ttot
  tempo[c(i:j)]=reg1[c(i:j)]-reg1[i-1]
  tempo[c(j+1:h)]=reg1[c(j+1:h)]-reg1[j]
  #tempo[c(p:(i-1))]=reg1[c(p:(i-1))]
}
tempo=rep(tempo[1:ttot],n1)

x11()
plot(tempo[1:ttot],yy[1:ttot])
points(tempo[1:ttot],yy[1:ttot])


## Priors selection --------------

# The idea is to take the last opened wells and use them to extract a prior for theta

# tau2--->cov(rate production)
RP3=cbind(rp[,c(36:44,46)])
varRP3=rep(0,10)
for(i in 1:10){
  varRP3[i]=var(RP3[76:117,i])
}
sigma2_betaPR3=mean(varRP3)
sigma2_betaPR3
meanRP3=colMeans(RP3[76:117,])
fit3=lm(meanRP3 ~ +avdist2[41:50])
summary(fit3)
fit3=lm(meanRP3 ~ -1+avdist2[41:50])
summary(fit3)
BETA.LS3=fit3$coef
BETA.LS3
S2.LS3=summary(fit3)$sigma^2
S2.LS3
# In the model choose the parameters theta and sigma2_beta acordingly to the frequentist linear regression
# For theta choose ~ N(1,1)
# For sigma_beta ~ InvGamma(a,b) s.t the mean mu=b/(a-1) is around s2.LS3


## Fit of the STAN model for the first 25 wells, the one starting to extract oil at time t=0 ---------

data_POZZI_end <-list(N = N, 
                      ngr = n1*3,  
                      Y = yy, 
                      X = as.vector(tempo),
                      D = as.vector(avgdist),
                      G = as.matrix(full_dummy))

# Initial conditions
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
                               chains = 2, 
                               iter = 20000, 
                               warmup = 2000, 
                               thin= 10, 
                               seed = 42, 
                               init = inits4end,
                               algorithm = 'NUTS')

#save(Normal_REGR_Spatialbis, file="Normal_REGR_Spatialbis_est.dat")

paramsend=rstan::extract(Normal_REGR_Spatial_end, c("beta", "theta"), perm = T)
jparend=cbind(paramsend[[1]], paramsend[[2]])

# Traceplots and density plots
x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "sigma2_beta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end, pars = "beta", inc_warmup = FALSE)
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
plot_post4 <- Normal_REGR_Spatial_end %>% 
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

## Diagnostic for model fitted with the first 25 wells -------

coda_chain4 <- As.mcmc.list(Normal_REGR_Spatial_end, pars = c("beta", "theta"))
summary(coda_chain4)

# Gelman and rubin's convergence diagnostic
gelman.diag(coda_chain4, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain4, frac1=0.1, frac2=0.5)

# autocorrelation and plot
x11()
acfplot(coda_chain4, lag.max = 30)


## Fit of the STAN model for the first successive 14 wells, the one starting to extract oil at time t=t1 ---------

# In order to fit the model also for the second slot of wells I extract the new y (response) and distance vectors
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

# Initial conditions
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
                                chains = 2, 
                                iter = 20000, 
                                warmup = 2000, 
                                thin= 10, 
                                seed = 42, 
                                init = inits4end2,
                                algorithm = 'NUTS')

paramsend2=rstan::extract(Normal_REGR_Spatial_end2, c("beta", "theta"), perm = T)
jparend2=cbind(paramsend2[[1]], paramsend2[[2]])

x11()
rstan::traceplot(Normal_REGR_Spatial_end2, pars = "sigma2_beta", inc_warmup = TRUE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end2, pars = "beta", inc_warmup = FALSE)
x11()
rstan::traceplot(Normal_REGR_Spatial_end2, pars = "theta", inc_warmup = FALSE)

x11()
plot_post4 <- Normal_REGR_Spatial_end2 %>% 
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
plot_post4 <- Normal_REGR_Spatial_end2 %>% 
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

## Diagnostic for model fitted with the second 14 wells -------

coda_chain42 <- As.mcmc.list(Normal_REGR_Spatial_end2, pars = c("beta", "theta"))
summary(coda_chain42)

# Gelman and rubin's convergence diagnostic
gelman.diag(coda_chain42, confidence = 0.95,  autoburnin = TRUE, multivariate=TRUE)

# Geweke's convergence diagnostic
geweke.diag(coda_chain42, frac1=0.1, frac2=0.5)

# autocorrelation and plot
x11()
acfplot(coda_chain42, lag.max = 30)


## Tests --------

#Test for equal thetas
theta1=jparend[,dim(jparend)[2]]
theta2=jparend2[,dim(jparend2)[2]]
ttestBF(x=theta1,y=theta2)

# Checking if the total downgrade in terms of amount of extracted oil by the first 25 wells is absorbed
# by the the extraction of the new 14 wells
betamean=colMeans(jparend)
b2vect=NA
for(i in 0:24){
  b2vect[i+1]=betamean[3*i + 2]
}
betamean2=colMeans(jparend2)
b2vect2=NA
for(i in 0:13){
  b2vect2[i+1]=betamean2[3*i + 2]
}
bayes.t.test(x=b2vect,y=14/25*b2vect2)


# Plotting densities of thetas to discover in a qualitative way the possible differences
d1 = density(theta1)
x11()
plot(d1,xlim=c(0.4,1.2),main = expression(paste("Distributions of ", theta, "'s")))
d2 = density(theta2)
points(d2,type="l")
tr=seq(0,15,by=1)
yt=rep(0.87,16)
points(yt,tr,col="red",type="l")

v1=seq(0,15,by=0.01)
v2=rep(0.87,1501)
points(v2,v1, type = "l", col="red")

betaavg=colMeans(jparend)
betaavg=betaavg[-length(betaavg)]
betaavg2=colMeans(jparend2)
betaavg2=betaavg2[-length(betaavg2)]

s1=rep(c(1,2,3),25)
s2=rep(c(1,2,3),14)
x11()
par(mfrow=c(1,2))
boxplot(betaavg~s1, xlab = expression(paste(beta, " index")), ylab = expression(beta), main = expression(paste(beta, "'s of the first 25 wells" )))
boxplot(betaavg2~s2, xlab = expression(paste(beta, " index")), ylab = expression(beta), main = expression(paste(beta, "'s of the second 14 wells" )))

#Checking if beta_1,j of the second model are zeros
beta12 <- as.matrix(paramsend2$beta)
dim(beta12)
# compute the 95% posterior credible interval for beta12
CI_beta12 = apply(beta12, 2, quantile, c(0.025, 0.975)) 
CI_beta12

idx_cov_BL = NULL
for(l in 1:42){
  if(CI_beta12[1,l]<0 && CI_beta12[2,l]>0)
  {
    cat("*** variable ", labels(beta)[[2]][l], " excluded \n")
  }
  else
  {
    cat("*** variable ", labels(beta)[[2]][l], " included \n")
    idx_cov_BL = c(idx_cov_BL, l)
  }
  
}

mean_beta12_post <- apply(beta12, 2, "mean")
mean_beta12_post
idx=c(1,4,7,10,13,16,19,22,25,28,31,34,37,40)
lidx=length(idx)

x11()
plotCI(x = idx, y = mean_beta12_post[idx], liw = (-CI_beta12[1,idx] + mean_beta12_post[idx]),  
       uiw = (CI_beta12[2,idx]- mean_beta12_post[idx]),
       type = "p", lwd = 1.5, main=expression(paste("Decision intervals for ", beta, "1")), ylab = "", xlab = "")
abline(h = 0, col = "blue")

