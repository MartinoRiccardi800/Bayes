# Negative Spatial Correlation for extracted oil

## Description
Our project consists in a spatial analysis of a dataset which is geographically referenced and temporally correlated. A bayesian model will be built in order to study the negative correlation between cumulative oil extractions, taking into account the geographical closeness of the wells. Our goal is trying to see if the amount of oil, gathered up in a certain location, depends on the activation of other wells nearby.

## Installation
All the packages used are:<br/>

install.packages("rstan")<br/>
install.packages("mice")<br/>
install.packages("readr")<br/>
install.packages("coda")<br/>
install.packages("ggplot2")<br/>
install.packages("tidyr")<br/>
install.packages("dplyr")<br/>
install.packages("purrr")<br/>
install.packages("ggsci")<br/>
install.packages("car")<br/>
install.packages("BayesFactor")<br/>
install.packages("devtools")<br/>
library(devtools)<br/>
install_github("rasmusab/bayesian_first_aid")<br/>

## Data and code
The repository contains 3 datasets (cumulative oil extracted, rates and coordinates of wells) and a code to be run.

## Model
Y_ij = beta1j\*t_ij + beta2j\*(t_ij-t1) + beta3j\*(t_ij-t2) +  epsilon_ij;<br/>
epsilon_ij~ N(0, tau^2)
## Results
Observing the boxplots of betas, the final results seem to show that  the "extraction coefficient" seems to be affected by opening other wells nearby, causing a decrease of the production
