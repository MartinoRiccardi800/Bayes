# Bayesian Statistics
Negative Spatial Correlation for extracted oil

## Description
Our project consists in a spatial analysis of a dataset which is geographically referenced and temporally correlated. A bayesian model will be built in order to study the negative correlation between cumulative oil extractions, taking into account the geographical closeness of the wells. Our goal is trying to see if the amount of oil, gathered up in a certain location, depends on the activation of other wells nearby.

## Installation
All the packages used are:

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

## Data and code
The repository contains two datasets (cumulative oil extracted and rates) and a code to be run.

## Model
Y_ij = beta1j\*t_ij + beta2j\*(t_ij-t1) + beta3j\*(t_ij-t2) +  epsilon_ij; 
epsilon_ij~ N(0, tau^2)
## Results
Observing the boxplots of betas, the final results seem to show that  the "extraction coefficient" seems to be affected by opening other wells nearby, causing a decrease of the production
