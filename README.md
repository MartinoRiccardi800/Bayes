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
Y_{i,j} = beta_{1,j}*t_{i,j} + beta_{2,j}*(t_{i,j}-t_{1}) + beta_{3,j}*(t_{i,j}-t_{2}) +  epsilon_{i,j} 
epsilon_{i,j}~ N(0, tau^{2})
## Results
Thr final results seem to show that (thanks to boxplot) the "extraction coefficient" seems to be affected by opening other wells nearby, causing a decrease of the production
