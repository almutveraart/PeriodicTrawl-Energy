#In this supplementary file we explore a GMM estimation 
#of the trawl function parameters.

#Load all relevant packages
library(readxl)
library(ambit)

source("PeriodicFunctions.R")

#Read in the data
Day_ahead_prices_DELU <- read_excel("Day-ahead_prices_DELU.xlsx")

#Remove the header
dailyprices <-Day_ahead_prices_DELU[10:1563,]
my_prices <-as.numeric(unlist(dailyprices[,2]))


#########"Calm" period 01.10.2018-31.12.2020: 823 observations
series1 <- my_prices[1:823]

#########"Volatile" period 01.01.2021-01.01.2023: 731 observations
series2 <- my_prices[824:1554]


#Estimating the periodic trawl kernel parameters
#Calm regime
fit_exp_periodic_trawl(series1, m=1, tau=7, GMMlag=50, Delta=1, plotacf=TRUE,lags=100)
fit_LM_periodic_trawl(series1, m=1, tau=7, GMMlag=50, Delta=1, plotacf=TRUE,lags=100)
#Volatile regime
fit_exp_periodic_trawl(series2, m=1, tau=7, GMMlag=50, Delta=1, plotacf=TRUE,lags=100)
fit_LM_periodic_trawl(series2, m=1, tau=7, GMMlag=50, Delta=1, plotacf=TRUE,lags=100)

