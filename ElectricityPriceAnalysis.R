#Analysis of Electricity Day-ahead prices	
#Region: DE-LU	Germany/Luxembourg
#Time period: Oct 1, 2018, 12:00 AM - Jan 1, 2023, 11:59 PM	
#Last update to data: Feb 21, 2023, 11:45 AM	
#Data source: (c) Bundesnetzagentur | SMARD.de	
#Resolution: Daily base load prices
#https://www.smard.de/en

#Load all relevant packages
library(readxl)
library(vctrs)
library(ggplot2)
library(ambit)
theme_update(text = element_text(size=30)) #Increase font size in labels



#Read in the data
Day_ahead_prices_DELU <- read_excel("Day-ahead_prices_DELU.xlsx")
#View(Day_ahead_prices_DELU)

#Remove the header and plot the time series
dailyprices <-Day_ahead_prices_DELU[10:1563,]
my_prices <-as.numeric(unlist(dailyprices[,2]))

#create dataframe
df <- data.frame(date = as.Date("2018-10-01") + 0:1553,
                 price = my_prices)

#create time series plot
p <- ggplot(df, aes(x=date, y=price)) +
  geom_line()+xlab("Time")+ylab("Price")
#display time series plot
p

postscript("data.eps")
p
dev.off()


#########"Calm" period 01.10.2018-31.12.2020: 823 observations
series1 <- my_prices[1:823]
df <- data.frame(date = as.Date("2018-10-01") + 0:822,
                 price = series1)
#create time series plot
p <- ggplot(df, aes(x=date, y=price)) +
  geom_line()+xlab("Time")+ylab("Price")
#display time series plot
p

postscript("data1.eps")
p
dev.off()

#Histogram
ggplot(df, aes(x=price)) + 
  geom_histogram(color="black", fill="darkgray")+
  xlab("Price")+ylab("Counts")


ggplot(df, aes(x=price)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray12")+
  geom_density(alpha=.6, size=1, fill="gray65")+
  xlab("Price")+ylab("Density")


p<-ggplot(df, aes(x=price)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray65")+
  geom_density(alpha=.6, size=1)+
  xlab("Price")+ylab("Density")
p


postscript("hist1.eps")
p
dev.off()

#Acf plot
my_acf <- acf(series1, lag=100, plot = FALSE)
n <- length(series1)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

p <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab("Autocorrelation")
p
postscript("acf1.eps")
p
dev.off()



#########"Volatile" period 01.01.2021-01.01.2023: 731 observations
series2 <- my_prices[824:1554]
df <- data.frame(date = as.Date("2021-01-01") + 0:730,
                 price = series2)

#create time series plot
p <- ggplot(df, aes(x=date, y=price)) +
  geom_line()+xlab("Time")+ylab("Price")
#display time series plot
p

postscript("data2.eps")
p
dev.off()

#Histogram
ggplot(df, aes(x=price)) + 
  geom_histogram(color="black", fill="darkgray")+
  xlab("Price")+ylab("Counts")


#Scaled histogram/density
ggplot(df, aes(x=price)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray12")+
  geom_density(alpha=.6, size=1, fill="gray65")+
  xlab("Price")+ylab("Density")


p<-ggplot(df, aes(x=price)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray65")+
  geom_density(alpha=.6, size=1)+
  xlab("Price")+ylab("Density")
p

postscript("hist2.eps")
p
dev.off()

#Acf plot
my_acf <- acf(series2, lag=100, plot = FALSE)
n<- length(series2)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

p <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab("Autocorrelation")
p
postscript("acf2.eps")
p
dev.off()

#################
#Estimating the periodic trawl process
#Calm regime

########
#Exponential case
tau <- 7
T <- tau+1
rho1 <- acf(series1)$acf[(1+1)]
rhoT <- acf(series1)$acf[(1+T)]
lambda_hat <- 1/(1-8)*log(rhoT/rho1)

#Report the estimated lambda parameter
round(lambda_hat, digits=3)

lags <- 100
tt <- (0:lags)

acfs1 <- acf(series1, lag=100)
graphics::lines(tt, acf_Exp(tt, lambda_hat), lty =1,col=2, lwd=2)

#acfs1[1] #Gives the first lag of the acf
#acfs1$acf[1] #Gives the 0 lag of the acf
#Estimating the function c for the first tau lags
c<-numeric(tau)
#c[1] <-1 #0-lag in acf
for(i in 0:(tau-1)){
  c[i+1]<-acfs1$acf[(1+i)]*exp(lambda_hat*i)
}

plot(c, type="l")
#Report the estimated c(l*Delta) parameters
round(c, digits=3)

#Produce the vector of the 7-periodic c function
c_long <- vec_rep(c,15)[1:101]
plot(c_long, type="l")

#Compute the estimated acf:
rho_hat <- numeric(101)
for(i in 0:100){
  rho_hat[i+1]<-c_long[i+1]*exp(-lambda_hat*i)
}

plot(rho_hat, type="l")
acf(series1, lag=100)
graphics::lines((0:100), rho_hat, lty =1,col=2, lwd=2)

#Plotting with ggplot
#Confidence limits
n <- length(series1)
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

my_acf <- acf(series1, lag=100, plot = FALSE)
n<- length(series1)
my_acfdf <- with(my_acf, data.frame(lag, acf,rho_hat))

p<-ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf))+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue')+
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_line(aes(x = lag, y = rho_hat),colour="red", size=1) +
  xlab("Lag")+
  ylab("Autocorrelation")
p

postscript("acf-exp-1.eps")
p
dev.off()


#Long memory
alpha <- 1.12 #From GMM estimation
H_hat <- 1+ log(rhoT/rho1)/log((alpha+1)/(alpha+T))

round(H_hat, digits=3)
acf(series1, lag=100)
graphics::lines(tt, acf_LM(tt, alpha, H_hat), lty =1,col=2, lwd=2)

#Estimating the function c for the first tau lags
c<-numeric(tau)
#c[1] <-1 #0-lag in acf
for(i in 0:(tau-1)){
  c[i+1]<-acfs1$acf[(1+i)]*(1+i/alpha)^(H_hat-1)
}
plot(c, type="l")

round(c, digits=3)

c_long <- vec_rep(c,15)[1:101]
plot(c_long, type="l")

#Estimated acf:
rho_hat <- numeric(101)
for(i in 0:100){
  rho_hat[i+1]<-c_long[i+1]*(1+i/alpha)^(1-H_hat)
}


plot(rho_hat, type="l")
acf(series1, lag=100)
graphics::lines(tt, rho_hat, lty =1,col=2, lwd=2)


#Plotting with ggplot
#Confidence limits
n <- length(series1)
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

my_acf <- acf(series1, lag=100, plot = FALSE)
n<- length(series1)
my_acfdf <- with(my_acf, data.frame(lag, acf,rho_hat))

p<-ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf))+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue')+
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_line(aes(x = lag, y = rho_hat),colour="red", size=1) +
  xlab("Lag")+
  ylab("Autocorrelation")
p

postscript("acf-lm-1.eps")
p
dev.off()

########
#Estimating the periodic trawl process
#Volatile regime

########
#Exponential case
tau <- 7
T <- tau+1
rho1 <- acf(series2)$acf[(1+1)]
rhoT <- acf(series2)$acf[(1+T)]
lambda_hat <- 1/(1-8)*log(rhoT/rho1)

#Report the estimated lambda parameter
round(lambda_hat, digits=3)

lags <- 100
tt <- (0:lags)

acfs2 <- acf(series2, lag=100)
graphics::lines(tt, acf_Exp(tt, lambda_hat), lty =1,col=2, lwd=2)

#acfs2[1] #Gives the first lag of the acf
#acfs2$acf[1] #Gives the 0 lag of the acf
#Estimating the function c for the first tau lags
c<-numeric(tau)
#c[1] <-1 #0-lag in acf
for(i in 0:(tau-1)){
  c[i+1]<-acfs2$acf[(1+i)]*exp(lambda_hat*i)
}

plot(c, type="l")
#Report the estimated c(l*Delta) parameters
round(c, digits=3)

#Produce the vector of the 7-periodic c function
c_long <- vec_rep(c,15)[1:101]
plot(c_long, type="l")

#Compute the estimated acf:
rho_hat <- numeric(101)
for(i in 0:100){
  rho_hat[i+1]<-c_long[i+1]*exp(-lambda_hat*i)
}


plot(rho_hat, type="l")
acf(series2, lag=100)
graphics::lines(tt, rho_hat, lty =1,col=2, lwd=2)

#Plotting with ggplot
#Confidence limits
n <- length(series2)
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

my_acf <- acf(series2, lag=100, plot = FALSE)
n<- length(series2)
my_acfdf <- with(my_acf, data.frame(lag, acf,rho_hat))

p<-ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf))+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue')+
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_line(aes(x = lag, y = rho_hat),colour="red", size=1) +
  xlab("Lag")+
  ylab("Autocorrelation")
p

postscript("acf-exp-2.eps")
p
dev.off()


#Long memory

alpha <- 5.297
H_hat <- 1+ log(rhoT/rho1)/log((alpha+1)/(alpha+T))

round(H_hat, digits=3)
acf(series2, lag=100)
graphics::lines(tt, acf_LM(tt, alpha, H_hat), lty =1,col=2, lwd=2)

#Estimating the function c for the first tau lags
c<-numeric(tau)
#c[1] <-1 #0-lag in acf
for(i in 0:(tau-1)){
  c[i+1]<-acfs2$acf[(1+i)]*(1+i/alpha)^(H_hat-1)
}
plot(c, type="l")
round(c, digits=3)

c_long <- vec_rep(c,15)[1:101]
plot(c_long, type="l")

#Estimated acf:
rho_hat <- numeric(101)
for(i in 0:100){
  rho_hat[i+1]<-c_long[i+1]*(1+i/alpha)^(1-H_hat)
}

plot(rho_hat, type="l")
acf(series2, lag=100)
graphics::lines(tt, rho_hat, lty =1,col=2, lwd=2)


#Plotting with ggplot
#Confidence limits
n <- length(series2)
a <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + a)/2)/sqrt(n)

my_acf <- acf(series2, lag=100, plot = FALSE)
n<- length(series2)
my_acfdf <- with(my_acf, data.frame(lag, acf,rho_hat))

p<-ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf))+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(yintercept=conf.lims, lty=2, col='blue')+
  geom_segment(mapping = aes(xend = lag, yend = 0), size=0.8)+
  geom_line(aes(x = lag, y = rho_hat),colour="red", size=1) +
  xlab("Lag")+
  ylab("Autocorrelation")
p

postscript("acf-lm-2.eps")
p
dev.off()

