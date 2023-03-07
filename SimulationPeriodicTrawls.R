#This file contains the code to create the simulated
#sample paths from a periodic trawl and trawls with
#additive or multiplicative seasonality as presented
#in the article.

library(ambit)
##Comparing periodic trawl, seasonal+trawl, seasonal*trawl
###########################################################
##Gaussian trawl with exponential trawl and sine periodic function
n<-500

Delta<-0.1


#kernel function
p <- function(x){
  tau <-3
  return(sin(2*pi*x/tau))
}

#identity function
one <- function(x){
  return(1)
}

#Exponential trawl function
trawlfct<-"Exp"
trawlfct_par <-0.5

#Gaussian trawl
distr<-"Gauss"
distr_par<-c(0,1)

#simulate 2*n observations to account for burn-in period of length n
set.seed(1)
periodic_trawl <-sim_weighted_trawl((2*n), Delta, trawlfct, trawlfct_par, distr, distr_par, p)$path[(n+1):(2*n)]
set.seed(1)
trawl <-sim_weighted_trawl((2*n), Delta, trawlfct, trawlfct_par, distr, distr_par, one)$path[(n+1):(2*n)]


a_trawl <-numeric(n)
m_trawl <-numeric(n)
seasonalweights <-numeric(n)

for(i in 1:n){
  seasonalweights[i]<-p((i-1)*Delta)
  a_trawl[i] <- seasonalweights[i]+trawl[i]
  m_trawl[i] <- seasonalweights[i]*trawl[i]
}

par(mfrow=c(4,2))
plot(trawl, type ="l", main="Gaussian")
acf(trawl)

plot(a_trawl, type ="l", main="Gaussian")
acf(a_trawl)

plot(m_trawl, type ="l", main="Gaussian")
acf(m_trawl)


plot(periodic_trawl, type ="l", main="Gaussian")
acf(periodic_trawl)

plot(seasonalweights, type ="l", main="Gaussian")
acf(seasonalweights)

par(mfrow=c(1,1))
####Creating the ggplots
library(lubridate)
library(ggplot2)
library(dplyr)
library(latex2exp)

theme_update(text = element_text(size=30)) #Increase font size in labels
#Plot as time series with correct times
time_seq<-seq(1,500,1)


##Trawl
trawl_df <- data.frame(time=time_seq,
                       value=trawl)
p <- ggplot(trawl_df, aes(x=time,y=value))+
  geom_line()+
  xlab(TeX("$t$"))+
  ylab(TeX("$X_t$"))
p
ggsave("X.eps", width = 20, height = 20, units = "cm")
#Plot the acf
my_acf <- acf(trawl, lag=30, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab(TeX("ACF of $X_t$"))
q
ggsave("ACF_X.eps", width = 20, height = 20, units = "cm")

##Trawl + seasonality
trawl_df <- data.frame(time=time_seq,
                       value=a_trawl)
p <- ggplot(trawl_df, aes(x=time,y=value))+
  geom_line()+
  xlab(TeX("$t$"))+
  ylab(TeX("$X_t^a$"))
p
ggsave("Xa.eps", width = 20, height = 20, units = "cm")
#Plot the acf
my_acf <- acf(a_trawl, lag=30, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab(TeX("ACF of $X^a_t$"))
q
ggsave("ACF_Xa.eps", width = 20, height = 20, units = "cm")


####
##Trawl * seasonality
trawl_df <- data.frame(time=time_seq,
                       value=m_trawl)
p <- ggplot(trawl_df, aes(x=time,y=value))+
  geom_line()+
  xlab(TeX("$t$"))+
  ylab(TeX("$X_t^m$"))
p
ggsave("Xm.eps", width = 20, height = 20, units = "cm")
#Plot the acf
my_acf <- acf(m_trawl, lag=30, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab(TeX("ACF of $X^m_t$"))
q
ggsave("ACF_Xm.eps", width = 20, height = 20, units = "cm")


##periodic Trawl
trawl_df <- data.frame(time=time_seq,
                       value=periodic_trawl)
p <- ggplot(trawl_df, aes(x=time,y=value))+
  geom_line()+
  xlab(TeX("$t$"))+
  ylab(TeX("$Y_t$"))
p
ggsave("Y.eps", width = 20, height = 20, units = "cm")
#Plot the acf
my_acf <- acf(periodic_trawl, lag=30, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab(TeX("ACF of $Y_t$"))
q
ggsave("ACF_Y.eps", width = 20, height = 20, units = "cm")


##periodic fct
trawl_df <- data.frame(time=time_seq,
                       value=seasonalweights)
p <- ggplot(trawl_df, aes(x=time,y=value))+
  geom_line()+
  xlab(TeX("$t$"))+
  ylab(TeX("$q(t)$"))
p
ggsave("q.eps", width = 20, height = 20, units = "cm")
#Plot the acf
my_acf <- acf(seasonalweights, lag=30, plot = FALSE)
my_acfdf <- with(my_acf, data.frame(lag, acf))

#Confidence limits
alpha <- 0.95
conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(n)

q <- ggplot(data = my_acfdf, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(yintercept=conf.lims, lty=2, col='blue') +
  xlab("Lag")+
  ylab(TeX("ACF of $q(t)$"))
q
ggsave("ACF_q.eps", width = 20, height = 20, units = "cm")





