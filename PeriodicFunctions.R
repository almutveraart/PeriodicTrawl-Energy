##Compute the acf and drop lag-0
my_tsa_acf <-function(x, lag.max=NULL, plot=FALSE){
  n <- length(x)-1
  z <- stats::acf(x,type="correlation", lag.max, plot)
  return(z[1:n])
}


#'@title fourier_sin_cos
#'@details Computes the Fourier approximation of a function using
#' the sin-cos representation for given coefficients
#'  a_0,...a_n, b_1, ...b_n
#'@param x argument for which the acf will be computed (this can be a vector)
#'@param a vector of coefficients a_0,...a_n
#'@param b vector of coefficients b_1,...b_n
#'@param tau period of function
#'@return Fourier approximation or vector of Fourier approximations evaluated in x
#'@export
fourier_sin_cos<-function(x, a, b, tau){

  l <- length(x)
  fct_vec <- numeric(l)

  n<-base::length(b)
  if((n+1)!=base::length(a)){print("a and b do not have the correct length")}
  #a contains a_0, a_1,..a_{n-1}
  #b contains b_1, ...b_{n-1}

  for(j in 1:l){
    tmp <- a[1] #a_0
    for(k in 1:n){
      tmp <- a[k+1]*base::cos(2*pi*k*x[j]/tau)+b[k]*base::sin(2*pi*k*x[j]/tau)
    }
    fct_vec[j]<-tmp
  }

  return(fct_vec)

}

#'@title kernel_trawl_acf_exp_periodic
#'@param x argument for which the acf will be computed (this can be a vector)
#'@param lambda parameter in the trawl function
#'@param tau period in the sine function
#'@param a vector of Fourier coefficients
#'@param b vector of Fourier coefficients
#'@details We compute the acf of a weighted trawl process with exponential trawl
#'function and periodic function
#'@return acf or vector of acfs evaluated in x
#'@export


kernel_trawl_acf_exp_periodic<-function(x, lambda, a, b, tau){
  l <- length(x)
  acf_vec <- numeric(l)
  for(i in 1:l){
    acf_vec[i]<- acf_Exp(x[i], lambda)*fourier_sin_cos(x[i], a, b, tau)
  }
  return(acf_vec)
}



#'@title kernel_trawl_acf_LM_periodic
#'@param x argument for which the acf will be computed (this can be a vector)
#'@param alpha parameter in the trawl function
#'@param H parameter in the trawl function
#'@param tau period in the sine function
#'@param a vector of Fourier coefficients
#'@param b vector of Fourier coefficients
#'@details We compute the acf of a weighted trawl process with supGamma trawl
#'function and periodic function
#'@return acf or vector of acfs evaluated in x
#'@export


kernel_trawl_acf_LM_periodic<-function(x, alpha, H, a, b, tau){
  l <- length(x)
  acf_vec <- numeric(l)
  for(i in 1:l){
    acf_vec[i]<- acf_LM(x[i], alpha, H)*fourier_sin_cos(x[i], a, b, tau)
  }
  return(acf_vec)
}

#'Fits an exponential periodic trawl function to equidistant time series data
#'@name fit_exp_periodic_trawl
#'@param x vector of equidistant time series data
#'@param m the order of the Fourier approximation of the periodic function
#'@param tau the period, if not provided it will be estimated
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param GMMlag lag length used in the GMM estimation, the default is 50
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return lambda: the parameter \eqn{\lambda} in the exponential trawl,
#'tau: the period \eqn{\tau},
#'a: the vector of the Fourier coefficients,
#'b: the vector of the Fourier coefficients.
#'@details The trawl function is parametrised by the parameter \eqn{\lambda > 0}
#'  as follows: \deqn{g(x) = e^{\lambda x},  \mbox{ for }  x \le 0.}
#'  The parameter \eqn{\tau} is estimated using the smoothed periodogram,
#'  the  parameter \eqn{\lambda} is estimated by GMM.
#'@export
fit_exp_periodic_trawl <- function(x, m=1, tau, GMMlag=50, Delta=1, plotacf=FALSE,lags=100){

  if(missing(tau)){
    #Estimating tau using the smoothed periodogram
    p_smooth <- LSTS::smooth.periodogram(x, plot=TRUE)
    p_smooth_max_loc <- nnet::which.is.max(p_smooth$smooth.periodogram)
    tau <-2*pi/p_smooth$lambda[p_smooth_max_loc]*Delta
  }

  my_acf <- my_tsa_acf(x, lag.max=GMMlag)#TSA::acf(x,lag.max=GMMlag, plot=F)

  fit_exp_periodic_trawl_foroptim <- function(y){

    lambda <- y[1]


    a<-y[2:(2+m)]
    b<-y[(2+m+1):(2+2*m)]

    lag <- GMMlag
    lss <- 0


    for(i in 1:lag)
    {
      Cor <- kernel_trawl_acf_exp_periodic(Delta*i, lambda, a, b, tau)
      lss <- lss+(my_acf$acf[i] - Cor)^2
    }
    lss
  }


  o <- DEoptim::DEoptim(fit_exp_periodic_trawl_foroptim,c(0,(-100+numeric(2*m+1))),c(100,(100+numeric(2*m+1))),control=DEoptim::DEoptim.control(itermax = 1000, trace = FALSE))

  lambda <- o$optim$bestmem[1]
  a<-o$optim$bestmem[2:(2+m)]
  b<-o$optim$bestmem[(2+m+1):(2+2*m)]


  if(plotacf){
    tt <- (1:lags)
    #TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    my_tsa_acf(x, lag.max=lags, plot=TRUE)
    graphics::lines(tt, kernel_trawl_acf_exp_periodic(tt*Delta, lambda, a, b, tau), lty =1,col=2, lwd=2)
  }

  return(list("lambda"=lambda, "tau"=tau, "avector"=a, "bvector=", b))

}



#'Fits a supGamma periodic trawl function to equidistant time series data
#'@name fit_LM_periodic_trawl
#'@param x vector of equidistant time series data
#'@param m the order of the Fourier approximation of the periodic function
#'@param tau the period, if not provided it will be estimated
#'@param Delta interval length of the time grid used in the time series, the
#'  default is 1
#'@param GMMlag lag length used in the GMM estimation, the default is 50
#'@param plotacf binary variable specifying whether or not the empirical and
#'  fitted autocorrelation function should be plotted
#'@param lags number of lags to be used in the plot of the autocorrelation
#'  function
#'@return alpha: the parameter \eqn{\alpha} in the supGamma trawl,
#'H: the parameter \eqn{H} in the supGamma trawl,
#'tau: the period \eqn{\tau}, a: the vector of the Fourier coefficients,
#'b: the vector of the Fourier coefficients.
#'@details The trawl function is parametrised by the parameter \eqn{\lambda > 0}
#'  as follows: \deqn{g(x) = (1+x/\alpha)^H,  \mbox{ for }  x \le 0.}
#'  The parameter \eqn{\tau} is estimated using the smoothed periodogram,
#'  the other parameters are estimated by GMM.
#'@export
fit_LM_periodic_trawl <- function(x, m=1, tau, GMMlag=50, Delta=1, plotacf=FALSE,lags=100){

  if(missing(tau)){
    #Estimating tau using the smoothed periodogram
    p_smooth <- LSTS::smooth.periodogram(x, plot=TRUE)
    p_smooth_max_loc <- nnet::which.is.max(p_smooth$smooth.periodogram)
    tau <-2*pi/p_smooth$lambda[p_smooth_max_loc]*Delta
  }

  my_acf <- my_tsa_acf(x, lag.max=GMMlag)#TSA::acf(x,lag.max=GMMlag, plot=F)

  fit_LM_periodic_trawl_foroptim <- function(y){

    alpha <- y[1]
    H <- y[2]


    a<-y[3:(3+m)]
    b<-y[(3+m+1):(3+2*m)]

    lag <- GMMlag
    lss <- 0


    for(i in 1:lag)
    {
      Cor <- kernel_trawl_acf_LM_periodic(Delta*i, alpha, H, a, b, tau)
      lss <- lss+(my_acf$acf[i] - Cor)^2
    }
    lss
  }


  o <- DEoptim::DEoptim(fit_LM_periodic_trawl_foroptim,c(0,1, (-100+numeric(2*m+1))),c(100,10,(100+numeric(2*m+1))),control=DEoptim::DEoptim.control(itermax = 1000, trace = FALSE))

  alpha <- o$optim$bestmem[1]
  H <- o$optim$bestmem[2]
  a<-o$optim$bestmem[3:(3+m)]
  b<-o$optim$bestmem[(3+m+1):(3+2*m)]


  if(plotacf){
    tt <- (1:lags)
    #TSA::acf(x,lag.max=lags,main = "", ylab="ACF", xlab="Lags",drop.lag.0 = TRUE)
    my_tsa_acf(x, lag.max=lags)
    graphics::lines(tt, kernel_trawl_acf_LM_periodic(tt*Delta, alpha, H, a, b, tau), lty =1,col=2, lwd=2)
  }

  return(list("alpha"=alpha, "H"=H, "tau"=tau, "avector"=a, "bvector=", b))

}
