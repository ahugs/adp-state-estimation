devtools::install_github("nickpoison/astsa")
require(astsa)
require(penalized)
require(ggplot2)
require(reshape)

# Generate Data
set.seed(999)
num = 100
N = num+1
phi = 0.7
sigma_epsilon = 1
sigma_alpha = 1
x = arima.sim(n=N, list(ar = c(phi), sd=sigma_epsilon))
y = ts(x[-1] + rnorm(num, 0, sigma_alpha))

# set up inputs
transition_func <- function(y_t, y_tm1){
  return(dnorm(y_t, phi*y_tm1, sigma_epsilon, log=TRUE))
}
obs_func <- function(s, x) {
  return(dnorm(x, s, sigma_alpha, log=TRUE))
}

adp <- function(y, transition_func, obs_func, freqs=c(1,2,3), iter=100){
  
  # estimate AR(1) parameters and Kalman filter for different frequencies of data
  ar1_params = list()
  xs = list()
  Ps = list()
  
  for(i in 1:length(freqs)){
    y_freq = ts(y[seq(1, num, freqs[i])])
    ar1_params[[i]] = estimate_params(y_freq)
    kf = Kfilter0(length(y_freq),y_freq,1,mu0=0,
                  Sigma0=(ar1_params[[i]]['sigw', 'estimate']^2/(1-ar1_params[[i]]['phi', 'estimate']^2)),
                  ar1_params[[i]]['phi', 'estimate'], ar1_params[[i]]['sigw', 'estimate'],
                  ar1_params[[i]]['sigv', 'estimate'])
    
    xs[[i]] = rep(kf$xp, 1, each=freqs[i])[1:num]
    Ps[[i]] = rep(kf$Pp, 1, each=freqs[i])[1:num]

  }

    # initialize approximation with no data
  df = list(data.frame(v=numeric(0), ar1=numeric(0),  ar2=numeric(0), ar3=numeric(0)))
  V = rep(df, length(y))

  coefs = list(data.frame(ar1=numeric(0),  ar2=numeric(0), ar3=numeric(0)))
  coefs = rep(coefs, length(y))

  S = data.frame(matrix(0, length(y), iter))
  for(i in 1:iter){
    
    # Select random final state
    S[length(y), i] = rnorm(1, 0, 2)

    for(j in length(y):2){
      
      if(nrow(V[[j-1]]) < 10) {
        coefs[[j]][i, ] = c(NA, NA, NA)
        max_func <- function(s){
          -transition_func(S[j, i], s) - obs_func(S[j, i], y[j])
        }
      } else {
        coefs[[j]][i, ] = glm(v~ar1+ar2+ar3-1, family=quasibinomial, data=V[[j]])$coefficients

        max_func <- function(s){

          theta = coefs[[j]][i, 1] * dnorm(s, xs[[1]][j-1], sqrt(Ps[[1]][j-1])) +
                  coefs[[j]][i, 2] * dnorm(s, xs[[2]][j-1], sqrt(Ps[[2]][j-1])) +
                  coefs[[j]][i, 3] * dnorm(s, xs[[3]][j-1], sqrt(Ps[[3]][j-1]))
          -transition_func(S[j, i], s) - obs_func(S[j, i], y[j]) -log(exp(theta)/(1+exp(theta)))
        }
      }
      opt = optim(S[j, i], max_func, method='Brent', lower=-10, upper=10)
      S[j - 1, i] = opt$par
    
      V[[j]] = rbind(V[[j]], data.frame(v=exp(-opt$value),
                                        ar1=dnorm(S[j, i], xs[[1]][j-1], sqrt(Ps[[1]][j-1])),
                                        ar2=dnorm(S[j, i], xs[[2]][j-1], sqrt(Ps[[2]][j-1])),
                                        ar3=dnorm(S[j, i], xs[[3]][j-1], sqrt(Ps[[3]][j-1]))))
      
    }
  }
  return(list(V, xs, Ps, S, coefs))
}

# Initial Estimates
estimate_params <- function(y) {
  u = ts.intersect(y, lag(y,-1), lag(y,-2))
  varu = var(u)
  coru = cor(u)
  phi = min(coru[1,3]/coru[1,2],0.99)
  q = (1-phi^2)*varu[1,2]/phi
  r = max(varu[1,1] - q/(1-phi^2), 0.1)
  # print(paste('varu ', varu[1,1]))
  # print(paste('rrac ', q/(1-phi^2)))
  (init.par = c(phi, sqrt(q), sqrt(r)))

  # Function to evaluate the likelihood
  Linn=function(para){
    phi = para[1]; sigw = para[2]; sigv = para[3]
    Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
    kf = Kfilter0(length(y),y,1,mu0=0,Sigma0,phi,sigw,sigv)
    return(kf$like)
  }

  # Estimation
  # print(init.par)
  (est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE))
  SE = sqrt(diag(solve(est$hessian)))
  return(data.frame(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)))
}


results = adp(y, transition_func, obs_func)

results_df = data.frame(kf=results[[2]][[1]], adp=results[[4]]$X100, states=x[-1], t=seq(1, num, 1))
results_df = melt(results_df, id=c('t'))
ggplot(data=results_df, aes(x=t, y=value, col=variable)) + geom_line()

