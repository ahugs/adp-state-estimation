devtools::install_github("nickpoison/astsa")
require(astsa)
require(penalized)
require(ggplot2)
require(reshape)

# Generate Data
set.seed(999)
num = 1000
N = num+1
phi = 0.7
sigw = 1
sigv = 2
x = arima.sim(n=N, list(ar = c(phi, 0.2), sd=sigw))
y = ts(x[-1] + rnorm(num, 0, sigv))


model_averaging <- function(freqs=c(1,2,3)) {
  ar1_params = list()
  bic = list()
  states = list()
  
  for(i in 1:length(freqs)){
    ar1_params[[i]] = list()
    bic[[i]] = list()
    states[[i]] = list()
    for(j in 1:i){
      y_freq = ts(y[seq(j, num, freqs[i])])
      params_estimate = estimate_params(y_freq)
      ar1_params[[i]][[j]] = params_estimate[[1]]
      max_llike = -params_estimate[[2]]
      bic[[i]][[j]] = bic(max_llike, length(y_freq), 3)
      states[[i]][[j]] = rep(NA, length(y))
      states[[i]][[j]][seq(j, num, freqs[i])] = 
        c(Ksmooth0(length(y_freq), y_freq, A=1, mu0=0, Sigma0=1, ar1_params[[i]][[j]]['phi', 'estimate'], 
               ar1_params[[i]][[j]]['sigw', 'estimate'], ar1_params[[i]][[j]]['sigv', 'estimate'])$xs)

    }
  }
  bics = unlist(bic)
  sigvs = unlist(sapply(ar1_params, function(l) {unlist(sapply(l, function(l2) l2['sigv', 'estimate']))}))
  sigws = unlist(sapply(ar1_params, function(l) {unlist(sapply(l, function(l2) l2['sigw', 'estimate']))}))
  valid = (sigvs > 0) & (sigw > 0)
  states_df = data.frame(states)
  weighted_states = apply(states_df, MARGIN=1, FUN=function(row){sum(row * exp(-bics/2) * valid * !is.na(row), na.rm=T)/
      sum(exp(-bics/2) * valid * !is.na(row))})
  
  return (list(ar1_params, bic, states, weighted_states))
}

bic <- function(max_llike, n, k){
  #return(log(n)*k - 2*max_llike)
  return(-max_llike/n)
}

estimate_params <- function(y) {
  u = ts.intersect(y, lag(y,-1), lag(y,-2))
  varu = var(u)
  coru = cor(u)
  phi = min(coru[1,3]/coru[1,2],0.99)
  q = max((1-phi^2)*varu[1,2]/phi,0.1)
  r = max(varu[1,1] - q/(1-phi^2), 0.1)
  # print(paste('varu ', varu[1,1]))
  # print(paste('rrac ', q/(1-phi^2)))
  (init.par = c(phi, sqrt(q), sqrt(r)))
  print(init.par)
  
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
  return(list(data.frame(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)), est$value))
}

results = model_averaging()

plot_df = data.frame(model_averaging=results[[4]], true_states=x[2:length(x)], ksmooth=results[[3]][[1]][[1]],
                     t=seq(1, length(y)))
mse_model_averaging=sum((plot_df[, 'model_averaging'] - plot_df[,'true_states'])^2)/length(plot_df)
mse_ksmooth=sum((plot_df[, 'ksmooth'] - plot_df[,'true_states'])^2)/length(plot_df)

plot_df = plot_df[1:100,]
plot_df = melt(plot_df, id='t')
plot_df[, 'variable'] <- sapply(plot_df[, 'variable'], as.character)

plot_df[plot_df[, 'variable'] == 'model_averaging', 'variable'] = paste('model_averaging \n (mse=', round(mse_model_averaging, 2), ')', sep='')
plot_df[plot_df[, 'variable'] == 'ksmooth', 'variable'] = paste('ksmooth \n(mse=', round(mse_ksmooth, 2), ')', sep='')

ggplot(data=plot_df, aes(x=t, y=value, col=variable)) + geom_line() + ggtitle('AR(2) Model (ar = c(0.7, 0.2))')




