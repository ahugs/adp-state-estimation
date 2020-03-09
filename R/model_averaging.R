devtools::install_github("nickpoison/astsa")
require(astsa)
require(penalized)
require(ggplot2)
require(reshape)
require(mvtnorm)

# Generate Data
set.seed(999)
num = 1000
N = num+1
phi = 0.7
sigw = 1
sigv = 2
x = arima.sim(n=N, list(ar = c(phi, 0.2), sd=sigw))
y = ts(x[-1] + rnorm(num, 0, sigv))


SEIR <- function(S_0, E_0, I_0, R_0, Tn, alpha, gamma, beta, sig_g, sig_y){
  pop_N = S_0 + E_0 + I_0 + R_0
  S = rep(NA, Tn); S[1] = S_0
  E = rep(NA, Tn); E[1] = E_0
  I = rep(NA, Tn); I[1] = I_0
  R = rep(NA, Tn); R[1] = R_0
  g = rep(NA, Tn)
  y = rep(NA, Tn)
  
  for(i in 2:Tn){
    g[i] = alpha * E[i-1]/I[i-1] - gamma + rnorm(1, 0, sig_g)
    y[i] = g[i] + rnorm(1, 0, sig_y)
    I[i] = (1 + g[i])*I[i-1]
    S[i] = S[i-1] - beta*S[i-1]*I[i-1]/pop_N
    E[i] = E[i-1] + beta*S[i-1]*I[i-1]/pop_N - alpha*E[i-1]
    R[i] = R[i-1] + gamma*I[i-1]
  }
  return(list(S, E, I, R, g, y))
}

seir_states = SEIR(800, 0, 200, 0, 100, 0.5, 0.1, 1, sqrt(1.1), sqrt(1.1))

seir_df = data.frame(S=seir_states[[1]], E=seir_states[[2]], I=seir_states[[3]], R=seir_states[[4]],
                     t=seq(1, 100))
seir_df = melt(seir_df, id='t')
ggplot(seir_df, aes(y=value, x=t, col=variable)) + geom_line()

g_df =data.frame(g=seir_states[[5]], y=seir_states[[6]], t=seq(1,100))
g_df = melt(g_df, id='t')
ggplot(g_df, aes(y=value, x=t, col=variable))+geom_line()

model_averaging <- function(y, freqs=c(1,2,3)) {
  num=length(y)
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
      #max_llike = -params_estimate[[2]]
      bic[[i]][[j]] = posterior_prob(params_estimate, length(y_freq))#bic(max_llike, length(y_freq), 3)
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

posterior_prob <- function(params_estimate, n) {
  sigma = params_estimate[[2]]
  if(!isSymmetric(sigma)){
    sigma = apply(params_estimate[[2]], MARGIN=1:2, FUN=function(x) round(x, 6))
  }
  
  rparams_df = data.frame(rmvnorm(n=1000, mean=params_estimate[[1]][, 'estimate'], sigma=sigma))
  colnames(rparams_df) = rownames(params_estimate[[1]])
  Linn=function(para){
       phi = para[1]; sigw = para[2]; sigv = para[3]
       Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
       kf = Kfilter0(length(y),y,1,mu0=0,Sigma0,phi,sigw,sigv)
       return(kf$like)
   }
  return(mean(apply(rparams_df, MARGIN=1, FUN=Linn)/n))
  
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

  # Function to evaluate the likelihood
  Linn=function(para){
    phi = para[1]; sigw = para[2]; sigv = para[3]
    Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
    kf = Kfilter0(length(y),y,1,mu0=0,Sigma0,phi,sigw,sigv)
    return(kf$like)
  }
  
  # Estimation
  (est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE))
  inv_hessian = solve(est$hessian)
  SE = sqrt(diag(inv_hessian))
  return(list(data.frame(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)), inv_hessian))
}

results = model_averaging(ts(seir_states[[6]][2:100]))

y = seir_states[[6]][2:100]
x = seir_states[[5]][1:100]

plot_df = data.frame(model_averaging=results[[4]], true_states=x[2:length(x)], ksmooth=results[[3]][[1]][[1]],
                     t=seq(1, length(y)))
mse_model_averaging=sum((plot_df[, 'model_averaging'] - plot_df[,'true_states'])^2)/length(plot_df)
mse_ksmooth=sum((plot_df[, 'ksmooth'] - plot_df[,'true_states'])^2)/length(plot_df)

plot_df = plot_df[1:100,]
plot_df = melt(plot_df, id='t')
plot_df[, 'variable'] <- sapply(plot_df[, 'variable'], as.character)

plot_df[plot_df[, 'variable'] == 'model_averaging', 'variable'] = paste('model_averaging \n (mse=', round(mse_model_averaging, 2), ')', sep='')
plot_df[plot_df[, 'variable'] == 'ksmooth', 'variable'] = paste('ksmooth \n(mse=', round(mse_ksmooth, 2), ')', sep='')

ggplot(data=plot_df, aes(x=t, y=value, col=variable)) + geom_line()# + ggtitle('AR(2) Model (ar = c(0.7, 0.2))')




