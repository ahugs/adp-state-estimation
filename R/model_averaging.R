devtools::install_github("nickpoison/astsa")
require(astsa)
require(ggplot2)
require(reshape)
require(mvtnorm)
# Sparse AR(1) model
# try the linear interpolation of the states
# try with a stable set of parameters for the simulations (look at the paper for equations on seeting initial params)
# try solving for all the states after getting the kalman smoothed g states
# Correlatoin parameter changing over time? Sequentially for all data up to and include current, or slices of the data

# Try varying ratio of variances
# filter outlier RMSE runs
require(bsts)
require(invgamma)
require(zoo)

# Generate Data
set.seed(999)
num = 400
N = num+1

# AR Model
phi1 = 0.7
phi2 = 0.2
sigw = 1
sigv = 2
x = arima.sim(n=N, list(ar = c(phi1), sd=sigw))
y = ts(x[-1] + rnorm(num, 0, sigv))

# SEIR Model
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
  return(list(S=S, E=E, I=I, R=R, g=g, y=y))
}

num=100
pop = 1000
I_0 = 10
S_0 = pop - I_0
E_0 = 0
R_0 = 0
# alpha = max(rnorm(1, 0.5, 0.5),0)
# gamma = max(rnorm(1, 0.5, 0.5),0)
# beta = max(rnorm(1, 0.5, 0.5),0)
# sig_g = rinvgamma(1, 1.1, 0.005)
# sig_y = rinvgamma(1, 1.1, 0.05)
alpha = 1.5
gamma = 0.35
beta = 0.4
sig_g = 0.3
sig_y = 0.13
seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)


seir_df = data.frame(S=seir[['S']], E=seir[['E']], I=seir[['I']], R=seir[['R']])
seir_df['total'] = apply(seir_df, MARGIN=1, FUN=sum)
seir_df['t'] = seq(1, num)
seir_df = melt(seir_df, id='t')
ggplot(seir_df, aes(y=value, x=t, col=variable)) + geom_line()

g_df =data.frame(g=seir[['g']], y=seir[['y']], t=seq(1,num))
g_df = melt(g_df, id='t')
ggplot(g_df, aes(y=value, x=t, col=variable))+geom_line()

model_averaging <- function(y, init_vals, freqs=c(1,2,3)) {
  num=length(y)
  ar1_params = list()
  model_fit = list()
  g = list()
  for(i in 1:length(freqs)){
    ar1_params[[i]] = list()
    model_fit[[i]] = list()
    g[[i]] = list()

    for(j in 1:i){
      y_freq = ts(y[seq(j, num, freqs[i])])
      params_estimate = estimate_params(y_freq)
      ar1_params[[i]][[j]] = params_estimate[['mu']]
      model_fit[[i]][[j]] = posterior_predictive(params_estimate, y_freq)
      g[[i]][[j]] = rep(NA, length(y))
      g[[i]][[j]][seq(j, num, freqs[i])] =
        c(Ksmooth0(length(y_freq), y_freq, A=1, mu0=0, Sigma0=1, ar1_params[[i]][[j]]['phi', 'estimate'],
               ar1_params[[i]][[j]]['sigw', 'estimate'], ar1_params[[i]][[j]]['sigv', 'estimate'])$xs)

    }
    # for(j in 1:length(y)){
    #   seg_len = max(j, 30)
    #   y_freq = ts(y[seq(1, seg_len, freqs[i])])
    #   params_estimate = estimate_params(y_freq)
    #   ar1_params[[i]][[j]] = params_estimate[['mu']]
    #   model_fit[[i]][[j]] = posterior_predictive(params_estimate, y_freq)
    #   states = rep(NA, seg_len)
    #   states[seq(1, seg_len, freqs[i])] =
    #       c(Ksmooth0(length(y_freq), y_freq, A=1, mu0=0, Sigma0=1, ar1_params[[i]][[j]]['phi', 'estimate'],
    #                  ar1_params[[i]][[j]]['sigw', 'estimate'], ar1_params[[i]][[j]]['sigv', 'estimate'])$xs)
    #   states = na.locf(na.locf(na.approx(states, na.rm=FALSE)), fromLast = TRUE)
    #   g[[i]][[j]] = states[j]
    #
    # }

  }
  model_fits = unlist(model_fit)

  states_df = data.frame(g)
  weighted_states = apply(states_df, MARGIN=1, FUN=function(row){sum(row * exp(-model_fits*2) * !is.na(row), na.rm=T)/
      sum(exp(-model_fits*2)  * !is.na(row))})

  return (list(ar1_params=ar1_params, model_fit=model_fit, states=g, weighted_states=weighted_states))
}


posterior_predictive <- function(params_estimate, y) {
  n = length(y)
  sigma = params_estimate[['sigma']]
  if(!isSymmetric(sigma)){
    sigma = apply(params_estimate[[2]], MARGIN=1:2, FUN=function(x) round(x, 6))
  }

  rparams_df = data.frame(rmvnorm(n=1000, mean=params_estimate[['mu']][, 'estimate'], sigma=sigma))
  colnames(rparams_df) = rownames(params_estimate[['mu']])
  Linn=function(para){
       phi = para[1]; sigw = para[2]; sigv = para[3]
       Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
       kf = Kfilter0(length(y),y,1,mu0=0,Sigma0,phi,sigw,sigv)
       return(kf$like)
   }
  return(mean(apply(rparams_df, MARGIN=1, FUN=Linn)/n))

}

estimate_params <- function(y) {
  u = ts.intersect(y, stats::lag(y,-1), stats::lag(y,-2))
  varu = var(u)
  coru = cor(u)
  phi = min(coru[1,3]/coru[1,2],0.99)
  q = max((1-phi^2)*varu[1,2]/phi,0.1)
  r = max(varu[1,1] - q/(1-phi^2), 0.1)
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
  return(list(mu=data.frame(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)),
              sigma=inv_hessian))
}

results = model_averaging(ts(seir[['y']][2:num]))


mse_model_averaging=sum((plot_df[, 'model_averaging'] - plot_df[,'true_states'])^2)/length(plot_df)
mse_ksmooth=sum((plot_df[, 'ksmooth'] - plot_df[,'true_states'])^2)/length(plot_df)
mse_ar3 = sum((ar3_model$state.contributions[niter,,] - plot_df[,'true_states'])^2)/length(plot_df)
plot_df = melt(plot_df, id='t')
plot_df[, 'variable'] <- sapply(plot_df[, 'variable'], as.character)

plot_df[plot_df[, 'variable'] == 'model_averaging', 'variable'] = paste('model_averaging \n (mse=', round(mse_model_averaging, 2), ')', sep='')
plot_df[plot_df[, 'variable'] == 'ksmooth', 'variable'] = paste('ksmooth \n(mse=', round(mse_ksmooth, 2), ')', sep='')
plot_df[plot_df[, 'variable'] == 'ar3', 'variable'] = paste('ar3 \n(mse=', round(mse_ar3, 2), ')', sep='')

ggplot(data=plot_df, aes(x=t, y=value, col=variable)) + geom_line()


# Compare model performance across multiple SEIR models
mse_model_averaging_list = list()
mse_ksmooth_list = list()
mse_ar3_list = list()
params_list = list()

for(i in 1:100){
  pop = 1000
  I_0 = 10
  S_0 = pop - I_0
  E_0 = 0
  R_0 = 0
  # alpha = alpha#max(rnorm(1, 0.5, 0.5),0)
  # gamma = gamma#max(rnorm(1, 0.5, 0.5),0)
  # beta = beta#max(rnorm(1, 0.5, 0.5),0)
  # sig_g = sig_g#rinvgamma(1, 1.1, 0.005)
  # sig_y = sig_y#rinvgamma(1, 1.1, 0.05)
  alpha = 1.5
  gamma = 0.35
  beta = 0.4
  sig_g = 0.3
  sig_y = 0.13
  # alpha = 1.5/2
  # gamma = 0.3/2
  # beta = 0.8/2
  # sig_g = 0.3/2
  # sig_y = 0.13/2


  params_list[[i]]= list(alpha=alpha, gamma=gamma, beta=beta, sig_g=sig_g, sig_y=sig_y)
  seir_states = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)

  y = seir_states[['y']][2:num]
  x = seir_states[['g']][1:num]


  results <- tryCatch(
    model_averaging(ts(y), freqs=c(1,2,3)),
    error=function(e) print(e)
  )

  par = estimate_params(ts(y))[['mu']]
  ksmooth_states = tryCatch(
    c(Ksmooth0(length(y), y, A=1, mu0=0, Sigma0=1, par['phi', 'estimate'],
               par['sigw', 'estimate'],par['sigv', 'estimate'])$xs),
    error=function(e) print(e))

  ar3_model<-tryCatch(bsts(y, AddAr(y=y, lags=3), niter=500),
                      error=function(e) e)


  if(!inherits(results, "error")){
     results_df = data.frame(model_averaging=results[['weighted_states']], true_states=x[2:length(x)])
     mse_model_averaging_list[[i]]=sum((results_df[, 'model_averaging'] - results_df[,'true_states'])^2)/length(results_df)

  } else {
    mse_model_averaging_list[[i]] = NA
  }

  if(!inherits(ksmooth_states, "error")){
      results_df = data.frame(true_states=x[2:length(x)], ksmooth=ksmooth_states)
      mse_ksmooth_list[[i]]=sum((results_df[, 'ksmooth'] - results_df[,'true_states'])^2)/length(results_df)

  } else {
      mse_ksmooth_list[[i]] = NA
  }

  if(!inherits(ar3_model, "error")){
    mse_ar3_list[[i]] = sum((ar3_model$state.contributions[niter,,] - results_df[,'true_states'])^2)/length(results_df)
  }
  else {
    mse_ar3_list[[i]] = NA
  }
}


df = data.frame(ksmooth=unlist(mse_ksmooth_list), model_average=unlist(mse_model_averaging_list),
                ar3=unlist(mse_ar3_list))
mean_df = melt(apply(df, MARGIN=2, function(x) mean(x, na.rm=TRUE)))
mean_df['model'] = row.names(mean_df)
ggplot(mean_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)

ggplot(melt(df), aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ xlim(0, 2)
ggplot(melt(df) ,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25) + xlim(0, 2)

error_df = melt(apply(df, MARGIN=2, function(x) mean(is.na(x), na.rm=TRUE)))
error_df['model'] = row.names(error_df)
ggplot(error_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)


df = df[rowSums(is.na(df)) < 3, ]
df['winner'] = unlist(apply(df, MARGIN=1, FUN=function(row) colnames(df)[which(row==min(row, na.rm=TRUE))]))
ggplot(df, aes(x=winner, fill=winner)) + geom_bar(alpha=0.25)


