PATH = '~/github/adp-state-estimation/R/model_averaging'
source(paste(PATH, 'SEIR.R', sep='/'))
source(paste(PATH, 'model_averaging.R', sep='/'))

require(ggplot2)
require(bsts)
require(reshape2)


set.seed(0)


# Set Initial Parameters
num=100
pop = 1000
I_0 = 10
S_0 = pop - I_0
E_0 = 0
R_0 = 0

# Set SEIR Parameters (Dukic et al, Tracking epidemics with Google Flu trends data and a state-space SEIR model)
alpha = 1.5
gamma = 0.35
beta = 0.4
sig_g = 0.3
sig_y = 0.13

seir_params = list(S_0=S_0, E_0=E_0, I_0=I_0, R_0=R_0, N=num, alpha=alpha, gamma=gamma, beta=beta)

# Simulate SEIR models and plot
sim_n = 9
sim_plot_df = data.frame(S=numeric(0), E=numeric(0), I=numeric(0), R=numeric(0), run=numeric(0))
for(i in 1:sim_n){
  seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)
  sim_plot_df = rbind(sim_plot_df, data.frame(S=seir[['S']], E=seir[['E']], I=seir[['I']], R=seir[['R']], run=i))
}

sim_plot_df['total'] = apply(sim_plot_df, MARGIN=1, FUN=sum)
sim_plot_df['t'] = rep(seq(1, num), sim_n)
sim_plot_df = melt(sim_plot_df, id=c('t', 'run'))

ggplot(sim_plot_df, aes(y=value, x=t, col=variable)) + geom_line() + facet_wrap(~ run, ncol=3)
ggsave(sprintf("%s/%s/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s.pdf",
               PATH, 'simulations', num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
               beta, sig_g, sig_y, format(Sys.time(), "%Y%m%d%H%M")))

# Run model averaging for random SEIR models
model_averaging_results = list()
ar1_results = list()
ar2_results = list()
ar3_results = list()
seir_list = list()
timer_list = list()
params_list = list()


RunBSTS<- function(y, lags, seir_params, niter=500){
  ar_seir <- list(S=rep(NA, length(y)), E=rep(NA, length(y)), I=rep(NA, length(y)),
                  R=rep(NA, length(y)), g=rep(NA, length(y)))
  tryCatch(
  {
    ar <- bsts(y, AddAr(y=y, lags=lags), niter=niter)
    ar_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                               seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                               ar$state.contributions[which(max(ar$log.likelihood) == ar$log.likelihood),,])

  },
  error=function(e) {
    print(e)
  }
  )
  return(ar_seir)
}

run_params = 'b'
for(i in 1:500){
  num = 100
  pop = 1000
  I_0 = 10
  S_0 = pop - I_0
  E_0 = 0
  R_0 = 0
  sig_g = 0.2
  sig_y = 0.1

  if(run_params == 'a') {
    alpha_mean = (1+2)/2
    alpha_sd = (alpha_mean-1)/qnorm(0.95, 0, 1)
    alpha = max(rnorm(1, alpha_mean, alpha_sd), 0)

    gamma_mean = (0.1 + 0.6)/2
    gamma_sd = (gamma_mean-0.1)/qnorm(0.95, 0, 1)
    gamma = max(rnorm(1, gamma_mean, gamma_sd), 0)

    beta_mean = (0.15+0.9)/2
    beta_sd = (beta_mean-0.15)/qnorm(0.95, 0, 1)
    beta = max(rnorm(1, beta_mean, beta_sd), 0)
  }
  if(run_params == 'b') {
    alpha_mean = (0.75+2)/2
    alpha_sd = (alpha_mean-0.75)/qnorm(0.95, 0, 1)
    alpha = max(rnorm(1, alpha_mean, alpha_sd), 0)

    gamma_mean = (0.1 + 0.35)/2
    gamma_sd = (gamma_mean-0.1)/qnorm(0.95, 0, 1)
    gamma = max(rnorm(1, gamma_mean, gamma_sd), 0)

    beta_mean = (0.1+0.4)/2
    beta_sd = (beta_mean-0.1)/qnorm(0.95, 0, 1)
    beta = max(rnorm(1, beta_mean, beta_sd), 0)
  }

  params_list[[i]] = list(alpha=alpha, gamma=gamma, beta=beta)
  seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)
  seir_list[[i]] = seir

  y = seir[['y']][-1]

  seir_params = list(S_0=S_0, E_0=E_0, I_0=I_0, R_0=R_0, N=num, alpha=alpha, gamma=gamma, beta=beta)
  ma_time <- system.time(model_averaging <- ModelAveraging(y, seir_params, freq=c(1,2,3,4)))['elapsed']
  ar1_time <- system.time(ar1_seir <- RunBSTS(y, 1, seir_params))['elapsed']
  ar2_time <- system.time(ar2_seir <- RunBSTS(y, 2, seir_params))['elapsed']
  ar3_time <- system.time(ar3_seir <- RunBSTS(y, 3, seir_params))['elapsed']

  model_averaging_results[[i]] = model_averaging$seir
  ar1_results[[i]] = ar1_seir
  ar2_results[[i]] = ar2_seir
  ar3_results[[i]] = ar3_seir
  timer_list[[i]] = list(ma=ma_time, ar1=ar1_time, ar2=ar2_time, ar3=ar3_time)

}


### Plot results summary
dir.create(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/",
                   PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                   beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d")))
times_df = data.frame(model_averaging=sapply(timer_list, function(x) x[['ma']]),
                      ar1=sapply(timer_list, function(x) x[['ar1']]),
                      ar2=sapply(timer_list, function(x) x[['ar2']]),
                      ar3=sapply(timer_list, function(x) x[['ar3']]))
mean_times_df = melt(apply(times_df, MARGIN=2, function(x) mean(x, na.rm=TRUE)))
mean_times_df['model'] = row.names(mean_times_df)
mean_times_p <- ggplot(mean_times_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
mean_times_p <- mean_times_p + ggtitle(sprintf("Mean Run Time (s)"))
print(mean_times_p)
ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/mean_run_time.pdf",
               PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
               beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d")))


for (state in c('S', 'E', 'I', 'R', 'g')){
  true_states_list = lapply(seir_list, FUN = function(x) x[[state]])
  if(state=='g'){
    true_states_list = lapply(true_states_list, FUN= function(x) x[-1])
  }
  ar1_states_list = lapply(ar1_results, FUN = function(x) x[[state]])
  ar2_states_list = lapply(ar2_results, FUN = function(x) x[[state]])
  ar3_states_list = lapply(ar3_results, FUN = function(x) x[[state]])
  model_averaging_states_list = lapply(model_averaging_results, FUN = function(x) x[[state]])

  rmse_ar1 = sapply(seq(1, length(true_states_list)),
                    FUN = function(i) sqrt(mean((true_states_list[[i]] - ar1_states_list[[i]])^2) ))
  rmse_ar2 = sapply(seq(1, length(true_states_list)),
                    FUN = function(i) sqrt(mean((true_states_list[[i]] - ar2_states_list[[i]])^2) ))
  rmse_ar3 = sapply(seq(1, length(true_states_list)),
                    FUN = function(i) sqrt(mean((true_states_list[[i]] - ar3_states_list[[i]])^2) ))
  rmse_model_averaging = sapply(seq(1, length(true_states_list)),
                    FUN = function(i) sqrt(mean((true_states_list[[i]] - model_averaging_states_list[[i]])^2) ))

  df = data.frame(ar1=rmse_ar1, ar2=rmse_ar2, ar3=rmse_ar3,
                  model_averaging=rmse_model_averaging)
  df[!sapply(df,FUN = is.finite)] = NaN

  mean_df = melt(apply(df, MARGIN=2, function(x) mean(x, na.rm=TRUE)))
  mean_df['model'] = row.names(mean_df)
  mean_p <- ggplot(mean_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  mean_p <- mean_p + ggtitle(sprintf("Mean RMSE of State %s", state))
  print(mean_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/mean_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                 beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))

  median_df = melt(apply(df, MARGIN=2, function(x) median(x, na.rm=TRUE)))
  median_df['model'] = row.names(median_df)
  median_p <- ggplot(median_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  median_p <- median_p + ggtitle(sprintf("Median RMSE of State %s", state))
  print(median_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/median_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                 beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))

  error_df = melt(apply(df, MARGIN=2, function(x) mean(is.na(x), na.rm=TRUE)))
  error_df['model'] = row.names(error_df)
  error_p <- ggplot(error_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  error_p <- error_p + ggtitle(sprintf("Percent Run Errors for State %s", state))
  print(error_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/error_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                 beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))

  winner_df = data.frame(winner=unlist(apply(df, MARGIN=1,
                                             FUN=function(row) { if(sum(is.na(row)) == length(row)) {
                                                                    return('None')}
                                                                else {
                                                                  return(colnames(df)[which(row==min(row, na.rm=TRUE))])
                                                                }})))

  winners_p <- ggplot(winner_df, aes(x=winner, fill=winner)) + geom_bar(alpha=0.25)
  winners_p <- winners_p + ggtitle(sprintf("Count of Winning Models (by min RMSE) of State %s", state))
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/winner_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                 beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))
}


