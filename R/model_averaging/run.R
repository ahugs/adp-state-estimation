PATH = '~/github/adp-state-estimation/R/model_averaging'
source(paste(PATH, 'SEIR.R', sep='/'))
source(paste(PATH, 'model_averaging.R', sep='/'))

require(ggplot2)


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

# Run model averaging for 100 SEIR models
model_averaging_results = list()
ar1_results = list()
ar2_results = list()
ar3_results = list()
seir_list = list()

for(i in 1:10){
  seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)
  seir_list[[i]] = seir

  y = seir[['y']][-1]

  model_averaging <- ModelAveraging(y, seir_params, freq=c(1,2,3,4))
  ar1 <- bsts(y, AddAr(y=y, lags=1), niter=500)
  ar1_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                              seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                              ar1$state.contributions[which(max(ar1$log.likelihood) == ar1$log.likelihood),,])
  ar2 <- bsts(y, AddAr(y=y, lags=2), niter=500)
  ar2_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                              seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                              ar2$state.contributions[which(max(ar2$log.likelihood) == ar2$log.likelihood),,])
  ar3 <- bsts(y, AddAr(y=y, lags=3), niter=500)
  ar3_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                              seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                              ar3$state.contributions[which(max(ar3$log.likelihood) == ar3$log.likelihood),,])

  model_averaging_results[[i]] = model_averaging$seir
  ar1_results[[i]] = ar1_seir
  ar2_results[[i]] = ar2_seir
  ar3_results[[i]] = ar3_seir

}

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
  dir.create(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f",
                     PATH, num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
                     beta, sig_g, sig_y))
  mean_df = melt(apply(df, MARGIN=2, function(x) mean(x, na.rm=TRUE)))
  mean_df['model'] = row.names(mean_df)
  mean_p <- ggplot(mean_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  mean_p <- mean_p + ggtitle(sprintf("Mean RMSE of State %s", state))
  print(mean_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f/mean_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
                 beta, sig_g, sig_y, state))

  median_df = melt(apply(df, MARGIN=2, function(x) median(x, na.rm=TRUE)))
  median_df['model'] = row.names(median_df)
  median_p <- ggplot(median_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  median_p <- median_p + ggtitle(sprintf("Median RMSE of State %s", state))
  print(median_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f/median_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
                 beta, sig_g, sig_y, state))

  error_df = melt(apply(df, MARGIN=2, function(x) mean(is.na(x), na.rm=TRUE)))
  error_df['model'] = row.names(error_df)
  error_p <- ggplot(error_df, aes(x=model, y=value, fill=model)) + geom_bar(stat = 'identity', alpha=0.25)
  error_p <- error_p + ggtitle(sprintf("Percent Run Errors for State %s", state))
  print(error_p)
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f/error_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
                 beta, sig_g, sig_y, state))

  df['winner'] = unlist(apply(df, MARGIN=1, FUN=function(row) { if(sum(is.na(row)) == length(row)) {
                                                                    return('None')}
                                                                else {
                                                                  return(colnames(df)[which(row==min(row, na.rm=TRUE))])
                                                                }}))

  winners_p <- ggplot(df, aes(x=winner, fill=winner)) + geom_bar(alpha=0.25)
  winners_p <- winners_p + ggtitle(sprintf("Count of Winning Models (by min RMSE) of State %s", state))
  ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f/winner_%s.pdf",
                 PATH, num, pop, I_0, S_0, E_0, R_0, alpha, gamma,
                 beta, sig_g, sig_y, state))
}





