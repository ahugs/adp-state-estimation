PATH = '~/github/adp-state-estimation/R/model_averaging'
source(paste(PATH, 'SEIR.R', sep='/'))
source(paste(PATH, 'model_averaging.R', sep='/'))

require(ggplot2)
require(bsts)
require(reshape2)
require(raster)
require(gridExtra)
require(dlm)


set.seed(0)


# Set Initial Parameters
num=1000
pop = 1000
I_0 = 10
S_0 = pop - I_0
E_0 = 0
R_0 = 0

# Set SEIR Parameters (Dukic et al, Tracking epidemics with Google Flu trends data and a state-space SEIR model)
alpha = 1.5/10
gamma = 0.35/10
beta = 0.4/10
sig_g = 0.3/10
sig_y = 0.13/10

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

RunBSTS<- function(y, lags, seir_params, niter=500){
  ar_seir <- list(S=rep(NA, length(y)), E=rep(NA, length(y)), I=rep(NA, length(y)),
                  R=rep(NA, length(y)), g=rep(NA, length(y)))
  tryCatch(
  {
    ar <- bsts(y, AddAutoAr(list(), y=y, lags=lags), niter=niter)
    burn = SuggestBurn(0.1, ar)
    ar_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                               seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                               colMeans(ar$state.contributions[-(1:burn),,])
    )

  },
  error=function(e) {
    print(e)
  }
  )
  return(ar_seir)
}

RunDLM<- function(y, lags, seir_params){
  ar_seir <- list(S=rep(NA, length(y)), E=rep(NA, length(y)), I=rep(NA, length(y)),
                  R=rep(NA, length(y)), g=rep(NA, length(y)))
  tryCatch(
    {
      model <- dlmMLE(y=y,parm=c(rep(0, lags),1,1),
                   build=function(parm) dlmModARMA(ar=parm[1:lags],
                                                   sigma2=parm[lags+1],
                                                   dV=parm[lags+2]))
      states <- dlmSmooth(y=y, mod=dlmModARMA(ar=model$par[1:lags],
                                              sigma2=model$par[lags+1],
                                                  dV=model$par[lags+2]))$s
      if (dim(data.frame(states))[2] > 1) {
        states = states[-1,1]
      }
      else {
        states = states[-1]
      }
      ar_seir <- SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                                 seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                                 states
      )

    },
    error=function(e) {
      print(e)
    }
  )
  return(ar_seir)
}

PlotRun <- function(results) {
  sim_plot_df = data.frame(S=results[['S']],
                           E=results[['E']],
                           I=results[['I']],
                           R=results[['R']])
  sim_plot_df['total'] = apply(sim_plot_df, MARGIN=1, FUN=sum)
  sim_plot_df['t'] = rep(seq(1, nrow(sim_plot_df)), 1)
  sim_plot_df = sim_plot_df[seq(1, nrow(sim_plot_df), nrow(sim_plot_df)/200),]

  sim_plot_df = melt(sim_plot_df, id=c('t'))

  return(ggplot(sim_plot_df, aes(y=value, x=t, col=variable)) + geom_line())
}

run_params_list = list()
alpha_mean = (1+2)/2
alpha_sd = (alpha_mean-1)/qnorm(0.95, 0, 1)#divide by 10

gamma_mean = (0.1 + 0.6)/2
gamma_sd = (gamma_mean-0.1)/qnorm(0.95, 0, 1)

beta_mean = (0.15+0.9)/2
beta_sd = (beta_mean-0.15)/qnorm(0.95, 0, 1)

run_params_list[[1]] = list(mean=list(alpha=alpha_mean, gamma=gamma_mean,beta=beta_mean),
                       sd=list(alpha=alpha_sd, gamma=gamma_sd, beta=beta_sd))
alpha_mean = (0.75+2)/2
alpha_sd = (alpha_mean-0.75)/qnorm(0.95, 0, 1)

gamma_mean = (0.1 + 0.35)/2
gamma_sd = (gamma_mean-0.1)/qnorm(0.95, 0, 1)

beta_mean = (0.1+0.4)/2
beta_sd = (beta_mean-0.1)/qnorm(0.95, 0, 1)

run_params_list[[2]] = list(mean=list(alpha=alpha_mean, gamma=gamma_mean,beta=beta_mean),
                            sd=list(alpha=alpha_sd, gamma=gamma_sd, beta=beta_sd))

# models = list(ma=function(y, seir_params) ModelAveraging(y, seir_params, freqs=c(1,2,3,4))$seir,
#               ar1=function(y, seir_params) RunBSTS(y, 1, seir_params),
#               ar2=function(y, seir_params) RunBSTS(y, 2, seir_params),
#               ar3=function(y, seir_params) RunBSTS(y, 3, seir_params))
models = list(ma=function(y, seir_params) ModelAveraging(y, seir_params, freqs=c(1,2,3,4))$seir,
              ar1=function(y, seir_params) RunDLM(y, 1, seir_params),
              ar2=function(y, seir_params) RunDLM(y, 2, seir_params),
              ar3=function(y, seir_params) RunDLM(y, 3, seir_params))
for (run_params in run_params_list) {
  seir_list = list()
  timer_list = list()
  params_list = list()
  results = lapply(names(models), function(x) list())
  names(results) = names(models)

  pop = 1000
  num = 1000
  I_0 = 10
  S_0 = pop - I_0
  E_0 = 0
  R_0 = 0
  sig_g = 0.2/10
  sig_y = 0.1/10

  alpha_mean = run_params[['mean']][['alpha']]/10
  alpha_sd = run_params[['sd']][['alpha']]/10
  gamma_mean = run_params[['mean']][['gamma']]/10
  gamma_sd = run_params[['sd']][['gamma']]/10
  beta_mean = run_params[['mean']][['beta']]/10
  beta_sd = run_params[['sd']][['beta']]/10
  dir.create(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/",
                     PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                     beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d")))
  dir.create(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/bad_runs",
                     PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                     beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d")))
  dir.create(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/good_runs",
                     PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                     beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d")))
  for(i in 1:200){
    alpha = max(rnorm(1, alpha_mean, alpha_sd), 0)
    gamma = max(rnorm(1, gamma_mean, gamma_sd), 0)
    beta = max(rnorm(1, beta_mean, beta_sd), 0)

    params_list[[i]] = list(alpha=alpha, gamma=gamma, beta=beta)
    seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)
    total = seir[['S']] + seir[['E']] + seir[['I']] + seir[['R']]
    while((sum(total<=pop*0.8, na.rm=TRUE) > 0) | (sum(total>pop*1.2, na.rm=TRUE) > 0) |
          (sum(seir[['I']]<0, na.rm=TRUE) > 0)) {
      badrun_p <- PlotRun(seir)
      tryCatch ( {
        print(badrun_p)
        ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/bad_runs/%.2f_%.2f_%.2f.pdf",
                       PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                       beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"),
                       alpha, gamma, beta))
      }, error=function(e) {}
      )
      alpha = max(rnorm(1, alpha_mean, alpha_sd), 0)
      gamma = max(rnorm(1, gamma_mean, gamma_sd), 0)
      beta = max(rnorm(1, beta_mean, beta_sd), 0)

      params_list[[i]] = list(alpha=alpha, gamma=gamma, beta=beta)
      seir = SEIR(S_0, E_0, I_0, R_0, num, alpha, gamma, beta, sig_g, sig_y)
      total = seir[['S']] + seir[['E']] + seir[['I']] + seir[['R']]
    }
    goodrun_p <- PlotRun(seir)
    print(goodrun_p)
    ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/good_runs/%.2f_%.2f_%.2f.pdf",
                   PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                   beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"),
                   alpha, gamma, beta))

    seir_list[[i]] = seir

    y = seir[['y']][-1]

    seir_params = list(S_0=S_0, E_0=E_0, I_0=I_0, R_0=R_0,
                       N=num, alpha=alpha, gamma=gamma, beta=beta)


    timer_list[[i]] = list()
    for (model in names(models)){
      timer_list[[i]][[model]] = system.time(temp <- models[[model]](y, seir_params))['elapsed']
      results[[model]][[i]] = temp
    }
  }


  ### Plot results summary
  times_df=data.frame(
    matrix(unlist(lapply(names(models),
                         function(name) sapply(timer_list,
                                               function(x) x[[name]]))),
           nrow=length(timer_list), byrow=FALSE))
  colnames(times_df) = names(models)
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
    model_states_list = lapply(results,
                               function(model_results) lapply(model_results,
                                                      FUN = function(x) x[[state]]))
    rmse = lapply(model_states_list, function(model_states) {
      sapply(seq(1, length(true_states_list)),
             FUN = function(i) sqrt(mean((true_states_list[[i]] - model_states[[i]])^2) ))
    })
    df = data.frame(matrix(unlist(rmse), nrow=length(rmse[[1]])))
    colnames(df) = names(models)
    df[!sapply(df,FUN = is.finite)] = NaN

    clip_upper = quantile(df, 0.95, na.rm=TRUE)
    print(clip_upper)
    hist_df = melt(as.data.frame(apply(df, FUN=function(x) clamp(x, lower=-Inf, upper=clip_upper), MARGIN=2)))
    hist_p <- ggplot(hist_df) + geom_histogram(aes(x=value, fill=variable, bins=40)) + facet_wrap(~variable)
    hist_p <- hist_p + ggtitle(sprintf("RMSE Distribution of State %s", state))
    print(hist_p)
    ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/hist_%s.pdf",
                   PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                   beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))


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
    print(winners_p)
    ggsave(sprintf("%s/results/%d_%d_%d_%d_%d_%d_%.2f_%.2f_%.2f_%.2f_%.2f__%s/winner_%s.pdf",
                   PATH, num, pop, I_0, S_0, E_0, R_0, alpha_mean, gamma_mean,
                   beta_mean, sig_g, sig_y, format(as.Date(Sys.time()), "%Y%m%d"), state))
  }
}

