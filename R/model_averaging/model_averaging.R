source(paste(PATH, 'SEIR.R', sep='/'))
require(parallel)

ModelAveraging <- function(y, seir_params, freqs=c(1, 2, 3), niter=500) {
  N = length(y)
  model_results = mclapply(freqs,
                           FUN=function(freq) FitSingleFrequency(freq=freq, y=y, N=N, niter=niter),
                           mc.cores=4)

  model_fits = lapply(model_results, FUN=function(x) x[['model_fit']])
  g = lapply(model_results, FUN=function(x) x[['g']])

  model_fits = unlist(model_fits)

  states_df = data.frame(g)
  weighted_states = apply(states_df, MARGIN=1, FUN=function(row){sum(row * exp(model_fits/2))/sum(exp(model_fits/2))})
  seir = SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                         seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                         weighted_states)
  return(list(model_fit=model_fits, g=g, states=weighted_states, seir=seir))
}

FitSingleFrequency<- function(freq, y, N, niter){
  model_fit = 0
  g = rep(NA, N)

  tryCatch(
    for(j in 1:freq){
      y_sub = ts(y[seq(j, N, freq)])
      model <- bsts(y_sub, state.specification=AddAr(y=y_sub, lags=1), niter=niter)
      g[seq(j, N, freq)] = model$state.contributions[which(model$log.likelihood==max(model$log.likelihood)),,]
      model_fit = model_fit + max(model$log.likelihood)
    },
    error=function(e) {}
  )
  return(list(g=g, model_fit=model_fit))
}



