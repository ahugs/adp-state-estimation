source(paste(PATH, 'SEIR.R', sep='/'))

ModelAveraging <- function(y, seir_params, freqs=c(1, 2, 3), niter=500) {
  N = length(y)
  model_fit = list()
  g = list()

  for(i in 1:length(freqs)){
    model_fit[[i]] = 0
    g[[i]] = rep(NA, N)

    for(j in 1:i){
      y_sub = ts(y[seq(j, N, freqs[i])])
      model <- bsts(y_sub, state.specification=AddAr(y=y_sub, lags=1), niter=niter)
      g[[i]][seq(j, N, freqs[i])] = model$state.contributions[which(model$log.likelihood==max(model$log.likelihood)),,]
      model_fit[[i]] = model_fit[[i]] + max(model$log.likelihood)
    }
  }
  model_fits = unlist(model_fit)

  states_df = data.frame(g)
  weighted_states = apply(states_df, MARGIN=1, FUN=function(row){sum(row * exp(model_fits/2))/sum(exp(model_fits/2))})
  seir = SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']], seir_params[['I_0']], seir_params[['R_0']],
                         seir_params[['N']], seir_params[['alpha']], seir_params[['gamma']], seir_params[['beta']],
                         weighted_states)
  return(list(model_fit=model_fit, g=g, states=weighted_states, seir=seir))
}
