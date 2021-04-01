PATH = '~/github/adp-state-estimation/R/model_averaging'
source(paste(PATH, 'SEIR.R', sep='/'))
require(parallel)
require(dlm)

WeightFunc <- function(row, model_fits, N) {
  sum(row * exp(model_fits/N))/sum(exp(model_fits/N))
}

ModelAveraging <- function(y, seir_params, freqs=c(1, 2, 3), niter=500, weight_func=WeightFunc) {
  N = length(y)
  model_results = mclapply(freqs,
                           FUN=function(freq) FitSingleFrequency(freq=freq, y=y, N=N, niter=niter),
                           mc.cores=4)

  model_fits = lapply(model_results, FUN=function(x) x[['model_fit']])
  model_fits = unlist(model_fits)

  g = lapply(model_results, FUN=function(x) x[['g']])
  seir = lapply(g, function(x) SEIRFixedStates(seir_params[['S_0']], seir_params[['E_0']],
                                               seir_params[['I_0']], seir_params[['R_0']],
                                               seir_params[['N']], seir_params[['alpha']],
                                               seir_params[['gamma']], seir_params[['beta']],
                                               x))
  state_names = names(seir[[1]])
  seir = lapply(state_names, function(s) lapply(seir, FUN = function(x) x[[s]]))
  names(seir) = state_names
  states = lapply(seir, function(x) as.data.frame(matrix(unlist(x),nrow=length(x[[1]]),byrow=FALSE)))

  weighted_states = lapply(states,
                           function(x) apply(x, MARGIN=1,
                                             FUN=function(row) weight_func(row, model_fits, N)))
  return(list(model_fit=model_fits, g=g, states=weighted_states[['g']], seir=weighted_states))
}

# FitSingleFrequency<- function(freq, y, N, niter){
#   model_fit = 0
#   g = rep(NA, N)
#
#   tryCatch(
#     for(j in 1:freq){
#       y_sub = ts(y[seq(j, N, freq)])
#       model <- bsts(y_sub, AddAutoAr(list(), y=y_sub, lags=1), niter=500)
#       burn = SuggestBurn(0.1, model)
#       g[seq(j, N, freq)] = colMeans(model$state.contributions[-(1:burn),,])
#       model_fit = model_fit - mean(model$log.likelihood[-(1:burn)])
#     },
#     error=function(e) {print(e)}
#   )
#   return(list(g=g, model_fit=-(freq*log(length(y)) + 2*model_fit)))
# }


FitSingleFrequency<- function(freq, y, N, niter){
  model_fit = 0
  g = rep(NA, N)

  tryCatch(
    for(j in 1:freq){
      y_sub = ts(y[seq(j, N, freq)])
      model <- dlmMLE(y=y_sub,parm=c(0,1,1),
                      build=function(parm) dlmModARMA(ar=parm[1],
                                                      sigma2=parm[2],
                                                      dV=parm[3]))
      states <- dlmSmooth(y=y_sub, mod=dlmModARMA(ar=model$par[1],
                                              sigma2=model$par[2],
                                              dV=model$par[3]))$s[-1]
      g[seq(j, N, freq)] = states
      model_fit = model_fit + model$value
    },
    error=function(e) {print(e)}
  )
  return(list(g=g, model_fit=-(freq*log(length(y)) + 2*model_fit)))
}



