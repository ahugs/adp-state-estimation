source('R/adp.R')
require(ggplot2)
require(reshape)

set.seed(999)
N = 101
sigv = 1
x = cumsum(rcauchy(N, 0, 0.1))
y = ts(x[-1] + rnorm(N-1, 0, sigv))

transition_func <- function(y_t, y_tm1){
  return(dcauchy(y_t-y_tm1, 0, 0.1))
}

obs_func <- function(s, x) {
  return(dnorm(x, s, sigv, log=TRUE))
}

# Run ADP
start_time = Sys.time()
results = adp(y, transition_func, obs_func, rand_n=3, prep_iter=15, iter=100)
end_time = Sys.time()
print(end_time-start_time)

ar1_params = estimate_params(y)
kf = Kfilter0(length(y), y, 1, mu0=0,
              Sigma0=(ar1_params['sigw', 'estimate']^2/(1-ar1_params['phi', 'estimate']^2)),
              ar1_params['phi', 'estimate'], ar1_params['sigw', 'estimate'],
              ar1_params['sigv', 'estimate'])
results_df = data.frame(kf=c(kf$xf), adp=results[[6]][, ncol(results[[6]])], states=x[-1], t=seq(1, N-1, 1))
results_df = melt(results_df, id=c('t'))
ggplot(data=results_df, aes(x=t, y=value, col=variable)) + geom_line()

mse = c()
for(i in 1:ncol(results[[6]])){
  mse = c(mse, sum((results[[6]][, i] - x[-1])^2)/length(results[[6]][, i]))
}

mse_plot_df = data.frame(adp=mse, kf=rep(sum((c(kf$xf)- x[-1])^2)/length(c(kf$xf)), length(mse)),
                         mean=rep(sum((mean(y) - x[-1])^2)/length(y), length(mse)), iter=seq(1, length(mse)))
mse_plot_df = mse_plot_df[-seq(1, 20),]
mse_plot_df = melt(mse_plot_df, id='iter')
colnames(mse_plot_df) = c('iter', 'model', 'mse')
ggplot(data=mse_plot_df, aes(x=iter, y=mse, col=model))+geom_line()
