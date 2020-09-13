source('R/adp.R')
require(ggplot2)
require(reshape)

# Simulate data
set.seed(999)
N = 101
phi = 0.7
sigw = 1
sigv = 1
x = arima.sim(n=N, list(ar = c(phi), sd=sigw))
y = ts(x[-1] + rnorm(N-1, 0, sigv))

# Set up inputs

transition_func <- function(y_t, y_tm1){
  return(dnorm(y_t, phi*y_tm1, sigw, log=TRUE))
}

obs_func <- function(s, x) {
  return(dnorm(x, s, sigv, log=TRUE))
}

# Run ADP
results = adp(y, transition_func, obs_func, freqs=c(1, 2))

results_df = data.frame(kf=results[[2]][[1]], adp=results[[6]][, ncol(results[[6]])], states=x[-1], t=seq(1, N-1, 1))
results_df = melt(results_df, id=c('t'))
ggplot(data=results_df, aes(x=t, y=value, col=variable)) + geom_line()

mse = c()
for(i in 1:ncol(results[[6]])){
  mse = c(mse, sum((results[[6]][, i] - x[-1])^2)/length(results[[6]][, i]))
}

mse_plot_df = data.frame(adp=mse, kf=rep(sum((results[[2]][[1]]- x[-1])^2)/length(results[[2]][[1]]), length(mse)),
                         mean=rep(sum((mean(y) - x[-1])^2)/length(y), length(mse)), iter=seq(1, length(mse)))
mse_plot_df = mse_plot_df[-seq(1, 20),]
mse_plot_df = melt(mse_plot_df, id='iter')
colnames(mse_plot_df) = c('iter', 'model', 'mse')
ggplot(data=mse_plot_df, aes(x=iter, y=mse, col=model))+geom_line()
