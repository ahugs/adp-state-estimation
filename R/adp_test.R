source('R/adp.R')

sample_states <- function(xtt, Ptt, xttm1, Pttm1, x, N, log=F){
  sample = rep(NA, N)
  m  = rep(NA, N)
  V = rep(NA, N)
  m[N] = xtt[N]
  V[N] = Ptt[N]
  sample[N] = x#rnorm(1, m[N], sqrt(V[N]))

  llike = dnorm(sample[N], m[N], sqrt(V[N]), log=log)
  if(N-1 > 1){
    for(t in (N-1):1){
      m[t] = xtt[t] + (Ptt[t] * phi * 1/Pttm1[t+1]) * (sample[t+1] - xttm1[t+1])
      V[t] = Ptt[t] - (Ptt[t] * phi * 1/Pttm1[t+1]) ^ 2 * Pttm1[t+1]
      sample[t] = m[t]#rnorm(1, m[t], sqrt(V[t]))
      if(log){
        llike = llike + dnorm(sample[t], m[t], sqrt(V[t]), log=T)
      } else {
        llike = llike * dnorm(sample[t], m[t], sqrt(V[t]))
      }
    }
  }
  return(list(sample, m, V, llike))
}

prev_states_llike <- function(xtt, Ptt, xttm1, Pttm1, s, N, log=T, x){
  sample = rep(NA, N)
  m  = rep(NA, N)
  V = rep(NA, N)
  m[N] = xtt[N]
  V[N] = Ptt[N]
  #sample[N] = x#rnorm(1, m[N], sqrt(V[N]))
  sample[1:(N-1)]=s
  sample[N] = x
  llike = dnorm(sample[N], m[N], sqrt(V[N]), log=log)
  for(t in (N-1):1){
    m[t] = xtt[t] + (Ptt[t] * phi * 1/Pttm1[t+1]) * (sample[t+1] - xttm1[t+1])
    V[t] = Ptt[t] - (Ptt[t] * phi * 1/Pttm1[t+1]) ^ 2 * Pttm1[t+1]
    #sample[t] = rnorm(1, m[t], sqrt(V[t]))
    llike = llike + dnorm(sample[t], m[t], sqrt(V[t]), log=log)
  }
  return(llike)
}

max_prev_states_llike <- function(x){
  return (optim(par=xttm1[1:(length(xttm1)-2)], fn= function(s) -prev_states_llike(xtt, Ptt, xttm1,
                                                                                   Pttm1, s, 99, log=T, x),
                method='L-BFGS-B')$value)
}

kf = Kfilter0(length(y), y, 1, mu0=0, Sigma0=(sigw^2/(1-phi^2)), phi, sigw,sigv)
xtt = c(kf$xf); Ptt = c(kf$Pf); xttm1 = c(kf$xp); Pttm1 = c(kf$Pp)
x100 = 2
opt = optim(par=c(20), fn=function(s) max_prev_states_llike(s) - transition_func(x100, s), method='Brent', lower=-10, upper=10)
x99 = opt$par

sample_states(xtt, Ptt, xttm1, Pttm1, x99, 99, log=T)[[4]] + obs_likelihood(phi, sigv, sigw, y, log=T)[99] +
  transition_func(x100, x99) + obs_func(x100, y[100])

sample_states(xtt, Ptt, xttm1, Pttm1, x100, 100, log=T)[[4]] + obs_likelihood(phi, sigv, sigw, y, log=T)[100]

opt=optim(par=c(20), fn=function(s) -log(sample_states(xtt, Ptt, xttm1, Pttm1, s, 99)[[4]]*obs_likelihood(phi, sigv, sigw, y)[99]) -
            transition_func(x100, s) - obs_func(x100, y[100]), method='Brent', lower=-10, upper=10)

exp(-opt$value)

sample_states(xtt, Ptt, xttm1, Pttm1, x100, 100)[[4]]*obs_likelihood(phi, sigv, sigw, y)[100]
