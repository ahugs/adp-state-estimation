devtools::install_github("nickpoison/astsa")
require(astsa)
require(glmnet)

adp <- function(y, transition_func, obs_func, freqs=NULL, rand_n=NULL,
                covar_funcs=NULL, iter=100, prep_iter=10){

  N = length(y)

  # estimate AR(1) parameters and Kalman filter for different frequencies of data
  ar1_params = list()
  xtt = list()
  Ptt = list()
  xttm1 = list()
  Pttm1 = list()

  if(!is.null(freqs)){
    for(i in 1:length(freqs)){
      y_freq = ts(y[seq(1, N, freqs[i])])
      ar1_params[[i]] = estimate_params(y_freq)
      kf = Kfilter0(length(y_freq), y_freq, 1, mu0=0,
                    Sigma0=(ar1_params[[i]]['sigw', 'estimate']^2/(1-ar1_params[[i]]['phi', 'estimate']^2)),
                    ar1_params[[i]]['phi', 'estimate'], ar1_params[[i]]['sigw', 'estimate'],
                    ar1_params[[i]]['sigv', 'estimate'])
      xtt[[i]] = rep(c(kf$xf), 1, each=freqs[i])[1:N]
      Ptt[[i]] = rep(c(kf$Pf), 1, each=freqs[i])[1:N]
      xttm1[[i]] = rep(c(kf$xp), 1, each=freqs[i])[1:N]
      Pttm1[[i]] = rep(c(kf$Pp), 1, each=freqs[i])[1:N]
    }
  }
  else if(!is.null(rand_n)) {
      for(i in 1:rand_n){
        temp = t(data.frame(sigw=rgamma(1,1,1), sigv=rgamma(1,1,1), phi=rnorm(1,0,0.4)))
        colnames(temp) = 'estimate'
        ar1_params[[i]] = temp
        kf = Kfilter0(length(y), y, 1, mu0=0,
                      Sigma0=(ar1_params[[i]]['sigw', 'estimate']^2/(1-ar1_params[[i]]['phi', 'estimate']^2)),
                      ar1_params[[i]]['phi', 'estimate'], ar1_params[[i]]['sigw', 'estimate'],
                      ar1_params[[i]]['sigv', 'estimate'])
        xtt[[i]] = kf$xf
        Ptt[[i]] = kf$Pf
        xttm1[[i]] = kf$xp
        Pttm1[[i]] = kf$Pp
      }
    }

  # initialize approximation with no data
  v_N = length(ar1_params)
  if(!is.null(covar_funcs)){
    v_init = rep(list(v=numeric(0)), length(covar_funcs))
    v_names = names(covar_funcs)
  } else{
    v_init = rep(list(v=numeric(0)), v_N)
    v_names = paste(rep('v', v_N), seq(1, v_N), sep='')
  }
  names(v_init) <- v_names
  df = list(data.frame(v=numeric(0), v_init))
  V = rep(df, length(y))

  S = data.frame(matrix(0, length(y), iter))
  models=replicate(iter, list())
  obs_llike = obs_likelihood(phi, sigv, sigw, y, log=T)
  for(i in 1:iter){

    # Select random final state
    S[length(y), i] = rnorm(1, 0, 1)

    for(j in length(y):2){
      if(nrow(V[[j-1]]) < prep_iter) {
        models[[i]][[j-1]] = NA
        max_func <- function(s){
          -transition_func(S[j, i], s) - obs_func(S[j, i], y[j]) -
            sample_states(xtt[[1]], Ptt[[1]], xttm1[[1]], Pttm1[[1]], s, j-1, log=T)[[4]] -
            obs_llike[j-1]

        }
      } else {
        models[[i]][[j-1]] = cv.glmnet(x=data.matrix(V[[j-1]][, v_names]),
                                       y=V[[j-1]]$v,
                                       family='gaussian',
                                       nfolds=5)
        max_func <- function(s){

          covars = sapply(seq(1, v_N), function(v) {sample_states(xtt[[v]], Ptt[[v]], xttm1[[v]], Pttm1[[v]], s, j-1, log=T)[[4]] +
              obs_llike[j-1]})
          if(!is.null(covar_funcs)){
            covars = sapply(covar_funcs, function(func){func(covars)})
            names(covars) = names(covar_funcs)
          } else {
            names(covars) = v_names
          }
          predict_df = data.frame(as.list(covars))

          -transition_func(S[j, i], s) - obs_func(S[j, i], y[j])  - predict(models[[i]][[j-1]], newx=data.matrix(predict_df),
                                                                            type='response')

        }
      }
      opt = optim(S[j, i], max_func, method='Brent', lower=-20, upper=20)
      S[j - 1, i] = opt$par
      covars = sapply(seq(1, v_N), function(v) {sample_states(xtt[[v]], Ptt[[v]], xttm1[[v]],
                                                              Pttm1[[v]], S[j, i], j, log=T)[[4]] + obs_llike[j]})
      if(!is.null(covar_funcs)){
        covars = sapply(covar_funcs, function(func){func(covars)})
        names(covars) = names(covar_funcs)
      } else {
        names(covars) = v_names
      }
      V[[j]] = rbind(V[[j]], data.frame(v=c(-opt$value), as.list(covars)))

    }
  }
  return(list(V, xtt, Ptt, xttm1, Pttm1, S, models, ar1_params))
}

# Initial Estimates
estimate_params <- function(y) {
  u = ts.intersect(y, lag(y,-1), lag(y,-2))
  varu = var(u)
  coru = cor(u)
  phi = min(coru[1,3]/coru[1,2],0.99)
  q = (1-phi^2)*varu[1,2]/phi
  r = max(varu[1,1] - q/(1-phi^2), 0.1)
  # print(paste('varu ', varu[1,1]))
  # print(paste('rrac ', q/(1-phi^2)))
  (init.par = c(phi, sqrt(q), sqrt(r)))

  # Function to evaluate the likelihood
  Linn=function(para){
    phi = para[1]; sigw = para[2]; sigv = para[3]
    Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
    kf = Kfilter0(length(y),y,1,mu0=0,Sigma0,phi,sigw,sigv)
    return(kf$like)
  }

  # Estimation
  # print(init.par)
  (est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE))
  SE = sqrt(diag(solve(est$hessian)))
  return(data.frame(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)))
}

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

obs_likelihood <- function(phi, sigv, sigw, y, log=F) {
  mu = rep(NA, length(y))
  sigma2 = rep(NA, length(y))
  llike = rep(NA, length(y))

  mu[1] = phi * x[1]
  sigma2[1] = sigw^2
  llike[1] = dnorm(y[1], mu[1], sqrt(sigma2[1] + sigv^2), log=log)
  for(i in 2:length(y)){
    mu[i] = (phi * y[i - 1] * sigma2[i - 1] + mu[i - 1] * phi * sigv^2)/(sigma2[i - 1] + sigv^2)
    sigma2[i] = (sigw^2 * sigv^2 + sigw^2 * sigma2[i - 1] + sigv^2 * sigma2[i - 1] * phi^2)/(sigma2[i - 1] + sigv^2)
    if(log){
      llike[i] = llike[i-1] + dnorm(y[i], mu[i], sqrt(sigma2[i] + sigv^2), log=T)
    }else{
      llike[i] = llike[i-1] * dnorm(y[i], mu[i], sqrt(sigma2[i] + sigv^2))

    }

  }
  return(llike)
}



