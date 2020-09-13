SEIR <- function(S_0, E_0, I_0, R_0, n, alpha, gamma, beta, sig_g, sig_y){
  pop_N = S_0 + E_0 + I_0 + R_0
  S = rep(NA, n); S[1] = S_0
  E = rep(NA, n); E[1] = E_0
  I = rep(NA, n); I[1] = I_0
  R = rep(NA, n); R[1] = R_0
  g = rep(NA, n)
  y = rep(NA, n)

  for(i in 2:n){
    g[i] = alpha * E[i-1]/I[i-1] - gamma + rnorm(1, 0, sig_g)
    y[i] = g[i] + rnorm(1, 0, sig_y)
    I[i] = (1 + g[i])*I[i-1]
    S[i] = S[i-1] - beta*S[i-1]*I[i-1]/pop_N
    E[i] = E[i-1] + beta*S[i-1]*I[i-1]/pop_N - alpha*E[i-1]
    R[i] = R[i-1] + gamma*I[i-1]
  }
  return(list(S, E, I, R, g, y))
}

SEIR(999, 0, 2, 0, 100, 2, 1.5, 1, 1, 1)
