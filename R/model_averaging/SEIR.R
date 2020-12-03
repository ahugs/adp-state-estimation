SEIR <- function(S_0, E_0, I_0, R_0, N, alpha, gamma, beta, sig_g, sig_y){
  pop_0 = S_0 + E_0 + I_0 + R_0
  S = rep(NA, N); S[1] = S_0
  E = rep(NA, N); E[1] = E_0
  I = rep(NA, N); I[1] = I_0
  R = rep(NA, N); R[1] = R_0
  y = rep(NA, N)
  g = rep(NA, N)

  for(i in 2:N){
    g[i] = alpha * E[i-1]/I[i-1] - gamma + rnorm(1, 0, sig_g)
    y[i] = g[i] + rnorm(1, 0, sig_y)
    I[i] = (1 + g[i])*I[i-1]
    S[i] = S[i-1] - beta*S[i-1]*I[i-1]/pop_0
    E[i] = E[i-1] + beta*S[i-1]*I[i-1]/pop_0 - alpha*E[i-1]
    R[i] = R[i-1] + gamma*I[i-1]
  }
  return(list(S=S, E=E, I=I, R=R, g=g, y=y))
}

SEIRFixedStates <- function(S_0, E_0, I_0, R_0, N, alpha, gamma, beta, g){
  pop_0 = S_0 + E_0 + I_0 + R_0
  S = rep(NA, N); S[1] = S_0
  E = rep(NA, N); E[1] = E_0
  I = rep(NA, N); I[1] = I_0
  R = rep(NA, N); R[1] = R_0

  for(i in 2:N){
    I[i] = (1 + g[i-1])*I[i-1]
    S[i] = S[i-1] - beta*S[i-1]*I[i-1]/pop_0
    E[i] = E[i-1] + beta*S[i-1]*I[i-1]/pop_0 - alpha*E[i-1]
    R[i] = R[i-1] + gamma*I[i-1]
  }
  return(list(S=S, E=E, I=I, R=R, g=g))
}
