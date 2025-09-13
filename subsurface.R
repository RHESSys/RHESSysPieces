
# by Louis Graup
# 3/5/2024
# calculate soil drainage characteristics from parameter values
# using internal RHESSys functions

library(tidyverse)

K_SAT_Z = function(Ksat0, m, z) {
  # basic exponential profile with depth
  K_SAT_Z = Ksat0 * exp(-z/m)
}

TR = function(Ksat0, m) {
  # transmissivity defined in Tague & Band (2004)
  # as integration of K_sat_z over depth
  # see mathematical derivation in RHESSys pieces
  TR = Ksat0 * m
}

# basic assumptions
z = seq(0, 1, .01)

# saturated hydraulic conductivity (m/s)
Ksat_0_v = 3 # RHESSys default is 3
# and its decay with depth (m)
m_v = 0.12 # RHESSys default is 0.12

plot(K_SAT_Z(Ksat_0_v, m_v, z), -z)


Ksat0 = 10^seq(-1, 2, .01)
m = seq(0.01, 1, .01)

p = merge(Ksat0, m)
colnames(p) = c("K","m")

# transmissivity
p$tr = TR(p$K, p$m)
ggplot(p, aes(K, m, color=log10(tr)))+geom_point()+scale_x_log10()


# toy models

K = 10
m = .25
Tr = K * m

# 1-d model

S = c(10, rep(0, 4))

for (i in 2:5) {
  
  S[i] = max(S[i-1] - Tr, 0)
  
}


# 2-d model

U = c(10, rep(0, 9))
R = U
Uout = 0
Qout = 0

for (i in 2:100) {
  
  Rin = min(Tr, U[i-1])
  U[i] = U[i-1] - Rin
  Uout = Uout + Rin
  
  Rout = min(Tr, R[i-1]+Rin)
  R[i] = R[i-1] + Rin - Rout
  Qout = Qout + Rout
}

Q_diff = which(R==0)[1] - which(U==0)[1]


TOY_MODEL = function(K, m) {
  Tr = K * m
  
  U = c(10, rep(0, 9))
  R = U
  Uout = 0
  Qout = 0
  
  for (i in 2:100) {
    
    Rin = min(Tr, U[i-1])
    U[i] = U[i-1] - Rin
    Uout = Uout + Rin
    
    Rout = min(Tr, R[i-1]+Rin)
    R[i] = R[i-1] + Rin - Rout
    Qout = Qout + Rout
  }
  
  Q_diff = which(R==0)[1] - which(U==0)[1]
  
  return(list(Q_sum=Qout, Q_diff=Q_diff))
  
}

res = pmap_df(list(p$K, p$m), TOY_MODEL)
df = cbind(p, res)

ggplot(df, aes(K, Q_diff))+geom_point()

ggplot(df, aes(K, Q_sum))+geom_point()+scale_x_log10()
