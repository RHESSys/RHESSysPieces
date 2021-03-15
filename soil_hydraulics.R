# by Louis Graup
# 3/4/2021
# calculate soil characteristics from parameter values
# using internal RHESSys functions

library(tidyverse)

FIELD_CAP = function(pa, po, porosity) {
  # definition from Dingman 2nd ed., equation 6.19
  FIELD_CAP = porosity * (pa/3.4) ^ po
}

WILTING_PT = function(pa, po, porosity, psi_close) {
  # ecohydrologic definition derived from Campbell model of psi-theta relationship below
  # see mathemetical derivation in RHESSys pieces
  WILTING_PT = porosity * exp(-1*log(-100*psi_close/pa)*po)
}

PSI = function(pa, po, S) {
  # Campbell (1974) model from Dingman 3rd ed., table 7.2
  # S is degree of saturation
  # 1 m water tension = 10000 Pa = 0.01 MPa
  PSI = -.01*pa*S^(-1/po) # in MPa
}

m_LWP = function(psi, psi_open, psi_close) {
  # original Biome BGC linear approximation
  if (psi>psi_open)
    m_LWP = 1
  else if (psi<psi_close)
    m_LWP = 0
  else
    m_LWP = (psi - psi_close) / (psi_open - psi_close)
  return(m_LWP)
}

# m_LWP = function(psi, psi_open, psi_close, slope=.25, int=1, t=psi_open, exp=2) {
#   # RHESSys defaults assume psi_slope=.2, psi_intercept=1.0, psi_threshold=-1 MPa
#   # exp=1 is linear, too shallow, default function values are my guess
#   # only useful if you have measurements with which to fit a curve
#   if (psi>psi_open)
#     m_LWP = 1
#   else if (psi<psi_close)
#     m_LWP = 0
#   else {
#     if (psi<t)
#       m_LWP = (slope*(psi-t)+int)^exp
#     else
#       m_LWP = 1
#   }
# }

# basic assumptions
porosity = 0.45 # RHESSys default is 0.435
psi_open = -0.5 # MPa (RHESSys default is -.65)
psi_close = -2.5 # MPa (RHESSys default is -2.5)

# ranges from Dingman 3rd ed., table 7.4
psi_air_entry = seq(.1, 1, by=.01) # (m) (RHESSys default is 0.218)
pore_size_index = seq(.15, .3, by=.01) # (-) (RHESSys default is 0.204)

p = merge(psi_air_entry, pore_size_index, all=T)
colnames(p) = c("pa","po")

# field capacity
p$fc = mapply(FIELD_CAP, p$pa, p$po, MoreArgs=list(porosity=porosity))
ggplot(p, aes(pa, po, fill=fc))+geom_tile()

# wilting point
p$wp = mapply(WILTING_PT, p$pa, p$po, MoreArgs=list(porosity=porosity, psi_close=psi_close))
ggplot(p, aes(pa, po, fill=wp))+geom_tile()

# available water capacity (could be multiplied by rooting depth to obtain PAWSC)
p$awc = p$fc - p$wp
ggplot(p, aes(pa, po, fill=awc))+geom_tile()

# psi or soil water potential
p$psi = mapply(PSI, p$pa, p$po, MoreArgs=list(S=.2)) # assume storage is .2
ggplot(p, aes(pa, po, fill=log10(-psi)))+geom_tile()

# example for sandy loam (RHESSys defaults)
pa = .218
po = .204
fc = FIELD_CAP(pa, po, porosity)
wp = WILTING_PT(pa, po, porosity, psi_close)
#S = seq(.1, 1, by=.01) # relative saturation
S = seq(wp/.45, fc/.45, by=.001)
psi = PSI(pa, po, S=S) # MPa (convert to cm by *10000)
plot(S, psi)
#plot(S, log10(-psi)) # easier to visualize if S=(.1,1)

# leaf water potential multiplier on stomatal conductance from Jarvis model in RHESSys
m_gs = sapply(psi, FUN=m_LWP, psi_open, psi_close)
plot(psi, m_gs)
plot(S, m_gs)


# use pedotransfer functions to derive van Genutchten parameters from soil characteristics
# equations from Dingman 3rd ed., table 7.3 and box 8.1
# sand fraction = 70%, clay fraction = 5%, assumes organic carbon fraction = 1%, bulk density = 1.4 g/cm^3
SAND = 70
CLAY = 5
C = 1
RHO_D = 1.4

porosity_vg = 1 - (RHO_D/2.65) # from Dingman 3rd ed., eq. 7.5, just for reference

t_residual = 0.015 + 0.005*CLAY + 0.014*C # residual water content
phi_vg = 0.810 - 0.283*RHO_D + 0.001*CLAY # porosity parameter
alpha_pa = exp(-2.486 + 0.025*SAND - 0.352*C - 2.617*RHO_D - 0.023*CLAY) # pressure head parameter (1/cm)
n_vg = exp(0.053 - 0.009*SAND - 0.013*CLAY + 0.00015*SAND^2) # quasi-pore size index
K_sat = 0.000011574 * exp(20.62 - 0.96*log(CLAY) - 0.66*log(SAND) - 0.046*log(C) - 8.43*RHO_D) # saturated hydraulic conductivity (cm/s)
K_sat = K_sat*3600*24 # cm/d

t_vg = n_vg^(-.6*(2+log10(K_sat))) # relative field capacity
t_fc = t_vg*(phi_vg-t_residual)+t_residual # field capacity

t_pwp = t_residual + (phi_vg-t_residual)/((1+(alpha_pa*15000)^n_vg)^((n_vg-1)/n_vg)) # permanent wilting point

t_a = t_fc - t_pwp # available water capacity

# van Genuchten soil water characteristic curve
psi_vg = 10^seq(-2, 25, by=.1) # matric suction (kPa)
theta_vg = (1/(1+(alpha_pa*psi_vg)^n_vg))^(1-(1/n_vg)) # effective saturation
theta = theta_vg*(phi_vg-t_residual)+t_residual # volumetric water content
plot(log10(psi_vg), theta)