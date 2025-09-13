## by Louis Graup
## 10/10/2021
## conceptual model of subsurface flow in RHESSys

library(tidyverse)

## functions
KSAT_Z = function(z, Ksat0, m) {
  # saturated hydraulic conductivity (m/day) and its decay with depth (m)
  # from eq. 8 in Tague & Band (2004)
  # z = array of depth values (in m)
  KSAT_Z = Ksat0 * exp(-z/m)
}

TRANS = function(m, Ksat0, soil_depth, sat_z) {
  # analytical solution for transmissivity (m^2/day)
  # defined as saturated hydraulic conductivity (as above)
  # integrated over depth of saturated zone
  # derived from eq. 29 in Tague & Band (2004)
  # sat_z = depth to saturated zone (m)
  TRANS = Ksat0 * m * (exp(-sat_z/m) - exp(-soil_depth/m))
}


## basic assumptions
soil_depth = 10 # [m] RHESSys default is 200
porosity = .5 # [-] RHESSys default is .435

# ranges taken from Dingman 3rd ed.
m = 10^seq(-2, 1, by=0.1) # [m] RHESSys default is .12
Ksat0 = 10^seq(-1, 3, by=0.1) # [m/day] RHESSys default is 3

p = merge(m, Ksat0, all=T)
colnames(p) = c("m","K")

# transmissivity
# assumes water table is at soil surface (sat_z = 0)
p$tr = mapply(TRANS, p$m, p$K, MoreArgs=list(soil_depth=soil_depth, sat_z=0))
ggplot(p, aes(log10(m), log10(K), fill=log10(tr)))+geom_tile()


# example for sandy loam soil
z = seq(0, soil_depth, .001)
m = .5
Ksat0 = 100
Ksat_z = KSAT_Z(z, Ksat0, m)
plot(Ksat_z, -z)
plot(z, log10(Ksat_z))
plot(TRANS(m, Ksat0, soil_depth, z), -z)

## TOY MODELS
# soil assumptions defined from above example

# simple 1-d version

# assumes no input
Qin = rep(0, 100)
Qout = Qin*0
sat_def = Qin*0 # saturation deficit [m]
sat_def[1] = soil_depth*porosity*.5

for (i in 2:100) {
  
  # calculate transmissivity, limited by available saturated water
  Qout[i] = max(min(TRANS(m, Ksat0, soil_depth, sat_def[i-1]/porosity), soil_depth*porosity-sat_def[i-1]), 0)
  
  # update saturation deficit (inverse balance since we're tracking a deficit)
  sat_def[i] = sat_def[i-1] - Qin[i] + Qout[i]
}
# plots
plot(sat_def)
plot(Qout[2:100])


# more complicated 2-d version

# assumes upslope (U) patch only receives rainfall, 
# while riparian (R) patch receives both rainfall and lateral input from upslope
# initialize saturation deficit (sd [m]) to be 50% of saturation

df = read_csv("data/precip.csv") # rainfall data for one water year
df2 = data.frame(ind=1:length(df$Precip_m), Qin_U=df$Precip_m, Qout_U=0, Qin_R=0, Qout_R=0, sd_U=soil_depth*porosity*.5, sd_R=soil_depth*porosity*.5)

for (i in 2:length(df2$ind)) {
  
  # fluxes
  df2$Qout_U[i] = max(min(TRANS(m, Ksat0, soil_depth, df2$sd_U[i-1]/porosity), soil_depth*porosity-df2$sd_U[i-1]), 0) # upslope transmissivity
  df2$Qout_R[i] = max(min(TRANS(m, Ksat0, soil_depth, df2$sd_R[i-1]/porosity), soil_depth*porosity-df2$sd_R[i-1]), 0) # riparian transmissivity
  df2$Qin_R[i] = df2$Qin_U[i] + df2$Qout_U[i] # local and lateral inputs
  
  # update stores
  df2$sd_U[i] = max(df2$sd_U[i-1] - df2$Qin_U[i] + df2$Qout_U[i], 0) # not allowing negative deficit to supply return flow
  df2$sd_R[i] = max(df2$sd_R[i-1] - df2$Qin_R[i] + df2$Qout_R[i], 0) # normally as saturation excess overland flow
}
# plots
ggplot(df2, aes(x=ind, y=sd_U, color="U"))+geom_line()+geom_line(aes(x=ind, y=sd_R, color="R"))+labs(x="Day of WY",y="Sat Def (m)")
ggplot(df2, aes(x=ind, y=Qout_U, color="U"))+geom_line()+geom_line(aes(x=ind, y=Qout_R, color="R"))+labs(x="Day of WY", y="Q out (m/day)")


# most complex 3-d version

# same assumptions as 2-d but includes groundwater model for vertical drainage
gw1 = .1 # proportion of infiltrated water that bypasses saturated zone and recharges groundwater store as macropore flow
gw2 = .1 # proportion of groundwater store that is routed directly to riparian area

df3 = mutate(df, ind=1:length(df$Precip_m), Qin_U=df$Precip_m, Qout_U=0, Qdrain_U=0, Qin_R=0, Qout_R=0, Qdrain_R=0, 
             gw_store=0, gw_out=0, sd_U=soil_depth*porosity*.5, sd_R=soil_depth*porosity*.5, streamflow=0)

for (i in 2:length(df3$Precip_m)) {
  
  # fluxes
  df3$Qdrain_U[i] = df3$Qin_U[i] * gw1 # calculate amount of bypass flow for upslope
  df3$Qdrain_R[i] = df3$Qdrain_U[i] # bypass flow for riparian = upslope since same P input
  df3$Qout_U[i] = TRANS(m, Ksat0, soil_depth, df3$sd_U[i-1]/porosity) # upslope transmissivity
  df3$Qout_R[i] = TRANS(m, Ksat0, soil_depth, df3$sd_R[i-1]/porosity) # riparian transmissivity
  df3$gw_out[i] = df3$gw_store[i-1] * gw2 # calculate amount of deep groundwater flow
  df3$Qin_R[i] = df3$Precip_m[i] + df3$Qout_U[i] + df3$gw_out[i] # local, lateral and gw inputs
  df3$streamflow[i] = df3$Qout_R[i] # calculate streamflow
  
  # update stores
  df3$gw_store[i] = df3$gw_store[i-1] + df3$Qdrain_U[i] + df3$Qdrain_R[i] - df3$gw_out[i]
  df3$sd_U[i] = df3$sd_U[i-1] - df3$Qin_U[i] + df3$Qout_U[i] + df3$Qdrain_U[i]
  df3$sd_R[i] = df3$sd_R[i-1] - df3$Qin_R[i] + df3$Qout_R[i] + df3$Qdrain_R[i]
}
# plots
ggplot(df3, aes(x=ind, y=sd_U, color="U"))+geom_line()+geom_line(aes(x=ind, y=sd_R, color="R"))+labs(x="Day of WY", y="Saturation Deficit (m)")
ggplot(df3[2:365,], aes(x=ind, y=Qout_U, color="U"))+geom_line()+geom_line(aes(x=ind, y=Qout_R, color="R"))+scale_y_log10()+labs(x="Day of WY", y="Q out (m/day)")
ggplot(df3[2:365,], aes(x=ind, y=Qdrain_U, color="U"))+geom_line()+geom_line(aes(x=ind, y=Qdrain_R, color="R"))+labs(x="Day of WY", y="Recharge (m/day)")
ggplot(df3, aes(x=ind, y=gw_store))+geom_line()+labs(x="Day of WY", y="Groundwater Store (m)")
ggplot(df3, aes(x=ind, y=streamflow))+geom_line()+labs(x="Day of WY", y="Streamflow (mm/day)")

# water balance
df3$sdU_diff = c(0, -diff(df3$sd_U)) # flux into or out of upslope saturated zone
df3$sdR_diff = c(0, -diff(df3$sd_R)) # flux into or out of riparian saturated zone
df3$gw_diff = c(0, diff(df3$gw_store)) # flux into or out of groundwater store
df3$watbal = with(df3, 2*Precip_m - streamflow - gw_diff - sdU_diff - sdR_diff) # daily
# check for closure
all(abs(df3$watbal)<1e-6)

# visualized
df3_cs = cumsum(df3)
df3_cs = df3_cs %>% mutate(Precip=2*Precip_m, ind=df3$ind, sdU=sdU_diff, sdR=sdR_diff, GW=gw_diff)
df3_long = gather(select(df3_cs, ind, streamflow, GW, sdU, sdR), key="Flux", value="Val", -ind)
df3_long$Flux = factor(df3_long$Flux, levels=c("sdU", "sdR", "streamflow", "GW"))

ggplot()+geom_area(data=df3_long, aes(x=ind, y=Val, fill=Flux))+
  geom_line(data=df3_cs, aes(x=ind, y=Precip, color="Precip"))+scale_color_manual(name="", values=c("Precip"="black"))+
  labs(x="Day of WY", y="Cumulative Flux (m)", title="Annual Water Budget by Day of WY")+theme(plot.title=element_text(hjust=.5))


## sensitivity analysis on 2-d model

# define function to run model with input parameters
SUBSURF_2D = function (m, Ksat0, df, soil_depth, porosity) {
  # same as above
  df2 = mutate(df, ind=1:length(df$water_in), Qin_U=df$water_in, Qout_U=0, Qin_R=0, Qout_R=0, sd_U=soil_depth*porosity*.5, sd_R=soil_depth*porosity*.5)
  for (i in 2:length(df2$ind)) {
    df2$Qout_U[i] = max(min(TRANS(m, Ksat0, soil_depth, df2$sd_U[i-1]/porosity), soil_depth*porosity-df2$sd_U[i-1]), 0) # upslope transmissivity
    df2$Qout_R[i] = max(min(TRANS(m, Ksat0, soil_depth, df2$sd_R[i-1]/porosity), soil_depth*porosity-df2$sd_R[i-1]), 0) # riparian transmissivity
    df2$Qin_R[i] = df2$Qin_U[i] + df2$Qout_U[i] # local and lateral inputs
  
    df2$sd_U[i] = max(df2$sd_U[i-1] - df2$Qin_U[i] + df2$Qout_U[i], 0)
    df2$sd_R[i] = max(df2$sd_R[i-1] - df2$Qin_R[i] + df2$Qout_R[i], 0)
  }
  # calculate metrics for model intercomparison
  Qout_U = sum(df2$Qout_U)
  Qout_R = sum(df2$Qout_R)
  Qout_diff = 100*(df2$Qout_R-df2$Qout_U)/df2$Qout_U
  Qday = ifelse(sum(Qout_diff, na.rm=T)==0, 0, which(diff(Qout_diff[2:length(Qout_diff)])<1e-1)[1] + 2)
  CM_U = sum(df2$Qout_U*df2$ind)/sum(df2$Qout_U)
  CM_R = sum(df2$Qout_R*df2$ind)/sum(df2$Qout_R)
  Qcum_U = which(cumsum(df2$Qout_U[which(df$water_in>0)[1]:length(df2$Qout_U)])>sum(df2$Qin_U))[1]
  Qcum_R = which(cumsum(df2$Qout_R[which(df$water_in>0)[1]:length(df2$Qout_R)])>sum(df2$Qin_R))[1]
  return(c(Qout_U, Qout_R, Qday, CM_U, CM_R, CM_diff=CM_R-CM_U, Qcum_U, Qcum_R))
}

## run the model for all parameter combinations 
## either using the provided input data of a single pulse or artificial input
df4 = read_csv("data/pulse.csv") # pulse generated from snowmelt and rain
# lag = 30
# df4 = data.frame(water_in=c(rep(0,lag),.01, rep(0,lag), .01, rep(0,998-lag*2)))

Qout = data.frame(t(mapply(SUBSURF_2D, p$m, p$K, MoreArgs=list(df=df4, soil_depth=soil_depth, porosity=porosity))))
colnames(Qout) = c("Qout_U", "Qout_R", "Qday", "CM_U", "CM_R", "CM_diff", "Qcum_U", "Qcum_R")
sens_met = cbind(p, Qout)
sens_met$Qout_diff = with(sens_met, 100*(Qout_R-Qout_U)/Qout_U)
sens_met$Qcum_diff = with(sens_met, Qcum_R - Qcum_U)

# visualize parameter relationships
ggplot(sens_met, aes(log10(m), Qout_U, color="U"))+geom_point()+geom_point(aes(log10(m), Qout_R, color="R"))
ggplot(sens_met, aes(log10(K), Qout_U, color="U"))+geom_point()+geom_point(aes(log10(K), Qout_R, color="R"))
ggplot(sens_met, aes(log10(m), log10(K), fill=Qout_R))+geom_tile()

ggplot(sens_met, aes(log10(tr), Qout_U, color="U"))+geom_point()+geom_point(aes(log10(tr), Qout_R, color="R"))

# annual percent difference between upslope and riparian outflow
ggplot(sens_met, aes(log10(tr), Qout_diff, color=log10(m)))+geom_point()
ggplot(sens_met, aes(log10(m), Qout_diff, color=log10(K)))+geom_point()
ggplot(sens_met, aes(log10(K), Qout_diff, color=log10(m)))+geom_point()

# days to move upslope outflow into riparian area
ggplot(sens_met, aes(log10(tr), Qday, color=log10(m)))+geom_point()
ggplot(sens_met, aes(log10(m), Qday, color=log10(K)))+geom_point()
ggplot(sens_met, aes(log10(K), Qday, color=log10(m)))+geom_point()

# center of mass for upslope and riparian outflow
ggplot(sens_met, aes(log10(tr), CM_U, color="U"))+geom_point()+geom_point(aes(log10(tr), CM_R, color="R"))
ggplot(sens_met, aes(log10(m), CM_U, color="U"))+geom_point()+geom_point(aes(log10(m), CM_R, color="R"))
ggplot(sens_met, aes(log10(K), CM_U, color="U"))+geom_point()+geom_point(aes(log10(K), CM_R, color="R"))

# difference in center of mass between upslope and riparian outflow
ggplot(sens_met, aes(log10(tr), CM_diff, color=log10(m)))+geom_point()
ggplot(sens_met, aes(log10(m), CM_diff, color=log10(K)))+geom_point()
ggplot(sens_met, aes(log10(K), CM_diff, color=log10(m)))+geom_point()

# days for cumulative outflow to evacuate inflow
ggplot(sens_met, aes(log10(tr), Qcum_U, color="U"))+geom_point()+geom_point(aes(log10(tr), Qcum_R, color="R"))
ggplot(sens_met, aes(log10(m), Qcum_U, color="U"))+geom_point()+geom_point(aes(log10(m), Qcum_R, color="R"))
ggplot(sens_met, aes(log10(K), Qcum_U, color="U"))+geom_point()+geom_point(aes(log10(K), Qcum_R, color="R"))

# difference in days for cumulative outflow
ggplot(sens_met, aes(log10(tr), Qcum_diff, color=log10(K)))+geom_point()
ggplot(sens_met, aes(log10(m), Qcum_diff, color=log10(K)))+geom_point()
ggplot(sens_met, aes(log10(K), Qcum_diff, color=log10(m)))+geom_point()
