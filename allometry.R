
# back-of-the-envelope allometry calculations for ponderosa pine in Sierra Nevada, CA

# allometric equation from McDonald & Skinner (1989)
# to calculate height (ft) from dbh (in)

dbh = seq(1, 36, .5)
height = 10.0333 + 6.3269*dbh - .061*dbh^2
plot(dbh, height)

# calc dbh for tree heights
ht_m = seq(1, 30, .2)
ht_ft = ht_m*3.28083299 # conversion

# solve quadratic equation
a = -.061
b = 6.3269
c = 10.0333 - ht_ft
dbh_in = (-b + sqrt(b^2 - 4*a*c)) / (2*a)

dbh_m = dbh_in / 39.3700787 # conversion
plot(dbh_m*100, ht_m)

# allometry to calculate biomass in g from Law et al (2001)
stem_mass = .293 * ht_m * dbh_m^2 * 407270
stem_Cmass = stem_mass / 2
stem_Cmass_area = stem_Cmass*.02/1000 # kg / m^2

plot(dbh_m*100, stem_Cmass_area)
plot(stem_Cmass_area, ht_m)


# fit power law relationship to derive coefficients for height_to_stem calculation
stem_density = .02 # default
power_lm = lm(log(ht_m) ~ log(stem_Cmass_area/stem_density))
opt_x = coefficients(power_lm)[2]
opt_c = exp(coefficients(power_lm)[1])

Y = opt_c * (stem_Cmass_area / stem_density) ^ opt_x
lines(stem_Cmass_area, Y, col="red")


# determine optimal parameters for RHESSys to match allometric equations
Tree_StemC = read.csv("data/Tree_StemC.csv")$x # from RHESSys

# parameter ranges
coef=seq(.01,2.0,.01) # height_to_stem_coef
x=seq(.1,1.5,.02) # height_to_stem_exp

p = merge(coef, x, all=T)
colnames(p) = c("c","x")

# retrieve yearly heights from stem carbon time-series based on input parameters
RHESSYS_HT = function(c, x, stemc) {
  H = c*(stemc/.02)^x
  H1 = max(H[365:730])
  H5 = max(H[1825:2190])
  H10 = max(H[3650:4015])
  H50 = max(H[18250:18615])
  H100 = max(H[36500:36865])
  H150 = max(H[54422:54787])
  return(list(H1,H5,H10, H50, H100, H150))
}

h = mapply(RHESSYS_HT, p$c, p$x, MoreArgs=list(Tree_StemC))
o = cbind(p, t(h))
colnames(o) = c("c","x","H1","H5","H10","H50","H100","H150")

# filter by "behavioral" growth trajectories from Grulke & Retzlaff (2001)
bh = dplyr::filter(o, H1<1, H5<2, H10<5, H10>1, H50>5, H50<10, H100>10, H100<20, H150<30)

# compare parameter distributions
plot(ecdf(p$c))
lines(ecdf(bh$c), col="red")

plot(ecdf(p$x))
lines(ecdf(bh$x), col="red")

# compare behavioral and fitted parameters

H = .09*(Tree_StemC/stem_density)^1
plot(opt_c*(Tree_StemC/stem_density)^opt_x)
lines(H, col="red")

plot(Tree_StemC, H)
lines(Tree_StemC, opt_c*(Tree_StemC/stem_density)^opt_x, col="red")

plot(stem_Cmass_area, .09*(stem_Cmass_area/stem_density)^1)
lines(stem_Cmass_area, ht_m, col="red")


