source("pp_sim.r")
test = pp_sim.f(S = 50000, iter = 10, Rp_stock = 20000,alpha_stock = 20, mort_r = 0.25, mort_r_harv = 0.5)


S = 1000
N = 500
iter = 150
num.loci =5
plotorno = 1
num.alleles = 2
sd.alleles = 0.1
recomb=1
mort_r = 0.2
mort_r_harv = 0.22
mat_prop = c(0.3,0.4,0.5,1)
age_repr = c(2,3,4,5)
alpha_stock = 50
Rp_stock = 5000