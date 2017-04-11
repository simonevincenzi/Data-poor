## Simulation model for a Pacific salmon species

## we need to start from 
# log(Rt)=log(St)+ at –bSt + wt.
# R = recruitment
# S = spawners
# a = density-independent productivity
# b = strength of density dependence
# w = autocorrelated process error (environmental stochasticity)

# let's start with up to 5 age classes

S = 500
N = 500
iter = 50
selection = 0.05
sdopt = 0.005
num.loci =5
plotorno = 1
meanopt = 0
major=0
num.alleles = 2
sd.alleles = 0.1
recomb=1
mutation=0
mu.mut=0.0002
mut.alfa=0.05
catastr = 0
p.pre = 0.05
p.post = 0.07
catastr.mort = 0.3
age.sex.mat = 2
avg.surv = 0.9
sons.mean = 2
first.phase.opt = 150
var.amb.init = 0.1
mort_r = 0.2
mort_r_harv = 0
mat_prop = c(0.5,0.7,0.9,1)
age_repr = c(2,3,4,5)
a_stock = 0.9
Rp_stock = 150  
sex_repro = 1
