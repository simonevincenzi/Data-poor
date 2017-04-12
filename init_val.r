# source("pp_sim.r")
test = pp_sim.f(S = 50000, iter = 100, Rp_stock = 20000,alpha_stock = 20, mort_r = 0.2, mort_r_harv = 0.6)

# Default parameter values
S = 1000 # carrying capacity and initial number of salmon with ages from 0 to 5
N = 10000 # Currently not used
iter = 150 # number of years of simulation
num.loci = 5 # loci
num.alleles = 2 # alleles for loci, 2 means SNP-type
sd.alleles = 0.1 # standard deviation of distribution of alleles
recomb = 1 # 1 means full recombination for offspring
mort_r = 0.2 # average natural mortality rate in the ocean
mort_r_harv = 0.27 # average harvest rate in the ocean, after natural mortality
mat_prop = c(0.05,0.5,0.93,1) # probability of sexual maturity at age 2, 3, 4, 5
age_repr = c(2,3,4,5) # age at which fish can be sexually mature
alpha_stock = 50 # alpha and Rp are parameter values for the Ricker's S-R function. See stock.r for plots of S-R function with different values of alpha and Rp
Rp_stock = 5000