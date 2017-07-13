# source("pp_sim.r")
# 
# # test = pp_sim.f(S = 5000, iter = 20, Rp_stock = 20000,alpha_stock = 20, mort_r = 0.2, mort_r_harv = 0.6)
# 
# test = pp_sim.f(S = 1000, iter = 100, Rp_stock = 1000,alpha_stock = 20, mort_r = 0.2, mort_r_harv = 0.55, num.loci = 20)
# 
# # Default parameter values
# S = 1000 # carrying capacity and initial number of salmon with ages from 0 to 5
# N = 10000 # Currently not used
# iter = 150 # number of years of simulation
# num.loci = 5 # loci
# num.alleles = 2 # alleles for loci, 2 means SNP-type
# sd.alleles = 0.1 # standard deviation of distribution of alleles
# recomb = 1 # 1 means full recombination for offspring
# mort_r = 0.2 # average natural mortality rate in the ocean
# mort_r_harv = 0.27 # average harvest rate in the ocean, after natural mortality
# mat_prop = c(0.05,0.5,0.93,1) # probability of sexual maturity at age 2, 3, 4, 5
# age_repr = c(2,3,4,5) # age at which fish can be sexually mature
# alpha_stock = 50 # alpha and Rp are parameter values for the Ricker's S-R function. See stock.r for plots of S-R function with different values of alpha and Rp
# Rp_stock = 5000 #
# L_inf = 3750
# k_vb = 0.057
# t0_vb = -0.21
# alpha_w = 1.44
# beta_w = 3.12
# sd_e = 0.2 # standard deviation of the escapement observation model
# sex_repro = 0





source("demo_sim.r")

# test = pp_sim.f(S = 5000, iter = 20, Rp_stock = 20000,alpha_stock = 20, mort_r = 0.2, mort_r_harv = 0.6)

test = demo_sim.f(S = 10000, N = 1000,iter = 100,mort_r = 0.2,
                  mort_r_harv = 0.55, strategy = 3, mult_targ = 5
                  , alpha_stock = 10,Rp_stock = 500, 
                  L_inf = 3750, k_vb = 0.057, t0_vb = -0.21, 
                  alpha_w = 1.44, beta_w = 3.12, sd_= 0.22)
  
  
S = 1000 # carrying capacity and initial number of salmon with ages from 0 to 5
N = 500 # Currently not used
iter = 50 # number of years of simulation
mort_r = 0.2 # average natural mortality rate in the ocean
mort_r_harv = 0.27 # average harvest rate in the ocean, after natural mortality
mat_prop = c(0.05,0.5,0.93,1) # probability of sexual maturity at age 2, 3, 4, 5
age_repr = c(2,3,4,5) # age at which fish can be sexually mature
alpha_stock = 50 # alpha and Rp are parameter values for the Ricker's S-R function. See stock.r for plots of S-R function with different values of alpha and Rp
Rp_stock = 5000 #
L_inf = 3750
k_vb = 0.057
t0_vb = -0.21
alpha_w = 1.44
beta_w = 3.12
sd_e = 0.2 # standard deviation of the escapement observation model



S.vett = 10000
N.vett = 1000
iter.vett = 50
mort_r.vett = 0.2
mort_r_harv.vett = 0.27
alpha_stock.vett = c(0.5,1,1.5,2,5,10,50) # alpha and Rp are parameter values for the Ricker's S-R function. See stock.r for plots of S-R function with different values of alpha and Rp

Rp_stock.vett = c(100,500,1000,1500,2000,2500,5000) #
L_inf.vett = 3750
k_vb.vett = 0.057
t0_vb.vett = -0.21
alpha_w.vett = 1.44
beta_w.vett = 3.12
sd_e.vett = 0.2 



###### remember that the columns must be named!!!

dataforpar = expand.grid(S = S.vett, N = N.vett ,iter = iter.vett , mort_r = mort_r.vett, mort_r_harv = mort_r_harv.vett, alpha_stock = alpha_stock.vett,Rp_stock = Rp_stock.vett, L_inf = L_inf.vett, k_vb = k_vb.vett, t0_vb = t0_vb.vett,alpha_w = alpha_w.vett,beta_w = beta_w.vett,sd_e= sd_e.vett)

dataforpar = as.matrix(do.call("rbind", rep(list(dataforpar), 10)))


y = list()

for (i in 1:nrow(dataforpar)) {
  y[[i]] = as.list(dataforpar[i,])
}



prova = mclapply(1:length(y),function (x) do.call(demo_sim.f,y[[x]]),
                 mc.cores=4)

 
stock_sim.df = readRDS("stock_sim.RDS")  
dataforpar = as.tibble(dataforpar)
for (i in 1:length(stock_sim.df)) {
  stock_sim.df[[i]]$alpha_stock = dataforpar$alpha_stock[i]
  stock_sim.df[[i]]$Rp_stock = dataforpar$Rp_stock[i]
  
}


stock_sim_filt.df = list.filter(stock_sim.df, extinct == 0)

stock_par.df = tibble(alpha = rep(0,length(stock_sim_filt.df)), Rp = 0)

for (i in 1:length(stock_sim_filt.df)) {
  stock_par.df$alpha[i] = stock_sim_filt.df[[i]]$alpha_stock 
  stock_par.df$Rp[i] = stock_sim_filt.df[[i]]$Rp_stock
  
}
  
stock_par.df %>%
  ggplot(aes(x = alpha, y = Rp)) +
  geom_point()

table(stock_par.df)
    Rp
alpha 100 500 1000 1500 2000 2500 5000
  2    0   0    1    1    3    1    3
  5   10  10   10   10   10   10   10
  10  10  10   10   10   10   10   10
  50   8  10   10   10   10   10   10