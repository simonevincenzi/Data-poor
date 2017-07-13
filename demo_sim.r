
# check if required pacakges are present, otherwise install them

list.of.packages <- c("Rlab", "MASS", "parallel", "tidyverse","MCMCpack","rlist")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rlab)
library(MASS)
library(parallel)
library(tidyverse)
library(MCMCpack)
library(rlist)

################################


demo_sim.f <- function (S = 10000,N = 1000,iter = 150,mort_r = 0.2,mort_r_harv = 0.27 , strategy = 2, mult_targ = 3, alpha_stock = 50,Rp_stock = 5000, L_inf = 3750, k_vb = 0.057, t0_vb = -0.21, alpha_w = 1.44, beta_w = 3.12, sd_e = 0.2)
  
  # S = maximum (and initial number of individuals in the population)
  # N = not used
  # iter = total number of simulation steps/year
  # num.loci = number of genes
  # plotorno = 1 if plot of phenotype/optimum/population
  # num.alleles = num alleles for each locus (randomly drawn at the start of the simulation)
  # sd.alleles = standard deviation of the allelic effects (not used, as it is defined in the function)
  # recomb = 1 for full recombination, 0 for no recombination
  
  
{
  
  extinct = 0
  
  options(error=recover)
  
  options(warn = 0) ### stop when there is warning when warn = 2
  
 
  mat_prop = c(0.05,0.5,0.93,1) # probability of sexual maturity at age 2, 3, 4, 5
  age_repr = c(2,3,4,5) 
  
  mort_r_harv_r = runif(iter,min = (mort_r_harv - mort_r_harv*0.5), max = mort_r_harv + (mort_r_harv*0.5))
  if (mort_r_harv < 0.01) {
    mort_r_harv_r = rep(0,iter)
  }
  mort_r_r = runif(iter,min = (mort_r - mort_r*0.5), max = mort_r + (mort_r*0.5))
  
  if (strategy == 4) {
    
    stock = seq(from = 0, to = 50000, by = 50)
    sons = rep(0, length(stock))
  
    
    for (i in 1:length(stock)) {
      sons[i] = alpha_stock*stock[i]*exp(-a*(stock[i]/(Rp_stock*exp(1))))
    }
    
    stock_df = tibble(stock = stock, sons = sons)
    
    esc_opt = stock_df$stock[which(stock_df$sons == max(stock_df$sons))]
    rm(stock_df)
    rm(stock)
    rm(sons)
  }
  
  
  ############ Define values of some model parameters not included in the function
  
  
  ## the population stays on a matrix (area.pop), each column is an individual. Here below I define
  # the position of entities in the column
  
  pheno.pos = 1  # this means that the first row is the phenotype of the individuals
  age.pos = 2  # second row is the age of the individual
  matur.pos = 3 # third row is maturity status
  linf.pos = 4
  k.pos = 5
  t0.pos = 6
  size.pos = 7 # length of fish
  weight.pos = 8 
  
  length_prop = 8 #  length of entities of individual vector excluding loci (phenotipic value, age, maturity status, and size are the top rows)
  
  
  area.pop <- matrix(0,length_prop,S) ##  matrix space area.pop of S elements, 
  
  
  
  ## define initial ages for the S individuals in the population at time = 1
  
  # initial_ages = round(runif(S,min = 0, max = 5))
  initial_ages = sample(x = 0:5, size = N,replace = T)
  
  area.pop[age.pos,1:N] = initial_ages # assign ages to individuals
  
  area.pop[pheno.pos,1:N] = 1 # 1 is alive
  
  area.pop[linf.pos,1:N] = rnorm(n = N, mean = L_inf, sd = L_inf/3)
  area.pop[k.pos,1:N] = rnorm(n = N, mean = k_vb, sd = k_vb/3)
  area.pop[t0.pos,1:N] = rnorm(n = N, mean = t0_vb, sd = -t0_vb/3)
  
  
  ### assign length and weight at age based on von Bertalanffy's growth parameters ####
  
  area.pop[size.pos,which(area.pop[pheno.pos,]!=0)] = round(sapply(which(area.pop[pheno.pos,]!=0), function(x) {
    area.pop[size.pos,x] = area.pop[linf.pos,x] * (1 - exp(-area.pop[k.pos,x]*(area.pop[age.pos,x]-area.pop[t0.pos,x])))
  }))
  
  area.pop[weight.pos,which(area.pop[pheno.pos,]!=0)] = sapply(which(area.pop[pheno.pos,]!=0), function(x) {
    area.pop[weight.pos,x] = alpha_w * 10^(-8) * area.pop[size.pos,x]^(beta_w)
  })
  
  
  #####create vectors and matrix for later use
 
  mean.age = rep(0,iter) # vector with mean age at each time step
  recruitment = rep(0,iter) # vector with escapement for each year
  escapement_real = rep(0,iter) # vector of escapement (real)
  escapement_obs = rep(0,iter) # vector of escapement (observed)
  
  vettad = rep(0,iter)    #######population size post_mortality at each time step
  
  mean.var.age = matrix(0,iter,2)  ## matrix with iter row and 2 cols (1 = mean age, 2 = sd of age every year)
  
 
  test.pop = list() # list of area matrix for testing
  
  harvest = rep(0,iter) # harvest in numbers of individuals 
  
  harvest_mean_size = rep(0,iter) # mean size of harvested fish 
  harvest_total_weight = rep(0,iter) # total weight of harvested fish 
  harvest_mean_age = rep(0,iter) # mean age of harvested fish
  stock_v = rep(0,iter) # stock
  escapement_mean_size = rep(0,iter)
  escapement_mean_age = rep(0,iter)
  
  ###########
  
 
  
  ##################### START FOR TIME LOOP ########################
  
  
  for(i in 1:iter)  { ##time-steps or cohorts or year
    
    print(i)  # prints time-steps
    
    
    
    #####################################
    
    #########check if the population went extinct####
    
    
    if(length(area.pop[pheno.pos,area.pop[pheno.pos,]!=0])<2){  ##is population extinct?  
      
      print(paste("Extinct at year", i)) #prints the year of extinction
      
      extinct = 1 # when an individual is dead the first column is 0
      
      break}   #close for to check for extinction
    
    ######################################  
    
    
    ##### compute mean age of ocean salmon, only 2 years old and older ######
    
    mean.age.who = which(area.pop[age.pos,]>=2 & area.pop[pheno.pos,]!=0)
    
    mean.age[i] = mean(area.pop[age.pos,mean.age.who],na.rm = T)
    
    
    ############ increase age by 1 if t>1 ##########
    
    if (i>1) {
      
      area.pop[age.pos,area.pop[pheno.pos,]!=0] = area.pop[age.pos,area.pop[pheno.pos,]!=0] + 1 
      recruitment[i] = min(recruit,length(empty))
      
      ########### assign length and weight ###########
      
      area.pop[size.pos,which(area.pop[pheno.pos,]!=0)] = round(sapply(which(area.pop[pheno.pos,]!=0), function(x) {
        area.pop[size.pos,x] = area.pop[linf.pos,x] * (1 - exp(-area.pop[k.pos,x]*(area.pop[age.pos,x]-area.pop[t0.pos,x])))
      }))
      
      area.pop[weight.pos,which(area.pop[pheno.pos,]!=0)] = sapply(which(area.pop[pheno.pos,]!=0), function(x) {
        area.pop[weight.pos,x] = alpha_w * 10^(-8) * area.pop[size.pos,x]^(beta_w)
      })
    }
    ###########
    
    
    ####vector of optimum and mean and sd of phenotype in the population ######
  
    
    
    
    #### close vector of optimum and mean and sd of phenotype in the population ######
    
  
    ############### MEAN AND SD OF AGE #####
    
    who = which(area.pop[pheno.pos,]!=0)  # where the individuals are (because empty columns have zero in first position)
   
    mean.var.age[i,1] = mean(area.pop[age.pos,who]) #mean age
    mean.var.age[i,2] = sd(area.pop[age.pos,who]) #sd of age
    
    
    
    ############mortality##############
    
    
    who_tag = as.matrix(which(area.pop[pheno.pos,]!=0 & area.pop[age.pos,]>=2)) 

    
    ## mortality is sequential, natural and then harvest
    
    ### random number for each adult 
    
    if (length(who_tag) > 0) {
      who_tag = as.tibble(cbind(who_tag,runif(nrow(who_tag))))
    r = 1 - (exp(-mort_r_r[i]))
    who_tag <- who_tag[which(who_tag[,2] < r),]
    if (length(who_tag) > 0) {
      
      area.pop[,who_tag$V1] <- 0
    
    }
    }
    
    ### harvest after natural mortality
    who_tag = 0
    
    who_tag = as.matrix(which(area.pop[pheno.pos,]!=0 & area.pop[age.pos,]>=2)) 

    
    if (strategy == 1) {
    
    if (length(who_tag) > 0) {
      who_tag = as.tibble(cbind(who_tag,runif(nrow(who_tag))))
      r_h = 1 - (exp(-mort_r_harv_r[i]))
      who_tag <- who_tag[which(who_tag[,2] < r_h),]
      if (length(who_tag) > 0) {
        harvest[i] = nrow(who_tag)
        harvest_mean_size[i] = mean(area.pop[size.pos,who_tag$V1], na.rm = T)
        harvest_total_weight[i] = sum(area.pop[weight.pos,who_tag$V1], na.rm = T)
        harvest_mean_age[i] = mean(area.pop[age.pos,who_tag$V1], na.rm = T)
        area.pop[,who_tag$V1] <- 0
        
      }
    } 
      
    } else if (strategy == 2) {
      
      who_tag = as.tibble(cbind(who_tag,runif(nrow(who_tag))))
      target_n = mult_targ * round(mean(escapement_obs[(i-1):max(1,(i-4))],na.rm = T))
      if (i >5) {
      who_tag <- who_tag[sample(1:min(target_n,nrow(who_tag)),size = min(target_n,nrow(who_tag)), replace = F) ,] } else {
        who_tag = who_tag[sample((1:nrow(who_tag)), size = round(0.1*nrow(who_tag)), replace = F),]
        }
      
      if (length(who_tag) > 0) {
        harvest[i] = nrow(who_tag)
        harvest_mean_size[i] = mean(area.pop[size.pos,who_tag$V1], na.rm = T)
        harvest_total_weight[i] = sum(area.pop[weight.pos,who_tag$V1], na.rm = T)
        harvest_mean_age[i] = mean(area.pop[age.pos,who_tag$V1], na.rm = T)
        area.pop[,who_tag$V1] <- 0
        
      }
    } else if (strategy == 3) {
      
      who_tag = as.tibble(cbind(who_tag,runif(nrow(who_tag))))
      target_n = mult_targ * round(min(escapement_obs[(i-1):max(1,(i-4))],na.rm = T))
      if (i >5) {
        who_tag <- who_tag[sample(1:min(target_n,nrow(who_tag)),size = min(target_n,nrow(who_tag)), replace = F) ,] } else {
          who_tag = who_tag[sample((1:nrow(who_tag)), size = round(0.1*nrow(who_tag)), replace = F),]
        }
      
      if (length(who_tag) > 0) {
        harvest[i] = nrow(who_tag)
        harvest_mean_size[i] = mean(area.pop[size.pos,who_tag$V1], na.rm = T)
        harvest_total_weight[i] = sum(area.pop[weight.pos,who_tag$V1], na.rm = T)
        harvest_mean_age[i] = mean(area.pop[age.pos,who_tag$V1], na.rm = T)
        area.pop[,who_tag$V1] <- 0
        
      }
    } else if (strategy == 4) {
      
      who_tag = as.tibble(cbind(who_tag,runif(nrow(who_tag))))
      target_n = max(1,(nrow(who_tag) - 1.8*esc_opt))
      if (i >5) {
        who_tag <- who_tag[sample(1:min(target_n,nrow(who_tag)),size = min(target_n,nrow(who_tag)), replace = F) ,] } else {
          who_tag = who_tag[sample((1:nrow(who_tag)), size = round(0.1*nrow(who_tag)), replace = F),]
        }
      
      if (length(who_tag) > 0) {
        harvest[i] = nrow(who_tag)
        harvest_mean_size[i] = mean(area.pop[size.pos,who_tag$V1], na.rm = T)
        harvest_total_weight[i] = sum(area.pop[weight.pos,who_tag$V1], na.rm = T)
        harvest_mean_age[i] = mean(area.pop[age.pos,who_tag$V1], na.rm = T)
        area.pop[,who_tag$V1] <- 0
        
      }
    } 
    
    
    vettad[i] = length(area.pop[pheno.pos,(area.pop[pheno.pos,]!=0 & area.pop[age.pos,]>=2)]) # populations size in the ocean post-natural mortality and post-harvest
    
   
    ########## close mortality ##############
    
    
    ####################   REPRODUCTION ##########################
    
    ### 
    
    area.pop[matur.pos,] = 0  # set 0 maturity for every fish
    
    for (r_rep in 1:length(mat_prop)) {
      pot_rep = 0
      pp_mat = 0
      pot_rep = which(area.pop[age.pos,]==age_repr[r_rep]) # where are fish of age age_repr[r_rep]
      repr = runif(length(pot_rep))
      pp_mat <- which(repr < mat_prop[r_rep]) # which are the mature fish
      area.pop[matur.pos,pot_rep[pp_mat]] <- 1 ## mature fish get one
    }
    
    
    
    reproductor = which(area.pop[matur.pos,] == 1)       ##### vector of position in the matrix of reproductors
    stock = length(reproductor)  # stock
    stock_v[i] = stock
    
    recruit = round(alpha_stock*stock*exp(-alpha_stock*(stock/(Rp_stock*exp(1))))) # recruitment
    
    # print(recruit)
    
    escapement_real[i] = stock  # stock = escapement = mature fish (real) 
    escapement_mean_size[i] = mean(area.pop[size.pos,reproductor], na.rm = T)
    escapement_mean_age[i] = mean(area.pop[age.pos,reproductor], na.rm = T)
    
    # observed escapement comes from an log-normal observation model log(Et)∼Normal(log(St),rs)
    escapement_obs[i] = round(exp(rnorm(n = 1, mean = log(escapement_real[i] + 1), sd = sd_e)))
    
  #  test.pop[[i]] = area.pop  # save area.pop matrix, mostly for testing
    
    
    if (recruit >1 & stock >=2)   { #open if for reproduction, that is there is at least one recruit
      
      area.pop[,reproductor] = 0  ## all mature fish die (semelparity). reproductors is a vector of columns in the area.pop matrix
      
      empty = which(area.pop[pheno.pos,]==0)  ### identify the available spots for newborns in the matrix. If there are more recruits than empty spots. some recruits die (the number of columns in the matrix is basically the carrying capacity of the system)
      if (length(empty) > recruit) {
        area.pop[pheno.pos, empty[1:recruit]] = 1
        area.pop[linf.pos,empty[1:recruit]] = rnorm(n = length(recruit), mean = L_inf, sd = L_inf/3)
        area.pop[k.pos,empty[1:recruit]] = rnorm(n = length(recruit), mean = k_vb, sd = k_vb/3)
        area.pop[t0.pos,empty[1:recruit]] = rnorm(n = length(recruit), mean = t0_vb, sd = -t0_vb/3)
        } else {
        area.pop[pheno.pos, empty] = 1
        area.pop[linf.pos,empty] = rnorm(n = length(empty), mean = L_inf, sd = L_inf/3)
        area.pop[k.pos,empty] = rnorm(n = length(empty), mean = k_vb, sd = k_vb/3)
        area.pop[t0.pos,empty] = rnorm(n = length(empty), mean = t0_vb, sd = -t0_vb/3)}
    }
      ######### open control for odd number of reproductors
      

      
      year = i  # counter for plot, I cannot use iter because year < iter if the population went extinct   
      
      ####################   CLOSE REPRODUCTION ##########################
      
    }  ##close time
    
  
    ##################### CLOSE FOR TIME LOOP ########################
    
    
    
    ### PLOT
    
    size.title = 15
    line.lwd = 1.2
    size.label.x = 18
    size.text.x = 14
    size.point = 2
    size.label.y = 18
    size.text.y = 14
    size.legend.text = 15
    size.legend.title = 20
    unit.legend.h = 1.8
    unit.legend.w = 1.8
    size.ann = 10
    colour.axis = "gray20"
    colour.theme = "black"
    colour.axis.line = "gray20"
    colour.line = "gray50"
    label.T = "Heterozygosity"
    max_size_dot = 8
    
    ## Theme to be used for all plots
    
    theme.pop =  theme(plot.title = element_text(lineheight=.8, face="bold", size = size.title,hjust = 0.5), 
                       plot.background = element_blank()
                       ,panel.grid.major = element_blank()
                       ,panel.grid.minor = element_blank()
                       ,panel.border = element_blank()
                       ,panel.background = element_blank(),
                       axis.line = element_line(color = 'black'),
                       plot.margin = unit(c(1,2,1,1), "cm"),
                       axis.title.x = element_text(size=size.label.x,vjust=-80),
                       axis.text.x  = element_text(size=size.text.x, vjust = 0.5),
                       axis.title.y = element_text(size=size.label.x, vjust = 2),
                       axis.text.y  = element_text(size=size.text.x),
                       legend.title = element_blank(),
                       legend.text = element_text(size = size.legend.text),
                       legend.position = c(0.15, 0.95),
                       legend.key = element_rect(fill = "white")) 
    
    ## create the data.frame
    
    pop_df = data.frame(abund = c(vettad[1:year],
                                  escapement_obs[1:year],
                                  escapement_real[1:year],
                                  harvest[1:year]), 
                        type = rep(c("ocean","escap_obs","escap_real","harvest"),each = year), year = 1:year) %>%
      filter(.,year>=10)
    
    
    if (year > 10) {
      pop_gg = ggplot(dplyr::filter(pop_df, type %in% c("ocean","escap_obs","escap_real","harvest")), aes(x = year, y = abund, group = type, shape = type, linetype = type)) +
        geom_point(alpha = 0.4, size = size.point) +
        geom_line() +
        theme.pop +
        guides(size = guide_legend(override.aes = list(alpha = 0.2))) +
        scale_y_continuous(limits = c(0,(max(pop_df$abund) + max(pop_df$abund)*0.3))) +
        scale_x_continuous(limits = c(0,year+2)) +
        labs(y = "#") +
        labs(x = "Year") 
      
    } else {
      pop_gg = "Not enough simulation years for plot (must be > 10)"} 
    
  
    ### For the exticnt populations, I want only the vectors up to year of exinction ####
    
    harvest_mean_size = harvest_mean_size[1:year]
    harvest_total_weight = harvest_total_weight[1:year]
    harvest_mean_age = harvest_mean_age[1:year]
    escapement_obs = escapement_obs[1:year]
    escapement_real = escapement_real[1:year]
    harvest = harvest[1:year]
    test.pop = test.pop[1:year]
    stock_v = stock_v[1:year]
    harvest_mean = mean(harvest[10:length(harvest)],na.rm = T)
    harvest_sd = sd(harvest[10:length(harvest)],na.rm = T)
    harvest_cv = harvest_sd/harvest_mean  
    
    ris.list = list("extinct"=extinct, 
                    "yearextinct" = i-1,
                    "ocean_size"=vettad,
                    #"test.pop" = test.pop,
                    "escapement_real" = escapement_real,
                    "escapement_obs" = escapement_obs,
                    "plot_pop" = pop_gg,
                    "alpha_stock" = alpha_stock,
                    "Rp_stock" = Rp_stock,
                    "recruit" = recruitment,
                    "harvest" = harvest,
                    "harvest_weight" = harvest_total_weight,
                    "harvest_size" = harvest_mean_size,
                    "harvest_age" = harvest_mean_age,
                    "harvest_mean_num" = harvest_mean,
                    "harvest_sd_num" = harvest_sd,
                    "harvest_cv_num" = harvest_cv,
                    "escapement_size" = escapement_mean_size,
                    "escapement_age" = escapement_mean_age,
                    "pop" = pop_df,
                    "stock" = stock_v)
    
    return(ris.list)
    
    dev.off()
    
    #########other potential items to be returned######
    
    
  }# close all