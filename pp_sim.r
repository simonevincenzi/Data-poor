
# check if required pacakges are present, otherwise install them

list.of.packages <- c("Rlab", "MASS", "parallel", "tidyverse","MCMCpack")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rlab)
library(MASS)
library(parallel)
library(tidyverse)
library(MCMCpack)

################################


pp_sim.f <- function (S = 1000,N = 10000,iter = 150,num.loci =5,
                         plotorno = 1,num.alleles = 2,sd.alleles = 0.1,recomb=1, mort_r = 0.2,mort_r_harv = 0.27, mat_prop = c(0.05,0.5,0.93,1), age_repr = c(2,3,4,5),alpha_stock = 50,Rp_stock = 5000)
  
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
  
  ### I define the standard deviation of distribution of single alleles, which depends on heritability at time 1 and environmental variance this can be adapted to other traits that have means != 0 the correct way is to proceed from a phenotypic variance and heritability and then obtain the starting additive genetic variance
  
  heritability = 0.3
  
  mort_r_harv_r = runif(iter,min = (mort_r_harv - mort_r_harv*0.5), max = mort_r_harv + (mort_r_harv*0.5))
  mort_r_r = runif(iter,min = (mort_r - mort_r*0.5), max = mort_r + (mort_r*0.5))
  
  starting_genetic_variance = (heritability * var.env)/(1-heritability)
  
  # sd.alleles is the standard deviation of the normal distribution of allelic effects 
  
  sd.alleles = sqrt(starting_genetic_variance/(2*num.loci))
  
  ############ Define values of some model parameters not included in the function
  
  
  ## the population stays on a matrix (area.pop), each column is an individual. Here below I define
  # the position of entities in the column
  
  pheno.pos = 1  # this means that the first row is the phenotype of the individuals
  age.pos = 2  # second row is the age of the individual
  matur.pos = 3 # third row is maturity status
  
  length_prop = 3 #  length of entities of individual vector excluding loci (phenotipic value, age, maturity status are the top rows)
  
  heteroz.mat = matrix(0,num.loci*num.alleles,(iter%/%50))  ##### define the matrix for heterozigosyis, same as all.freq100 but with a col less
  
  area.pop <- matrix(0,length_prop+(num.loci*2),S) ##  matrix space area.pop of S elements, 
  # for each col the first row is the pheno value, second is age, third is maturity status 
  #the other num.loci*2 rows are the values for each gene
  #the first num.loci (2-num.loci+1)for the first chromosome, the second num.loci (num.loci+2-(2*num.loci)+2) for the second chromosome.
  
  fchrom = seq((length_prop+1),num.loci+length_prop,1) ##identifier of rows for the 1st chromosome
  schrom = seq(num.loci+(length_prop+1),num.loci*2+(length_prop),1) ##identifier of rows for the 2nd chromosome
  
  
  ## define initial ages for the S individuals in the population at time = 1
  
  initial_ages = round(runif(S,min = 0, max = 5))
  
  
  area.pop[age.pos,] = initial_ages # assign ages to individuals
  
  #####CREATE THE MATRIX OF ALLELES IN THE POPULATION##########
  
  alleles.mat = matrix(rnorm(num.loci*num.alleles,0,sd.alleles), 
                       num.loci,num.alleles)  
  
  #in matrix alleles.mat I extract randomly, for each of the num.loci, the num.alleles number of alleles 
  
  alleles.mat.keep = alleles.mat  
  
  alleles.mat.gen = matrix(0, num.loci,num.alleles)
  
  for (nl in 1:num.loci) {
    
    alleles.mat.gen[nl,] = seq(1:num.alleles)
    
  } ###I assign a number from 1 to num.alleles to each allele (for Fst since I need the genotype)
  

  #double the matrix alleles.mat for ease of extraction of alleles for the diploid organism 
  
  doublealleles.mat = rbind(alleles.mat,alleles.mat)
  
  ##########assign alleles to individuals at the start of simulation
  
  for (a in 1:(num.loci*2)) { #assign alleles to loci in the population matrix area.pop
    
    area.pop[length_prop+a,] = sample(doublealleles.mat[a,],ncol(area.pop),replace=T) } #end assignment 
  
  area.pop[pheno.pos,] = colSums(area.pop[(length_prop+1):nrow(area.pop),]) #the first row of the population matrix is the sum of the genetic values of each alleles 
  #of an individual
  ############################################ 
  pheno = area.pop[pheno.pos,] + rnorm(ncol(area.pop),0,sqrt(var.env)) 
  # phenotypic value = genotypic value + environmental deviate from N(0,sqrt(var.env)) (all stored in pheno)
  
  
  #####create vectors and matrix for later use
  mediaphen = rep(0,iter) # vector with mean of the phenotype (one value each year)
  sdphen = rep(0,iter)  # vector with standard deviation of the phenotype (one value each year)
  mean.age = rep(0,iter) # vector with mean age at each time step
  recruitment = rep(0,iter) # vector with escapement for each year
  escapement = rep(0,iter) # vector of escapement 
  
  vettad = rep(0,iter)    #######population size post_mortality at each time step
  vettad_pre = rep(0,iter)  ###### population size pre_mortality at each time step
  addgenvar = rep(0,iter)   ####additive genetic variance of the whole population at each time step
  addgenmean = rep(0,iter)  ### mean genetic value in the whole population at each time step
  phenomean = rep(0,iter)  #### mean of the phenotype pre_selection, including the newborns of the year i-1
  mean.var.age = matrix(0,iter,2)  ## matrix with iter row and 2 cols (1 = mean age, 2 = sd of age every year)
  
  w_gen_mean = rep(0,iter)  ###### mean of the phenotype in a new generation at each time step
  w_gen_var = rep(0,iter)		###### sampling variance of the phenotype in a new generation at each time step
  
  num.sons.check = matrix(0,2,iter) #check the number of offspring recruited in the population
  
  test.pop = list() # list of area matrix
  harvest = rep(0,iter)
  
  ###########
  
  all.freq100 = matrix(0,num.loci*num.alleles,(iter%/%50)+3) ## initialize matrix with 
  #frequencies of alleles during simulation, every 50 years I record the value
  dim(alleles.mat) = c(num.loci*num.alleles,1) ##unravel matrix with alleles (alleles.mat) in order to have a vector
  
  heter_mat_year = matrix(0,num.loci*num.alleles,iter) # saves allelic frequency each year
  heter_vect_mean_year = rep(0,iter) #saves mean heterozygosity each year
  heter_vect_sd_year = rep(0,iter) #saves sd of heterozygosity each year
  
  col.freq = 2 # counter for column number when checking for frequencies of alleles at different time steps  
  
  a = 0
  
  ##################### START FOR TIME LOOP ########################
  
  
  for(i in 1:iter)  { ##time-steps or cohorts or year
    
    print(i)  # prints time-steps
    
    ###### Compute mean and variance of the phenotype for the first generation
    
    if (i==1) {
      w_gen_mean[i] = mean(pheno)
      w_gen_var[i]  = var(pheno)
    }
    
    #####################################
    
    #########check if the population went extinct####
    
    
    if(length(area.pop[1,area.pop[1,]!=0])<2){  ##is population extinct?  
      print(paste("Extinct at year", i)) #prints the year of extinction
      extinct = 1 # when an individual is dead the first column is 0
      
      if(iter>=50) {
        
        all.freq100[,dim(all.freq100)[2]] = all.freq100[,col.freq] - all.freq100[,2]  #in case of population going extinct, it fills the 
        #last col of all.freq100 with the difference between 
        heteroz = 0		
      }#frequency at time 50 and frequency at the last year to be divided
      #by 50 before extinction
      break}   #close for to check for extinction
    
    ######################################  
    
    ### here below I create the matrix of allele frequencies every 50 years starting from year 1
    
    if (i==1) { #if first year of the simulation
      
      all.freq100[,1] = alleles.mat  #first column is the vector of alleles
      
      for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
        
        all.freq100[zz,col.freq] = length(area.pop[area.pop==alleles.mat[zz]])/
          (2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
        #population by 2 times the number of individuals alive
        # (diploid organism)
      }
      
    }
    
    
    if (i%%50==0) {  #if year is a multiple of 50
      
      col.freq=col.freq+1  #update counter of rows

      for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
        
        all.freq100[zz,col.freq] = length(area.pop[area.pop==alleles.mat[zz]])/(2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
        #population by 2 times the number of individuals alive
      }
      
    }
    
    
    
    if (i==iter & iter>=50){ #if last year of the simulation
      
      all.freq100[,dim(all.freq100)[2]] = all.freq100[,(dim(all.freq100)[2]-1)] - all.freq100[,2] #it fills the 
      #last col of all.freq100 with the difference between 
      #allelic frequency at time 50 and frequency at the last year to be divided
      #by 50 before extinction
      
      heteroz.mat = all.freq100[,(dim(all.freq100)[2]-1)]   #matrix for heterozigosity (no mutants included)
      
      
      
    }
    
    a = a+1   ## counter for the creation of matrix of allelic values
    
    
    #### allelic frequency each year
    
    for (zz in (1:(num.loci*num.alleles))) {  #for each alleles in the population (present at time 1)
      
      
      heter_mat_year[zz,i] = length(area.pop[area.pop==alleles.mat[zz]])/(2*length(area.pop[1,area.pop[1,]!=0])) #divide the number of alleles in the 
      #population by 2 times the number of individuals alive
      #print(i)																							# (diploid organism)
    }
    
    
    
    ##### close before and after climate change #######
    
    
    ##### compute mean age of ocean salmon ######
    
    mean.age.who = which(area.pop[age.pos,]>=2 & area.pop[1,]!=0)
    
    mean.age[i] = mean(area.pop[2,mean.age.who],na.rm = T)
    
    
    ############ increase age by 1 if t>1 ##########
    
    if (i>1) {
      area.pop[age.pos,area.pop[pheno.pos,]!=0] = area.pop[age.pos,area.pop[pheno.pos,]!=0] + 1 
    }
    ###########
    
    
    ####vector of optimum and mean and sd of phenotype in the population ######
    
    
    mediaphen[i] = mean(area.pop[pheno.pos,area.pop[pheno.pos,]!=0])
    
    sdphen[i] = sd(area.pop[pheno.pos,area.pop[pheno.pos,]!=0])
    
    recruitment[i] = recruit
    
    #### close vector of optimum and mean and sd of phenotype in the population ######
  
    
    addgenvar[i] <- var(area.pop[pheno.pos,area.pop[pheno.pos,]!=0])
    addgenmean[i] <- mean(area.pop[pheno.pos,area.pop[pheno.pos,]!=0])
    phenomean[i] <- mean(pheno[area.pop[pheno.pos,]!=0])
    varphen[i] = var(pheno[pheno!=0])
    
    
    vettad_pre[i] = length(area.pop[pheno.pos,area.pop[pheno.pos,]!=0]) ########population size pre_selection
    
    ############### MEAN AND SD OF AGE #####
    
    who = which(area.pop[pheno.pos,]!=0)  # where the individuals are (because empty columns have zero in first position)
    mean.var.age[i,1] = mean(area.pop[age.pos,who]) #mean age
    mean.var.age[i,2] = sd(area.pop[age.pos,who]) #sd of age
    
    
    
    ############mortality##############
    
  
    who = which(area.pop[pheno.pos,]!=0) # where the individuals are (because empty columns have zero in first position) 

    vettad_pre[i] = length(area.pop[1,area.pop[pheno.pos,]!=0]) ########population size pre_mortality including juveniles
    
    ## mortality is sequential, natural and then harvest
    
    ### random number for each adult 
    
    f = rep(0,S) 
    p = rep(0,S)
    p[who] = runif(length(p[who]))
    r = 1 - (exp(-mort_r_r[i]))
    dead <- which(p<r)   #### since the p of absent columns are also empty, they are always "dead"
    area.pop[,who[dead]] <- 0 ##remove who dies, ie all columns gets 0
    pheno[dead] <- 0  #phenotypes of the dead (thus empty cols) are 0

    ### harvest after natural mortality
    dead = 0
    p = rep(0,S)
    who = which(area.pop[pheno.pos,]!=0)
    empty_b = length(which(area.pop[pheno.pos,]==0))
    p[who] = runif(length(p[who]))
    r_h = 1 - (exp(-mort_r_harv_r[i]))
    dead <- which(p<r_h)   #### since the p of absent columns are also empty, they are always "dead"
    if (length(dead) > 0) {
      area.pop[,dead] <- 0 ##remove who dies, ie all columns gets 0
    pheno[dead] <- 0 }  #phenotypes of the dead (thus empty cols) are 0
    empty_a = length(which(area.pop[pheno.pos,]==0))
    harvest[i] = empty_a - empty_b
    ####
    
    vettad[i] = length(area.pop[pheno.pos,(area.pop[pheno.pos,]!=0 & area.pop[age.pos,]>=2)])
    
    ########population size post_mortality includig juveniles  
    
    
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
    
    
    nsonscont = 0
    
    cosons = 0
    
    vettsons = rep(100,S*2)  # define with 100s and then delete those 100s when manipulating the vector
    #matsons = rbind(rep(100,Rp_stock*2),matrix(0,(length_prop+(num.loci*2)-1),Rp_stock*2))# define with 100s the first and then delete the remaining 100s 
    
    matsons = matrix(0,(length_prop+(num.loci*2)),Rp_stock*2)
    matsons[1,] = 100
    
    #when manipulating the vector
    pheno_parents = rep(100,S*2)  ##empty vector that will contain the mean phenotype of parents
    #sex.from.parents = sample(c(0,1),S*2,replace=T)
    #matsons[(sex),] = sex.from.parents 
    matsons[age.pos,] = 0  # age of newborns is 0, so next year they will be 0 
    
    reproductor = which(area.pop[matur.pos,] == 1)       ##### vector of position in the matrix of reproductors
    stock = length(reproductor)  # stock
    
    recruit = round(alpha_stock*stock*exp(-alpha_stock*(stock/(Rp_stock*exp(1))))) # recruitment

    # print(recruit)
    
    escapement[i] = stock

    test.pop[[i]] = area.pop
    
    
    ## Sexual reproduction
    
    
    gene.from.parents = rbind(rep(100,Rp_stock*2),matrix(0,(num.loci*2*2)-1,Rp_stock*2))
    
    number.of.sons = rep(100,S*2)
    
    if (recruit >1 & stock >=2)   { #open if for reproduction
      
      
      ######### open control for odd number of reproductors
      
      if (length(reproductor)%%2 != 0) {		#open and close control for odd number of reproductors, since I assume that a female mates with only one male
        
        reproductor <- reproductor[-sample(1:length(reproductor),1)]
      }  # close for odd number of reproductors 
      
      
      ############## close control for odd number of reproductors #########
      
      repro.matrix <- matrix(sample(reproductor,length(reproductor)),2,length(reproductor)/2) 
      ## produce the matrix of reproductors, each column is a couple 
      
      ## only 20% of pairs reproduce
      
      if(length(reproductor) > 5) { ## just by chance it can happen that when numbers are low, we then have zero reproductors
      who_repr = rbern(n = length(reproductor)/2, prob = 0.2) 
      if (sum(who_repr) == 0) {# at least 2 reproduce
      who_repr[sample(1:length(who_repr),1)] = 1
      }
      repro.matrix = repro.matrix[,which(who_repr==1),drop=FALSE] # keep the matrix!
      }
      rel_repr = rgamma(n = ncol(repro.matrix), alpha = 1, beta = 4) # gamma distribution of recruit per pair
      if (sum(rel_repr) == 0) {# at least 2 reproduce
        rel_repr[sample(1:length(rel_repr),1)] = 1
      }
      # print(ncol(repro.matrix))
      sons_rel = round((rel_repr/sum(rel_repr)) * recruit) # number of recruit per pair
      # at least one fish is produced
      # print(sons_rel)
      ######### begin loop for reproductors ###########
      
      
      for(rr in 1:ncol(repro.matrix)) {   #begin loop for of reproduction, go on for all the couples
        
        nsonscouple <- sons_rel[rr]   ### this is the number of newborns for a couple, 
        #it is a deviate from a poisson distribution with sons.mean mean and variance
        
        
        if(nsonscouple>0) {  # open if for couples with at least one offspring
          
          
          
          #matsons[sector,(nsonscont+1):(nsonscont + nsonscouple)] = sector.from.parents #assign the sector to offspring
          
          
          #nsonscont = nsonscont + nsonscouple
          
          gene.from.father = matrix(0,num.loci*2,nsonscouple)  # alleles coming from father (empty) used in case of recombination
          gene.from.mother = matrix(0,num.loci*2,nsonscouple) # alleles coming from mother (empty) used in case of recombination
          
          gene.sons = matrix(0,num.loci*2,nsonscouple)  #create a matrix with cols = number of offspring of the couple and rows = num.loci * 2
          
          ####### open recombination####
          #close for recombination = 0
          if (recomb==0) {  #whether to recombine or not
            for (aa in 1:nsonscouple) {   #cicle for number of offspring
              
              sampleseq = sample(c(1,2),2,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              
              if (sampleseq[1]==1) {gene.sons[1:num.loci,aa] = area.pop[fchrom,repro.matrix[1,rr]]}  ## if 1 for father I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              else if (sampleseq[1]==2) {gene.sons[1:num.loci,aa] = area.pop[schrom,repro.matrix[1,rr]]} # if 2 for father I pick the second chromosome
              
              if (sampleseq[2]==1) {gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[fchrom,repro.matrix[2,rr]]} ## if 1 for mother I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              else if (sampleseq[2]==2) {gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[schrom,repro.matrix[2,rr]]} # if 2 for mother I pick the second chromosome
              
            }
          }  else if (recomb==1) { #### recombination ==1
            
            
            for (aa in 1:nsonscouple) {   #cicle for number of offspring
              
              rchrom.f = 0
              rchrom.m = 0
              
              
              sampleseq = sample(c(1,2),num.loci,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              rchrom.f = diag(cbind(fchrom,schrom)[,sampleseq])
              gene.sons[1:num.loci,aa] = area.pop[rchrom.f,repro.matrix[1,rr]]  ## if 1 for father I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              
              
              gene.from.father[,aa] = area.pop[rchrom.f,repro.matrix[1,rr]] 
              
              
              sampleseq = sample(c(1,2),num.loci,replace=T)    #draw the chromosome to be passed to offspring (First for father, second for mother)
              rchrom.m = diag(cbind(fchrom,schrom)[,sampleseq])
              gene.sons[(num.loci+1):(num.loci*2),aa] = area.pop[rchrom.m,repro.matrix[2,rr]] ## if 1 for mother I pick the first chromosome (fchrom = 2:num.loci+1, schrom = num.loci + 2:max) 
              # if 2 for mother I pick the second chromosome
              
              gene.from.mother[,aa] = area.pop[rchrom.m,repro.matrix[2,rr]]
              
            }
          }  
          
        
          ##### close recombination options #########
          
          ########
          
           matsons[((length_prop+1):nrow(matsons)),(cosons+1):(cosons + nsonscouple)] = gene.sons
          
          matsons[1,(cosons+1):(cosons + nsonscouple)] <- 0
          pheno_parents[(cosons+1):(cosons + nsonscouple)] = mean(pheno[repro.matrix[1,rr]],pheno[repro.matrix[1,rr]])
          
          gene.from.parents[,(cosons+1):(cosons + nsonscouple)] = rbind(gene.from.father,gene.from.mother)
          
          number.of.sons[(cosons+1):(cosons + nsonscouple)] = nsonscouple
          
          cosons <- length(matsons[1,matsons[1,]!=100])      
          
        }  # close if  for couples with at least one offspring
      }   #close loop for reproduction in one sector
      
    } ### close for if (length(reproductor) >1)
    
    area.pop[,reproductor] = 0  ## all escapees die
    
    empty = which(area.pop[pheno.pos,]==0)  ### identify the available spots for newborns in the matrix
  
    
    ##### open if there are offspring produced ########
    
    
    if (length(matsons[1,matsons[1,]!=100]) >1) ##if there are offspring produced
    { 
      
      #print(1.5)
      matsons = matsons[,matsons[1,]!=100]  ## delete the empty places in matsons
      pheno_parents = pheno_parents[pheno_parents!=100]  ##delete the empty places in pheno_parents
      number.of.sons = number.of.sons[number.of.sons!=100]  ##delete the empty places in number.of.sons
      gene.from.parents = gene.from.parents[,gene.from.parents[1,]!=100] ##delete the empty places in gene.from.parents
      

      
      sons = 1  ## this is related to if (length(matsons[1,matsons[1,]!=100]) >1)
    } else # if there are no newborns produced
    { matsons[,] == 0  
      sons=0
    }
    
    ####### end of if (length(matsons[1,matsons[1,]!=100]) >1)  (including the else)
    
    
    
    
    vettpheno = rep(50,S*2) #### I prepare the empty vector of newborn phenotypes (50 means there is nothing)
    
    
    #### vector of offspring is empty if sons == 0 or no space is available in the matrix #####
    
    if(sons==0 | length(empty)==0){
      vettpheno=0
      ####### close of offspring is empty if sons == 0 or no space is available #####
    } else if (sons==1){
      
      if(length(empty)>=ncol(matsons) ){
        
        area.pop[1:nrow(area.pop),empty[1:ncol(matsons)]] = matsons
        
        area.pop[1,empty[1:ncol(matsons)]] = colSums(area.pop[(length_prop+1):nrow(area.pop),
                                                                empty[1:ncol(matsons)]])
        vettpheno[1:ncol(matsons)] <- area.pop[1,empty[1:ncol(matsons)]] + 
          rnorm(ncol(matsons),0,sqrt(var.env)) 
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty[1:ncol(matsons)]] <- vettpheno
        
        pheno_parents_unique = unique(pheno_parents)
        
        off.mean = rep(0,length(pheno_parents_unique))
        
        sons.for.couple = rep(0,length(pheno_parents_unique))
        
        
        
        for (hh in 1:length(pheno_parents_unique)) {
          off.mean[hh] = mean(vettpheno[pheno_parents == pheno_parents_unique[hh]]) 
        
        sons.for.couple[hh] = min(number.of.sons[pheno_parents == pheno_parents_unique[hh]])
        
        
        if (i==50) {
          
          check.parents  = rbind(gene.from.parents,matsons)
          
          
        }
        
        
        }
        
        
      
        
        num.sons.check[1,i] = 1   ##this tells me that there were more empty spots than offspring (when 0 is more offspring than spots)
        num.sons.check[2,i] = length(vettpheno) ##this tells me how many offspring were introduced
        #herit.coef[i] = summary(lm(vettpheno  ~ pheno_parents))$coef[2]
        #print(10)
        
        
        
      } else if(length(empty)==1){   #this is special case because I cannot use colSums with only one column
  
        area.pop[1:nrow(area.pop),empty]<-matsons[,1:length(empty)]
        area.pop[1,empty] = sum(area.pop[(length_prop+1):nrow(area.pop),empty])
        
        vettpheno[1:length(empty)] = area.pop[1,empty] + rnorm(length(empty),0,sqrt(var.env))
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty]<-vettpheno
        
      } else if(length(empty)<= ncol(matsons)){
        
        
        
        area.pop[1:nrow(area.pop),empty]<-matsons[,1:length(empty)]
        area.pop[1,empty] = colSums(area.pop[(length_prop+1):nrow(area.pop),empty])
        
        vettpheno[1:length(empty)] = area.pop[1,empty] + rnorm(length(empty),0,sqrt(var.env))
        
        vettpheno = vettpheno[vettpheno!=50]
        
        pheno[empty]<-vettpheno
        
        pheno_parents_unique = unique(pheno_parents[1:length(empty)])
        
        off.mean = rep(0,length(pheno_parents_unique))
        sons.for.couple = rep(0,length(pheno_parents_unique))
        
        
        
        ##############
        
        num.sons.check[2,i] = length(vettpheno) ##this tells me how many offspring were introduced
        #herit.coef[i] = summary(lm(vettpheno  ~ pheno_parents[1:length(empty)]))$coef[2]
        #print(13)
        
        
        
      }
      
    }
    
    
    ####here I compute the mean and the variance of the new generation  phenotypes####
    if (i>1 & sons==1 & length(empty)>0) {
      w_gen_mean[i] = mean(vettpheno)
      w_gen_var[i]  = var(vettpheno)
    }
    
      
  year = i  # useful for plot  
    
    ####################   CLOSE REPRODUCTION ##########################
    
  }  ##close time
  
  ### heteroz  vector for each locus at the end of simulation time
  if (i == iter & iter >=50 & extinct ==0) {
    
    dim(heteroz.mat) =  c(num.loci,num.alleles)
    heteroz = matrix(0,num.loci,1)
    heteroz = 1 - rowSums(heteroz.mat^2)
  } else {heteroz = NULL}
  
  ### mean and sd of heteroz vector for each year of simulation time
  
  for (hh in 1:i) {
    
    all_freq_year = heter_mat_year[,hh]
    dim(all_freq_year) =  c(num.loci,num.alleles)
    heter_vect_mean_year[hh] = mean((1- rowSums(all_freq_year^2)), na.rm =T)
    heter_vect_sd_year[hh] = sd((1- rowSums(all_freq_year^2)), na.rm =T)
  } 
  
  ##################### CLOSE FOR TIME LOOP ########################
  
  
  
  ### PLOT
  
  size.title = 15
  line.lwd = 1.2
  size.label.x = 18
  size.text.x = 14
  size.point = 4
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
  
  pop_df = data.frame(abund = c(vettad,escapement,harvest), type = rep(c("ocean","escap","harvest"),each = year), year = 1:year) %>%
    filter(.,year>=10)
  
  
  pop_gg = ggplot(pop_df, aes(x = year, y = abund, group = type, shape = type, linetype = type)) +
    geom_point(alpha = 0.4) +
    geom_line() +
    theme.pop +
    guides(size = guide_legend(override.aes = list(alpha = 0.2))) +
    scale_y_continuous(limits = c(0,(max(pop_df$abund) + 1000))) +
    scale_x_continuous(limits = c(0,year+2)) +
    labs(y = "#") +
    labs(x = "Year") 
  
  
  ### change all.freq100 to data.frame #####

  
  
    
  all.freq100 = as.data.frame(all.freq100)
  steps = c(1,seq(0,iter,50)[-1])
  #steps[length(steps)] = min(i,steps[length(steps)])  ### if the population goes extinct
  colnames(all.freq100)[1] = "Allelic_value"
  colnames(all.freq100)[2:(length(steps)+1)] = paste("Freq_time",steps,sep ="_")
  colnames(all.freq100)[ncol(all.freq100)] = "Freq_diff"
  all.freq100$Gene = rep(0,nrow(all.freq100))
  
  for (j in 1:num.loci) {
    all.freq100$Gene[seq(j,(num.loci*num.alleles),num.loci)] = j
    
    
  }
  ### Gene 0 is the mutant alleles
  
  all.freq100 = arrange(all.freq100, Gene)
  
  ris.list = list("extinct"=extinct, 
                  "yearextinct" = i-1,
                  "phenomean" = phenomean,
                  "heteroz" = heteroz,
                  "numloci"=num.loci,
                  "num.alleles"=num.alleles,
                  "sd.alleles" = sd.alleles, 
                  "recomb"= recomb, 
                  "fitness.mean" = fitness.mean, 
                  "popsize.post"=vettad,
                  "all.freq" = all.freq100,
                  "heter.mean.year" = heter_vect_mean_year,
                  "heter.sd.year" = heter_vect_sd_year,
                  "test.pop" = test.pop,
                  "escap" = escapement,
                  "plot" = pop_gg,
                  "harvest" = harvest,
                  "pop" = pop_df)
  
  return(ris.list)
  
  dev.off()
  
  #########other potential items to be returned######
  
  
}# close all

#### The function returns a list of results that can be accessed with $
# extinct = 1 if a population wen extinct, 0 otherwise
# addvar = vector of additive genetic variance, one value for each year
# phenvar = vector of phenotypic variance, one value for each year 
# yearextinct = year of exinction if the population went extinct, iter-1 otherwise
# phenomean = vector of mean phenotype, one value for each year
# heteroz = vector of mean heterozigosity, one value for each year
# sdopt = input sdopt
# meanopt = input meanopt
# numloci = input num.loci
# num.alleles = input num.alleles,
# sd.alleles = input sd.alleles
# recomb = input recomb 
# major = input major
# mumut = input mu.mut
# mutalfa = input mut.alfa, 
# fitness.mean = vector of mean fitness, one value for each year 
# fitness.variance = vector of fitness variance, one value for each year 
# optimum = vector of realized climate variable values, one value for each year
# catastr.vett = vector of point extremes, 1 if the point extreme occurred, 0 otherwise
# popsize.post = population size (after mortality, both normal and due to point extreme), one value for eacy year
# cat.freq = input p.post 
# age.sex.mat = input age.sex.mat
# avg.surv = input avg.surv
# sons.mean = input sons.mean
# cor.all.est = correlation (Pearson's r) between allele frequency of alleles at time = 1 and allele frequency 
##    at the end of simulation time
# cor.all.p" = p-value of the correlation
# n.mut = number of mutant alleles in the population at the end of simulation time


##### DATA

# S.vett = 500
# N.vett = 500
# iter.vett = 250
# selection.vett = c(0.08,0.11,0.13,0.15)
# #selection.vett = c(0.15)
# sdopt.vett = c(0,0.01,0.015)
# num.loci.vett = 10
# plotorno.vett = 0
# meanopt.vett = c(0,0.015)
# major.vett = 0
# num.alleles.vett = 10
# sd.alleles.vett = 0.05
# recomb.vett = 0
# mutation.vett = 0
# mut.alfa.vett = 0.3
# catastr.vett = 1
# p.pre.vett = 0.05
# p.post.vett = c(0.05,0.1,0.15)
# catastr.mort.vett = c(0.3,0.5,0.7)
# age.sex.mat.vett = c(1,2,3,4)
# avg.surv.vett = c(0.7)
# sons.mean.vett = c(1,1.5,2,2.5,3,3.5)    ##### with 1 every population goes extinct
# first.phase.opt.vett = 100
# var.amb.init.vett = 1
# 
# 
# ###### remember that the columns must be named!!!
# 
# dataforpar = expand.grid(S = S.vett, N = N.vett ,iter = iter.vett , selection = selection.vett ,sdopt = sdopt.vett ,num.loci = num.loci.vett,plotorno = plotorno.vett,meanopt = meanopt.vett,major = major.vett,num.alleles = num.alleles.vett,sd.alleles = sd.alleles.vett ,recomb = recomb.vett ,mutation = mutation.vett,mut.alfa = mut.alfa.vett,catastr = catastr.vett,p.pre = p.pre.vett,p.post = p.post.vett,catastr.mort = catastr.mort.vett,
#                          age.sex.mat = age.sex.mat.vett,avg.surv = avg.surv.vett,sons.mean = sons.mean.vett,first.phase = first.phase.opt.vett,var.amb = var.amb.init.vett) 
# 
# 
# 
# replicate.sim = 10
# 
# dataforpar = do.call(rbind, replicate(replicate.sim, dataforpar, simplify=FALSE))