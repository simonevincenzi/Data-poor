# We begin by specifying the names of four necessary data files that contain the following information:
#   
# 1. observed total number of adult spawners (escapement) by year;
# 2. observed age composition of adult spawners by year;
# 3. observed total harvest by year;
# 4. hatchery releases by year.
# 5. Let’s also define the following parameters, which will be referenced throughout the analysis.
# 
# n_yrs: number of calendar years of data
# A: number of age classes
# M: number of covariates

## 1. file with escapement data
## [n_yrs x 2] matrix of obs counts; 1st col is calendar yr
fn_esc <- "SkagitSthdEsc.csv"

## 2. file with age comp data
## [n_yrs x (1+A)]; 1st col is calendar yr
fn_age <- "SkagitSthdAge.csv"
## min & max ages
age_min <- 3
age_max <- 8
## years, if any, of age-comp to skip; see below
age_skip <- 2

## 3. file with harvest data
## [n_yrs x 2] matrix of obs catch; 1st col is calendar yr
fn_harv <- "SkagitSthdCatch.csv"

## 4. file with covariate data
## [n_yrs x (1+MM)]; 1st col is calendar yr
fn_hrel <- "SkagitHatchRel.csv"

## time lags (years) for covariates
flow_lag <- 1
marine_lag <- 2
hrel_lag <- 2

## number of years of forecasts
n_fore <- 1

## file where to save JAGS model
fn_jags <- "SkagitSthd_RR_JAGS.txt"

## upper threshold for Gelman & Rubin's (1992) potential scale reduction factor (Rhat).
Rhat_thresh <- 1.1

## URL for example data files
## set to NULL if using a local folder/directory
ex_url <- "https://raw.githubusercontent.com/mdscheuerell/ASSESSOR/master/"


# Loading the fish data
# Here we load in the first three data files and do some simple calculations and manipulations. First the spawner data:
  

## escapement
dat_esc <- read.csv(text = getURL(paste0(ex_url,fn_esc)))
## years of data
dat_yrs <- dat_esc$year
## number of years of data
n_yrs <- length(dat_yrs)
## get first & last years
yr_frst <- min(dat_yrs)
yr_last <- max(dat_yrs)
## log of escapement
ln_dat_esc <- c(log(dat_esc[,-1]),rep(NA,n_fore))


# Next the age composition data:
  
## age comp data
dat_age <- read.csv(text = getURL(paste0(ex_url,fn_age)))
## drop year col & first age_min+age_skip rows
dat_age <- dat_age[-(1:(age_min+age_skip)),-1]
## num of age classes
A <- age_max-age_min+1
## add row(s) of NA's for forecast years
dat_age <- rbind(dat_age,matrix(0,n_fore,A,dimnames=list(n_yrs+seq(n_fore),colnames(dat_age))))
## total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age,1,sum)
## row indices for any years with no obs age comp
idx_NA_yrs <- which(dat_age$sum<A,TRUE)
## replace 0's in yrs w/o any obs with NA's
dat_age[idx_NA_yrs,(1:A)] <- NA
## change total in yrs w/o any obs from 0 to A to help dmulti()
dat_age[idx_NA_yrs,"sum"] <- A
## convert class
dat_age <- as.matrix(dat_age)

# And then the harvest data:
  
## harvest
dat_harv <- read.csv(text = getURL(paste0(ex_url,fn_harv)))

## 1/5 of escapement

# dat_harv$catch = round(dat_esc$escapement/5)
# dat_harv[which(is.na(dat_harv$catch)),"catch"] = 500


## drop year col & first age_max rows
dat_harv <- c(dat_harv[,-1],rep(0,n_fore))

#dat_harv = 0
#dat_harv = rep(0,length(dat_harv))

dat_harv[10] = NA
