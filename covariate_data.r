# Loading the covariates
# This analysis uses 6 covariates as drivers of the population’s instrinic growth rate:
#   
# 1. Maximum river discharge in winter
# 2. Minimum river discharge in summer
# 3. Pacific Decadal Oscillation (PDO)
# 4. North Pacific Gyre Oscillation (NPGO)
# 5. Hatchery releases
# 6. River discharge
# 
# We begin by getting the daily flow data from the US Geological Service National Water Information System. We will use the direct link to the gage data from the Skagit River near Mount Vernon, WA (#12178100), beginning with the first year of fish data.
  
  ## flow site
  flow_site <- 12178100
  ## get URL for flow data from USGS
  flow_url <- paste0("http://waterdata.usgs.gov/nwis/dv",
                     "?cb_00060=on",
                     "&format=rdb",
                     "&site_no=",flow_site,
                     "&begin_date=",yr_frst,"-01-01")
  
  flow_url = ("https://waterdata.usgs.gov/nwis/dv?referred_module=sw&cb_00060=on&format=rdb&site_no=12178100&begin_date=1978-01-01") # the one above does not work

  # Next we will retrieve the raw data file and print its metadata.
  
  ## raw flow data from USGS
  flow_raw <- readLines(flow_url)
  ## lines with metadata
  hdr_flow <- which(lapply(flow_raw,grep,pattern="\\#")==1, arr.ind=TRUE)
  ## print flow metadata
  print(flow_raw[hdr_flow],quote=FALSE)
  
  
# Lastly, we will extract the actual flow data for the years of interest and inspect the file contents.
  
  ## flow data for years of interest
  dat_flow <-  read.csv(textConnection(flow_raw[-c(hdr_flow,max(hdr_flow+2))]),
                        header=TRUE, stringsAsFactors=FALSE, sep="\t")
  colnames(dat_flow) <- unlist(strsplit(tolower(flow_raw[max(hdr_flow)+1]), split="\\s+"))
  head(dat_flow)
  
  
  # The first 3 columns in the data file are the agency (agency_cd), location (site_no), and date (datetime). The daily flow measurements are in the 4th column (149343_00060_00003), so we will only keep datetime and 149343_00060_00003, and rename them to date and flow, respectively. We will also convert the units from “cubic feet per second” to “cubic meters per second”.
  
  ## keep only relevant columns
  dat_flow <- dat_flow[c("datetime",grep("[0-9]$",colnames(dat_flow),value=TRUE))]
  ## nicer column names
  colnames(dat_flow) <- c("date","flow")
  dat_flow$flow = as.numeric(dat_flow$flow)
  ## convert cubic feet to cubic meters
  dat_flow$flow <- dat_flow$flow/35.3147
  
  ## some days in 2017 were ice, so they are NAs. I change to 0
  
  dat_flow$flow[which(is.na(dat_flow$flow))] = 0
  ## flow by year & month
  
  # dat_flow[,"year"] <- as.integer(sub("^([0-9]{4})-([0-9]{2})-([0-9]{2})","\\1",
  #                                     dat_flow[,"date"]))
  # dat_flow[,"month"] <- as.integer(sub("^([0-9]{4})-([0-9]{2})-([0-9]{2})","\\2",
  #                                      dat_flow[,"date"]))
  # 
  # easier
  library(lubridate)
  dat_flow$year = year(dat_flow$date)
  dat_flow$month = month(dat_flow$date)
  
  dat_flow <- dat_flow[,c("year","month","flow")]
  
  
# Winter maximum
# 
# We are interested in the maximum of the daily peak flows from October through March during the first year that juveniles are rearing in streams. This means we need to combine flow values for the fall of year tt with those in the spring of year t+1t+1. Therefore, the flow time series will begin in yr_frst = 1978`; the last year of flow data will then beyr_last - age_min + n_fore + flow_lag = 2014`.
  
  ## autumn flows in year t
  flow_aut <- subset(dat_flow, (month>=10 & month<=12)
                     & year >= yr_frst & year <= yr_last-age_min+n_fore)
  ## spring flows in year t+1
  flow_spr <- subset(dat_flow, (month>=1 & month<=3)
                     & year >= yr_frst+flow_lag & year <= yr_last-age_min+n_fore+flow_lag)
  ## change spr year index to match aut
  flow_spr[,"year"] <- flow_spr[,"year"] - flow_lag
  ## combined flows indexed to brood year and calculate max flow over time period
  dat_flow_wtr <- aggregate(flow ~ year, data=rbind(flow_aut,flow_spr), max)
  ## change year index to brood year
  dat_flow_wtr[,"year"] <- dat_flow_wtr[,"year"] 
  ## for plotting purpose later
  colnames(dat_flow_wtr)[2] <- "Winter"
  
# Summer minimum
# 
# Now we will calculate the minimum flow juveniles would experience during their first summer (June through September).
  
## summer flows in year t
flow_sum <- subset(dat_flow, (month>=6 & month<=9)
                     & year >= yr_frst+flow_lag & year <= yr_last-age_min+n_fore+flow_lag)
## change year index to brood year
flow_sum[,"year"] <- flow_sum[,"year"] - flow_lag
## combined flows indexed to brood year and calculate max flow over time period
dat_flow_sum <- aggregate(flow ~ year, data=flow_sum, min)
## for plotting purpose later
colnames(dat_flow_sum)[2] <- "Summer" 


## URL for NPGO data
url_NPGO <- "http://www.o3d.org/npgo/npgo.php"
## raw NPGO data 
NPGO_raw <- readLines(url_NPGO)
## line with data headers
hdr_NPGO <- which(lapply(NPGO_raw,grep,pattern="YEAR")==1, arr.ind=TRUE)
## print PDO metadata
print(NPGO_raw[seq(hdr_NPGO)],quote=FALSE)

# Next, we will extract the actual NPGO indices for the years of interest and inspect the file contents. We also want the average NPGO annual index from January 1 through December 31 during the first year that the juvenile steelhead are in the ocean (i.e., during their second year of life). Therefore, we need NPGO values from yr_frst + marine_lag = 1980 through  yr_last - age_min + n_fore + marine_lag = 2015.

## NPGO data for years of interest
dat_NPGO<- read.table(url_NPGO, header=FALSE, stringsAsFactors=FALSE,
                      skip=hdr_NPGO + (yr_frst-1950)*12, nrows = n_yrs*12)
colnames(dat_NPGO) <- c("year","month","NPGO")

## select only years of interest indexed by brood year 
dat_NPGO <- dat_NPGO[dat_NPGO$year >= yr_frst+marine_lag &
                       dat_NPGO$year <= yr_last-age_min+n_fore+marine_lag,]
dat_NPGO <- aggregate(dat_NPGO$NPGO, by = list(year = dat_NPGO$year), FUN = mean)
dat_NPGO <- data.frame(year=seq(yr_frst,yr_last-age_min+n_fore), NPGO=dat_NPGO[,2])
colnames(dat_NPGO) <- c("year","NPGO")
dat_NPGO
  
# Spring Transition Index
# We calculated the spring transition index (STI) from the daily coastal upwelling index (CUI) provided by NOAA’s Pacific Fisheries Environmental Laboratory (PFEL); you can find more information here. We begin by downloading the raw CUI data and viewing the metadata.

## URL for CUI data
url_CUI <- "https://www.pfeg.noaa.gov/products/PFELData/upwell/daily/p06dayac.all"
## raw CUI data from PFEL
CUI_raw <- readLines(url_CUI)
## line with data headers
hdr_CUI <- which(lapply(CUI_raw,grep,pattern="YYYYMMDD")==1, arr.ind=TRUE)
## print CUI metadata
print(CUI_raw[seq(hdr_CUI-1)],quote=FALSE)


## get daily CUI data
dat_CUI <- read.table(url_CUI, header=TRUE, stringsAsFactors=FALSE, skip=hdr_CUI-1)
## extract year from date
dat_CUI$yr <- gsub("[0-9]{4}$","",dat_CUI$YYYYMMDD)
## select only years of interest
cui <- dat_CUI[dat_CUI$yr >= yr_frst+marine_lag & dat_CUI$yr <= yr_last-age_min+n_fore+marine_lag,]
## calculate cumulative upwelling by year
cum_CUI <- tapply(cui[,"Index"], cui$yr, cumsum)
## function to return spring transition index
get_STI <- function(x, day_max=200) {
  return(min(which(x==min(x[1:day_max]))))
}
## calc STI for each year
dat_STI <- data.frame(year=seq(yr_frst,yr_last-age_min+n_fore),STI=sapply(cum_CUI,get_STI))
 

# Hatchery releases
# The numbers of hatchery fish released each year is listed in a file on the project site. They have already been lagged 2 years (i.e., brood year + 2) to account for the potential competitive interactions during their juvenile life stage. (We will divide the release number by 1000 for plotting purposes.)

dat_hrel <- read.csv(text = getURL(paste0(ex_url,fn_hrel))) 
dat_hrel 


# Combine all covariates
# The last thing we will do is combine the covariates into one file and standardize them to have zero-mean and unit-variance.

## covariate(s)
dat_cvrs <- Reduce(function(...) merge(..., all=T),
                   list(dat_flow_wtr,dat_flow_sum,dat_NPGO,dat_STI,dat_hrel))
## drop year col
dat_cvrs <- dat_cvrs[,-1] 
## transform the covariates to z-scores
scl_cvrs <- scale(dat_cvrs) 
## total number of covariates
n_cov <- dim(scl_cvrs)[2] 

