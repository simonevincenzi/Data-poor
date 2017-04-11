if(!require("R2jags")) {
  install.packages("R2jags")
  library("R2jags")
}
if(!require("RCurl")) {
  install.packages("RCurl")
  library("RCurl")
}
if(!require("gsl")) {
  install.packages("gsl")
  library("gsl")
}

Re2prec <- function(x,map="round",prec=1) {
  ## 'map' can be round, floor, or ceiling
  ## 'prec' is nearest value (eg, 0.1 means to nearest tenth); default 1 gives normal behavior
  if(prec<=0) { stop("\"prec\" cannot be less than or equal to 0") }
  do.call(map,list(x/prec))*prec
}