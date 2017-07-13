test = matrix(0,nrow  = 6, ncol = 10)
L_inf = 3750; k_vb = 0.057; t0_vb = -0.21
t
pos.1 = 1
pos.2 = 2
pos.3 = 3
pos.4 = 4
pos.5 = 5

est[pos.1,] = runif(n = ncol (test),min = 0, max = 1)
test[pos.1,c(3,5,7)] = 0 
test[pos.2,] = sample(x = 1:5, size = ncol(test),replace = T)

test[pos.3,] = rnorm(n = ncol(test), mean = L_inf, sd = L_inf/3)
test[pos.4,] = rnorm(n = ncol(test), mean = k_vb, sd = k_vb/3)
test[pos.5,] = rnorm(n = ncol(test), mean = t0_vb, sd = -t0_vb/3)

test[6,] = ifelse(test[pos.1,] == 0, 0, L_inf * (1 - exp(-k_vb*(test[pos.5,]-t0_vb))))

calc = which(test[pos.1,]!=0)

test[6,calc] = round(sapply(calc, function(x) {
  test[6,x] = test[pos.3,x] * (1 - exp(-test[pos.4,x]*(test[pos.2,x]-test[pos.5,x])))
}))

vb_growth = function(test, pheno = postage = pos.1, linf = pos) {
  which(test[])
}