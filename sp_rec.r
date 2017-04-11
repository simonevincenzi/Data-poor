# Here is the relationship between spawner and subsequent recruits (a), assuming mean values for all covariates. Gray lines show the median relationship for each of the 38 years based on atat. Note that for plotting purposes only in (b) and (c), the density in the largest bin for each parameter contains counts for all values greater or equal to that. Vertical arrows under the x-axes in (b) and (c) indicate the 2.5th, 50th, and 97.5th percentiles.

layout(matrix(c(1,1,2,3),2,2),c(3,2),c(1,1))
CI_vec <- c(0.025,0.5,0.975)
offSet <- 0.06

## posterior of spawners
sDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,quantile,CI_vec)
sDat <- sDat[,1:(n_yrs-age_min+n_fore)]
## posterior of recruits
rDat <- exp(apply(mod_fit$BUGSoutput$sims.list$tot_ln_Rec,2,quantile,CI_vec))
## median values for a & b
aa <- apply(mod_fit$BUGSoutput$sims.list$ln_Rkr_a,2,median)
bb <- apply(mod_fit$BUGSoutput$sims.list$Rkr_b,2,median)

## empty plot space for spawner-recruit relationships
dd <- 3000
yM <- Re2prec(max(rDat),"ceiling",dd)
yM <- 30000
xM <- Re2prec(max(sDat),"ceiling",dd)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0,0,0))
plot(sDat[2,],rDat[2,], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3", type="n",
     xaxs="i", yaxs="i", ylab="Recruits (1000s)", xlab="Spawners (1000s)", cex.lab=1.2,
     xaxt="n", yaxt="n")
axis(1, at=seq(0,xM,dd*2), labels=seq(0,xM,dd*2)/1000)
axis(2, at=seq(0,yM,dd*2), labels=seq(0,yM,dd*2)/1000)
for(i in 1:length(aa)) { lines(seq(xM)*exp(aa[i]-bb*seq(xM)), col="darkgray") }
## add S-R estimates and medians
abline(a=0,b=1,lty="dashed")
nCB <- n_yrs-age_max
points(sDat[2,1:nCB],rDat[2,1:nCB], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3")
segments(sDat[2,1:nCB],rDat[1,1:nCB],sDat[2,1:nCB],rDat[3,1:nCB], col="blue3")
segments(sDat[1,1:nCB],rDat[2,1:nCB],sDat[3,1:nCB],rDat[2,1:nCB], col="blue3")
nTB <- dim(sDat)[2]
clr <- rgb(100, 0, 200, alpha=seq(200,100,length.out=age_max-age_min+n_fore), maxColorValue=255)
segments(sDat[2,(nCB+1):nTB],rDat[1,(nCB+1):nTB],sDat[2,(nCB+1):nTB],rDat[3,(nCB+1):nTB], col=clr)
segments(sDat[1,(nCB+1):nTB],rDat[2,(nCB+1):nTB],sDat[3,(nCB+1):nTB],rDat[2,(nCB+1):nTB], col=clr)
points(sDat[2,(nCB+1):nTB],rDat[2,(nCB+1):nTB],
       xlim=c(0,xM), ylim=c(0,yM), pch=16, col=clr)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(a)")

## posterior for exp(a)
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
par(mai=c(0.8,0.4,0.3,0.1))
## Ricker a
R_alpha_est <- mod_fit$BUGSoutput$sims.list$alpha
alphaCI <- quantile(R_alpha_est,c(0.025,0.5,0.975))
R_alpha_est[R_alpha_est>9] <- 9
hist(R_alpha_est,freq=FALSE,xlab="",main="",breaks=seq(0,9,0.2),
     col=clr, border="blue3", ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/12
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(paste("Instrinsic productivity ",(e^italic(a)))), 1, line=3, cex=1)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(b)")

## posterior for a/b
par(mai=c(0.8,0.4,0.3,0.1))
aa <- matrix(mod_fit$BUGSoutput$sims.array[,,"E_Rkr_a"],ncol=1)
bb <- matrix(mod_fit$BUGSoutput$sims.array[,,"Rkr_b"],ncol=1)
R_b_est <- aa/bb
R_b_est <- R_b_est[R_b_est > 0]
R_b_CI <- quantile(R_b_est,c(0.025,0.5,0.975))
R_b_est[R_b_est>2e4] <- 2e4
brks <- seq(Re2prec(min(R_b_est),"floor",2000),2e4,length.out=length(seq(0,9,0.2)))
hist(R_b_est, freq=FALSE, breaks=brks, col=clr, border="blue3",
     xlab="", xaxt="n", yaxt="n",
     main="", ylab="", cex.lab=1.2)
axis(1, at=seq(Re2prec(min(R_b_est),"floor",2000),1.8e4,4000))
aHt <- (par()$usr[4]-par()$usr[3])/12
arrows(R_b_CI,par()$usr[3],R_b_CI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(paste("Carrying capacity (",italic(a)/italic(b),")")), 1, line=3, cex=1)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(c)")