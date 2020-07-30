### history matching
library(LeMaRns)
load("data.Rdata")

calc_catch <- function(params,model_run,effort){
  catches <- get_CPG(params,model_run,effort=effort)
  com_catches_t <- apply(catches[,1:3,],c(2,3),sum)
  times <- (dim(model_run@N)[3] - 1)*params@phi_min
  sur_catch <-  com_catches <- matrix(0,times,3)
  for (i in 1:times){
    com_catches[i,] <- log(rowSums(com_catches_t[,10 * (i-1) + (1:10)])) - log(1e6)
    sur_catch[i,] <- log(rowSums(catches[,4,10 * (i-1) + (1:10)])) - log(1e6)
  }
  return(list(com=com_catches,sur=sur_catch))
}

######################## just for history matching
hist_mat <- function(input,params){
  params@other <- exp(input[1])
  ## recruits
  params@recruit_params[[1]][2] <- exp(input[2]);params@recruit_params[[2]][2] <- exp(input[3]);params@recruit_params[[3]][2] <- exp(input[4]); 
  model_run <- run_LeMans(params, years=100, effort=t(matrix(c(input[5:7],0),4,100)))
  model_run1 <- run_LeMans(params,years=1,effort=t(matrix(c(input[8:10],sur_ef[1]),4,1)),N0=model_run@N[,,dim(model_run@N)[3]])
  effort <- t(matrix(c(input[8:10],sur_ef[1]),4,1))
  #browser()
  return(unlist(calc_catch(params,model_run1,effort)))
}


library(randtoolbox) ## for sobol indicies
sam1 <- sobol(5000,3*3 + 1)

## other first - between 20 and 30
sam1[,1] <- sam1[,1] *20 +10
## Rmax between 0 and 20
sam1[,2:4] <- sam1[,2:4] * 20
## spin-up F betweem 0 and 1.5
sam1[,5:10] <- sam1[,5:10] * 1.5

#hist_mat(sam1[1,],params=params)


library(pbapply) ## do in parallel
round1 <- pbapply(sam1,1,hist_mat,params=params)
c_r1 <- round1 - c(obs_c[1,],obs_s[1,])
###
par(mfrow=c(2,3))
plot(sam1[,1],c_r1[1,],ylim=c(-10,10),xlab="other",ylab="sprat difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,1],c_r1[2,],ylim=c(-10,10),xlab="other",ylab="mackerel difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,1],c_r1[3,],ylim=c(-10,10),xlab="other",ylab="cod difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,1],c_r1[4,],ylim=c(-10,10),xlab="other",ylab="sprat difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,1],c_r1[5,],ylim=c(-10,10),xlab="other",ylab="mackerel difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,1],c_r1[6,],ylim=c(-10,10),xlab="other",ylab="cod difference",main="Survey catches")
abline(h=0,col="red")

###
par(mfrow=c(2,3))
plot(sam1[,2],c_r1[1,],ylim=c(-10,10),xlab="sprat b",ylab="sprat difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,3],c_r1[2,],ylim=c(-10,10),xlab="mackerel b",ylab="mackerel difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,4],c_r1[3,],ylim=c(-10,10),xlab="cod b",ylab="cod difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,2],c_r1[4,],ylim=c(-10,10),xlab="sprat b",ylab="sprat difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,3],c_r1[5,],ylim=c(-10,10),xlab="mackerel b",ylab="mackerel difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,4],c_r1[6,],ylim=c(-10,10),xlab="cod b",ylab="cod difference",main="Survey catches")
abline(h=0,col="red")

par(mfrow=c(2,3))
plot(sam1[,5],c_r1[1,],ylim=c(-10,10),xlab="Spin up F sprat",ylab="sprat difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,6],c_r1[2,],ylim=c(-10,10),xlab="Spin up F mackerel",ylab="mackerel difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,7],c_r1[3,],ylim=c(-10,10),xlab="Spin up F cod",ylab="cod difference",main="Commercial catches")
abline(h=0,col="red")
plot(sam1[,5],c_r1[4,],ylim=c(-10,10),xlab="Spin up F sprat",ylab="sprat difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,6],c_r1[5,],ylim=c(-10,10),xlab="Spin up F mackerel",ylab="mackerel difference",main="Survey catches")
abline(h=0,col="red")
plot(sam1[,7],c_r1[6,],ylim=c(-10,10),xlab="Spin up F cod",ylab="cod difference",main="Survey catches")
abline(h=0,col="red")

par(mfrow=c(2,3))
plot(sam1[,8],c_r1[1,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam1[,9],c_r1[2,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam1[,10],c_r1[3,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam1[,8],c_r1[4,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam1[,9],c_r1[5,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam1[,10],c_r1[6,],ylim=c(-10,10))
abline(h=0,col="red")

### seen enough 

sam2 <- sobol(5000,3*3 + 1,init = F)

## other first - between 20 and 30
sam2[,1] <- sam2[,1] *6 +24
## Rmax between 0 and 20
sam2[,2:4] <- t(t(sam2[,2:4]) * c(13,14,20) + c(7,6,0))
## spin-up F betweem 0 and 1.5
sam2[,5:7] <- t(t(sam2[,5:7]) * c(1.5,1.5,1) + c(0,0,0.0))
sam2[,8:10] <- sam2[,8:10] * 1.5

round2 <- pbapply(sam2,1,hist_mat,params=params)
c_r2 <- round2 - c(obs_c[1,],obs_s[1,])
###
par(mfrow=c(2,3))
plot(sam2[,1],c_r2[1,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,1],c_r2[2,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,1],c_r2[3,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,1],c_r2[4,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,1],c_r2[5,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,1],c_r2[6,],ylim=c(-10,10))
abline(h=0,col="red")

###
par(mfrow=c(2,3))
plot(sam2[,2],c_r2[1,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,3],c_r2[2,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,4],c_r2[3,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,2],c_r2[4,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,3],c_r2[5,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,4],c_r2[6,],ylim=c(-10,10))
abline(h=0,col="red")

par(mfrow=c(2,3))
plot(sam2[,5],c_r2[1,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,6],c_r2[2,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,7],c_r2[3,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,5],c_r2[4,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,6],c_r2[5,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,7],c_r2[6,],ylim=c(-10,10))
abline(h=0,col="red")

par(mfrow=c(2,3))
plot(sam2[,8],c_r2[1,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,9],c_r2[2,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,10],c_r2[3,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,8],c_r2[4,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,9],c_r2[5,],ylim=c(-10,10))
abline(h=0,col="red")
plot(sam2[,10],c_r2[6,],ylim=c(-10,10))
abline(h=0,col="red")

save.image("after_2_rounds.Rdata")
## done quite well
## ready for MCMC now