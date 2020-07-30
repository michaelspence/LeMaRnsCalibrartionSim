sv_ll<- likelihood(obs_c,obs_s,params,model_run_act,sigma_c,sigma_s,fish_rate,time=20)

start<-Sys.time()
for(i in 1:100) {test <- pda_mcmc(obs_c,obs_s,params,model_run,sigma_c,sigma_s,Fs,sur_ef,llike=sum(like),sigmasFF,change=c(3,7))}
Sys.time()-start

start<-Sys.time()
for(i in 1:10) {test <- mda_mcmc(params,model_run,theta,Fs,sur_ef,prop_sdR,prop_sdF,sigma_c,sigma_s,obs_c,obs_s,like)}
Sys.time()-start

start<-Sys.time()
for(i in 1:100) {test <- calderhead_move(params,model_run,theta,Fs,sur_ef,cd_calder_move,sigma_c,sigma_s,obs_c,obs_s,like, M=10)}
Sys.time()-start

start<-Sys.time()
for(i in 1:100) {test <- gibbs_move(params,model_run,theta,Fs,sur_ef)}
Sys.time()-start


##### Now for the MCMC part
## tuning pars
prop_sdR <- c(0.2,0.2,0.2)
prop_sdF <- c(0.1,0.1,0.1)
sigmasFF <- 0.2
#load("C:/Users/MS23/OneDrive - CEFAS/BrexitMSY/CelticSeaPaper/FishFisheries/LeMaRns example/prelim.Rdata")
#cd_calder_move <- chol(2*cov(thetas))

load("~/Run.Rdata")
cd_calder_move <- chol(cov(thetas)/4)


library(abind)
theta <- c(26,11,7,1,0.1,0.3,0.6)
theta <- c(24, 10,6,0, 0.5,0.5,0.5)
thetas <- theta
Fs <- matrix(runif(20 *3 ,0,1.5),20,3)
Fss <- array(Fs,dim=c(20,3,1))
test <- gibbs_move(params,model_run,theta,Fs,sur_ef)
sigma_c <- test[[1]]; sigma_s <- test[[2]]; like <- test[[3]]
sigmas_s <- sigma_s  <- c(0.5,0.5,0.5)
sigmas_c <- sigma_c <- c(0.5,0.5,0.5)

tmp <- get_likelihood(theta,Fs,params,sur_ef,sigma_c,sigma_s,obs_c=obs_c,obs_s=obs_s)
like <-  tmp[[1]]; model_run <- tmp[[2]]

for(ii in 1:1000){
  test <- pda_mcmc_spin_up(obs_c,obs_s,theta,params,model_run,sigma_c,sigma_s,Fs,sur_ef,llike=sum(like),sigmasFF=sigmasFF,change=which(sample(c(0,1),20,replace=T,prob=c(0.8,0.2))==1))
  like <- test[[1]]; Fs <- test[[2]]
 # test <- gibbs_move(params,model_run,theta,Fs,sur_ef)
#  sigma_c <- test[[1]]; sigma_s <- test[[2]]; like <- test[[3]]
 # sigmas_s <- sigma_s 
#  sigmas_c <- sigma_c 
}

tmp <- get_likelihood(theta,Fs,params,sur_ef,sigma_c,sigma_s,obs_c=obs_c,obs_s=obs_s)
like <-  tmp[[1]]; model_run <- tmp[[2]]
likes <- sum(like)
ii <- 1
N <- 5e4

library(parallel)
cl <- makeCluster(10)
clusterExport(cl=cl, c("get_likelihood", "run_LeMans","likelihood","get_CPG"))

while(ii < N){
  ## mda-mcmc
  test <- mda_mcmc(params,model_run,theta,Fs,sur_ef,prop_sdR,prop_sdF,sigma_c,sigma_s,obs_c,obs_s,like)
  theta <- test[[1]]; model_run <- test[[2]]; like <- test[[3]]
  
  ## pda-mcmc
  for(kk in 1:10){
    test <- pda_mcmc(obs_c,obs_s,theta,params,model_run,sigma_c,sigma_s,Fs,sur_ef,llike=sum(like),sigmasFF=sigmasFF,change=which(sample(c(0,1),20,replace=T,prob=c(0.3,0.7))==1))
    like <- test[[1]]; Fs <- test[[2]]
  }
  ## calderhead
  test <- calderhead_move(params,model_run,theta,Fs,sur_ef,cd_calder_move,sigma_c,sigma_s,obs_c,obs_s,like, M=9)
  theta <- test[[1]]; model_run <- test[[2]]; like <- test[[3]]
  
  ## gibbs move
  test <- gibbs_move(params,model_run,theta,Fs,sur_ef)
  sigma_c <- test[[1]]; sigma_s <- test[[2]]; like <- test[[3]]
  
  ## collecting
  thetas <- rbind(thetas,theta)
  Fss <- abind(Fss,Fs)
  sigmas_c <- rbind(sigmas_c,sigma_c)
  sigmas_s <- rbind(sigmas_s,sigma_s)
  likes <- c(likes,sum(like))
  if ((ii / 1000) == floor(ii / 1000)){
    plot.ts(thetas)
  }
  ii <- ii + 1
  print(ii)
}

plot.ts(t(Fss[1:10,1,]))
plot.ts(t(Fss[1:10,2,]))
plot.ts(t(Fss[1:10,3,]))
plot.ts(thetas)
