get_likelihood<- function(theta,Fs,params,sur_ef,sigma_c,sigma_s,obs_c=obs_c,obs_s=obs_s){
  params@other <- exp(theta[1])
  ## recruits
  params@recruit_params[[1]][2] <- exp(theta[2]);params@recruit_params[[2]][2] <- exp(theta[3]);params@recruit_params[[3]][2] <- exp(theta[4]); 
  tmp <- theta[5:7]
  model_run <- run_LeMans(params, years=100, effort=t(matrix(c(tmp,0),4,100)))
  fish_rate <- as.matrix(cbind(Fs,sur_ef))
  colnames(fish_rate) <- NULL
  model_run_act <- run_LeMans(params, years=20, effort=fish_rate,N0=model_run@N[,,dim(model_run@N)[3]])
  #browser()
  like<- likelihood(obs_c,obs_s,params,model_run_act,sigma_c,sigma_s,Fs,sur_ef)
  return(list(like=like,model_run=model_run))
}

### for the MDA move

mda_mcmc <- function(params,model_run,theta,Fs,sur_ef,prop_sdR,prop_sdF,sigma_c,sigma_s,obs_c,obs_s,like){
  prop_ <- parSapply(cl,1:3,mda_mcmc1,params=params,theta=theta,Fs=Fs,sur_ef=sur_ef,prop_sdR=prop_sdR,prop_sdF=prop_sdF,sigma_c=sigma_c,sigma_s=sigma_s,obs_c=obs_c,obs_s=obs_s,like=like) ## could be in parallel
  #prop_ <- sapply(1:3,mda_mcmc1,params=params,theta=theta,Fs=Fs,sur_ef=sur_ef,prop_sdR=prop_sdR,prop_sdF=prop_sdF,sigma_c=sigma_c,sigma_s=sigma_s,obs_c=obs_c,obs_s=obs_s,like=like) ## could be in parallel
  thetar <- matrix(rep(c(theta[1],prop_[2,],prop_[3,]),4),4,byrow=T)
  for (i in 2:4){
    thetar[i,c(i,i + 3)] <- c(prop_[c(4,5),i-1]) ### always rejected bit
  }
  like_ <- parRapply(cl,thetar,get_likelihood,Fs=Fs,params=params,sur_ef=sur_ef,sigma_c=sigma_c,sigma_s=sigma_s,obs_c=obs_c,obs_s=obs_s)
  tmp <- lapply(like_[-1],function(x,y){rowSums(x[[1]]-y)},y=like_[[1]][[1]])
  tmp1 <- matrix(unlist(tmp),nrow=3,byrow = T)
  tmp <- c(tmp1[1,1] + tmp1[1,4],tmp1[2,2] + tmp1[2,5],tmp1[3,3] + tmp1[3,6])
  tmp1 <- ifelse(tmp > 0,1,exp(tmp))
  acc_rej<- sum(like_[[1]][[1]]) - sum(like) - sum(log(ifelse(prop_[6,]==1,prop_[1,],1-prop_[1,]))) + sum(log(ifelse(prop_[6,]==1,tmp1,1-tmp1)))
  if (log(runif(1)) < acc_rej){
    print("accept")
    return(list(thetar[1,],like_[[1]][[2]],like_[[1]][[1]]))
  }
  return(list(theta,model_run,like))
}

mda_mcmc1 <- function(i,params,theta,Fs,sur_ef,prop_sdR,prop_sdF,sigma_c,sigma_s,obs_c,obs_s,like){
  theta_ <- rnorm(1,theta[i+1],prop_sdR[i])
  spin_upF_ <- abs(rnorm(1,theta[i+4],prop_sdF[i]))
  #browser()
  params@recruit_params[[1]][2] <- exp(theta[2]);params@recruit_params[[2]][2] <- exp(theta[3]);params@recruit_params[[3]][2] <- exp(theta[4]); 
  params@other <- exp(theta[1])
  params@recruit_params[[i]][2] <- exp(theta_)
  #params@recruit_params[[i]][2] <- exp(theta[i+1])
  tmp <- theta[5:7]
  tmp[i] <- spin_upF_
  model_run <- run_LeMans(params, years=100, effort=t(matrix(c(tmp,0),4,100)))
  fish_rate <- as.matrix(cbind(Fs,sur_ef))
  colnames(fish_rate) <- NULL
  model_run_act <- run_LeMans(params, years=20, effort=fish_rate,N0=model_run@N[,,dim(model_run@N)[3]])
  #browser()
  like_ <- sum(likelihood(obs_c,obs_s,params,model_run_act,sigma_c,sigma_s,Fs,sur_ef)[c(i,i+3),])
  alpha <- exp(min(0,like_  - sum(like[c(i,i+3),])))
  if (runif(1) <= alpha){
    return(c(alpha,theta_,spin_upF_,theta[i+1],theta[i+4],1))
  }
  return(c(alpha,theta[i+1],theta[i+4],theta_,spin_upF_,0))
}

### calderhead move
calderhead_move <- function(params,model_run,theta,Fs,sur_ef,cd_calder_move,sigma_c,sigma_s,obs_c,obs_s,like, M=10){
  theta1 <- as.numeric(theta + t(cd_calder_move) %*% rnorm(7)   )
  theta_ <- abs(cbind(theta,theta1 + t(cd_calder_move) %*% matrix(rnorm(7 * M),nrow=7,ncol=M) ))
  theta_[1,]<- ifelse(theta_[1,] > 30, 60 - theta_[1,],theta_[1,])
  like_ <- parCapply(cl,theta_,get_likelihood,Fs=Fs,params=params,sur_ef=sur_ef,sigma_c=sigma_c,sigma_s=sigma_s,obs_c=obs_c,obs_s=obs_s)
  tmp <- unlist(lapply(like_,function(x){sum(x[[1]])}))
  number <- sample(M+1,1,prob=exp(tmp - max(tmp)))
  print(number)
  return(list(as.numeric(theta_[,number]),like_[[number]][[2]],like_[[number]][[1]]))
}

### gibbs move
gibbs_move <- function(params,model_run,theta,Fs,sur_ef){
  params@other <- exp(theta[1])
  ## recruits
  params@recruit_params[[1]][2] <- exp(theta[2]);params@recruit_params[[2]][2] <- exp(theta[3]);params@recruit_params[[3]][2] <- exp(theta[4]); 
  tmp <- theta[5:7]
  fish_rate <- as.matrix(cbind(Fs,sur_ef))
  colnames(fish_rate) <- NULL
  model_run_act <- run_LeMans(params, years=20, effort=fish_rate,N0=model_run@N[,,dim(model_run@N)[3]])
  catches <- get_CPG(params,model_run_act,effort=fish_rate)
  com_catches_t <- apply(catches[,1:3,],c(2,3),sum)
  sur_catches <-  com_catches <- matrix(0,20,3)
  for (i in 1:20){
    com_catches[i,] <- log(rowSums(com_catches_t[,10 * (i-1) + (1:10)])) - log(1e6)
    sur_catches[i,] <- log(rowSums(catches[,4,10 * (i-1) + (1:10)])) - log(1e6)
  }
  sigma_c <- sqrt(1/rgamma(3,shape=10 + 0.1,rate=0.1+
    rowSums((t(obs_c)-t(com_catches))^2)/2))
  sigma_s <- sqrt(1/ rgamma(3,shape=10 +0.1,rate=0.1 +
    rowSums((t(obs_s)-t(sur_catches))^2)/2))
  return(list(sigma_c,sigma_s,rbind(dnorm(t(obs_c),t(com_catches),sigma_c,log=T),dnorm(t(obs_s),t(sur_catches),sigma_s,log=T))))
}