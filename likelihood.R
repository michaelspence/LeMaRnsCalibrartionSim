### The likelihood of the model after time runs

likelihood <- function(obs_c,obs_s,params,model_run,sigma_c,sigma_s,Fs,sur_ef,time=20){
  fish_rate <- cbind(Fs,sur_ef)
  colnames(fish_rate) <- NULL
  catches <- get_CPG(params,model_run,effort=fish_rate)
  com_catches_t <- apply(catches[,1:3,],c(2,3),sum)
  sur_catches <-  com_catches <- matrix(0,time,3)
  for (i in 1:time){
    com_catches[i,] <- log(rowSums(com_catches_t[,10 * (i-1) + (1:10)])) - log(1e6)
    sur_catches[i,] <- log(rowSums(catches[,4,10 * (i-1) + (1:10)])) - log(1e6)
  }
  return(rbind(dnorm(t(obs_c),t(com_catches),sigma_c,log=T),dnorm(t(obs_s),t(sur_catches),sigma_s,log=T)))
}


log_like <- function(model_run,i,obs_c,obs_s,sigma_c,sigma_s,fishing_rate){
  catches <- get_CPG(params,model_run,effort=fishing_rate)
  com_catches_t <- apply(catches[,1:3,],c(2,3),sum)
  com_catches <- log(rowSums(com_catches_t[,(1:10)])) - log(1e6)
  sur_catches <- log(rowSums(catches[,4,(1:10)])) - log(1e6)
  return(sum(dnorm(obs_c[i,],com_catches,sigma_c,log=T) + dnorm(obs_s[i,],sur_catches,sigma_s,log=T)))
}


### update the fishing levels

pda_mcmc <- function(obs_c,obs_s,theta,params,model_run,sigma_c,sigma_s,Fs,sur_ef,llike,sigmasFF,change=1:20){
  params@other <- exp(theta[1])
  ## recruits
  params@recruit_params[[1]][2] <- exp(theta[2]);params@recruit_params[[2]][2] <- exp(theta[3]);params@recruit_params[[3]][2] <- exp(theta[4]); 
  tmp <- theta[5:7]
  fish_rate <- cbind(Fs,sur_ef)
  colnames(fish_rate) <- NULL
  Q_0_ <-Q_0 <- M_0 <- model_run
  fish_rate_ <- fish_rate
  alpha <- rep(0)
  alpha1 <- rep(0)
  llike_ <- 0
  hast <- 0
    for (i in 1:20){
      if (any(change==i)){
        tmp <- rnorm(3,log(fish_rate[i,1:3]),sigmasFF)
        ## accept or reject
        if (all(exp(tmp) < 2)){
          M_0_ <- run_LeMans(params, years=1, effort=t(c(exp(tmp),fish_rate[i,4])),N0=M_0@N[,,dim(M_0@N)[3]])
          ll_p <- log_like(M_0_,i,obs_c,obs_s,sigma_c,sigma_s,t(c(exp(tmp),fish_rate[i,4])))
          Q_0_ <- run_LeMans(params, years=1, effort=t(c(exp(tmp),fish_rate[i,4])),N0=Q_0@N[,,dim(Q_0@N)[3]])
          alpha1b <-  log_like(Q_0_,i,obs_c,obs_s,sigma_c,sigma_s,t(c(exp(tmp),fish_rate[i,4])))
        }
        else{
          ll_p <- -Inf
          alpha1 <- Inf
          alpha1b <- -Inf
        }
        M_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=M_0@N[,,dim(M_0@N)[3]])
        Q_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=Q_0@N[,,dim(Q_0@N)[3]])
        ll_r <- log_like(M_0,i,obs_c,obs_s,sigma_c,sigma_s,t(fish_rate[i,]))
        alpha1a <- log_like(Q_0,i,obs_c,obs_s,sigma_c,sigma_s,t(fish_rate[i,]))
        alpha1 <- exp(alpha1a - alpha1b)
        alpha <- min(1,exp(ll_p - ll_r))
        
        if( runif(1) <= alpha){
          llike_ <- llike_ + ll_p
          M_0 <- M_0_
          fish_rate_[i,] <- c(exp(tmp),fish_rate[i,4])
          #print(alpha1)
          hast <- hast + log(min(alpha1,1)) - log(alpha)
        }
        else{
          llike_ <- llike_ + ll_r
          fish_rate_[i,] <- fish_rate[i,]
          #print(1/alpha1)
          hast <- hast + log(1 - min(1/alpha1,1)) - log(1 - alpha)
        }
      }else{
        M_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=M_0@N[,,dim(M_0@N)[3]])
        Q_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=Q_0@N[,,dim(Q_0@N)[3]])
        llike_ <- llike_ +  log_like(M_0,i,obs_c,obs_s,sigma_c,sigma_s,t(fish_rate[i,]))
      }
    }
  #browser()
  ### now do it for them all at once
  if(log(runif(1)) < (llike_ - llike + hast)){
    print("accept fish_rate")
    return(list(llike=llike_,fish_rate=fish_rate_[,1:3]))
  }
  else{
    return(list(llike=llike,fish_rate=fish_rate[,1:3]))
  }
}