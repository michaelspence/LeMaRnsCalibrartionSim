pda_mcmc_spin_up <- function(obs_c,obs_s,theta,params,model_run,sigma_c,sigma_s,Fs,sur_ef,llike,sigmasFF,change=1:20){
    params@other <- exp(theta[1])
    ## recruits
    params@recruit_params[[1]][2] <- exp(theta[2]);params@recruit_params[[2]][2] <- exp(theta[3]);params@recruit_params[[3]][2] <- exp(theta[4]); 
    tmp <- theta[5:7]
    fish_rate <- cbind(Fs,sur_ef)
    colnames(fish_rate) <- NULL
    Q_0_ <-Q_0 <- M_0 <- model_run
    fish_rate_ <- fish_rate
    alpha <- rep(0,20)
    alpha1 <- rep(0,20)
    llike_ <- 0
    hast <- 0
    for (i in 1:20){
      if (any(change==i)){
        tmp <- rnorm(3,log(fish_rate[i,1:3]),sigmasFF)
        if (all(exp(tmp) < 2)){
          M_0_ <- run_LeMans(params, years=1, effort=t(c(exp(tmp),fish_rate[i,4])),N0=M_0@N[,,dim(M_0@N)[3]])
          ll_p <- log_like(M_0_,i,obs_c,obs_s,sigma_c,sigma_s,t(c(exp(tmp),fish_rate_[i,4])))
        }
        else{
          ll_p <- -Inf
        }
        M_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=M_0@N[,,dim(M_0@N)[3]])
        ## accept or reject
        ll_r <- log_like(M_0,i,obs_c,obs_s,sigma_c,sigma_s,t(fish_rate[i,]))
        alpha <- min(1,exp(ll_p - ll_r))
        if( runif(1) <= alpha){
          llike_ <- llike_ + ll_p
          M_0 <- M_0_
          fish_rate_[i,] <- c(exp(tmp),fish_rate[i,4])
          #print(alpha1)
        }
        else{
          llike_ <- llike_ + ll_r
          fish_rate_[i,] <- fish_rate[i,]
          #print(1/alpha1)
        }
      }else{
        M_0 <- run_LeMans(params, years=1, effort=t(fish_rate[i,]),N0=M_0@N[,,dim(M_0@N)[3]])
        llike_ <- llike_ +  log_like(M_0,i,obs_c,obs_s,sigma_c,sigma_s,t(fish_rate[i,]))
        #print(llike_)
      }
    }
    ### now do it for them all at once
    if(log(runif(1)) < (llike_ - llike)){
      print(paste("accept fish_rate",change))
      return(list(llike=llike_,fish_rate=fish_rate_[,1:3]))
    }
    else{
      return(list(llike=llike,fish_rate=fish_rate[,1:3]))
    }
  }