library(LeMaRns)
### create a model with only 3 species
set.seed(14)
spec_par <- NS_par[c(1,11,20),]
tau <- NS_tau[c(1,11,20),c(1,11,20)]
### adjustment for re-parameterisation
spec_par$a <- exp(17.429 - 3.792 * log(spec_par$Linf))
spec_par$b <- exp(c(11,7.5,1))

toy_params <- LeMansParam(spec_par, tau=tau, eta=rep(0.25, 3), L50=spec_par$Lmat, other=exp(27))

# Define fishing effort
effort <- matrix(c(0.1,0.3,0.6), 100, dim(toy_params@Qs)[3],byrow = T)

# Run the model
model_run <- run_LeMans(toy_params, years=100, effort=effort)
plot_biomass(toy_params,model_run)

### need to add research vessel (from Walker)
custom_q <- matrix(0,32,3)
custom_q[,1] <- c(0.017254858,0.043039643,0.284392929,rep(0.430695848,29))
custom_q[,2] <- c(0.096113412,rep(0.24028353,31))
custom_q[,3] <- c(0.001761793,0.006166274,0.039181631,0.092776614,0.145237111,0.193340804,0.232928633,0.2564465,0.27551981,0.28422822,0.290723986,0.296690919,0.305572708,0.316083176,0.334233041,0.358307207,0.388331838,0.416365328,0.455782759,0.499352575,0.536897304,0.587032109,0.640318341,0.695650091,0.740381939,0.740381939,0.740381939,0.740381939,0.740381939,0.740381939,0.740381939,0.740381939
)
toy_params@Qs <- abind::abind(toy_params@Qs,custom_q)
dimnames(toy_params@Qs)[[3]][4] <- "survey"

### 20 year time series.
###

### then generate some data.

## This is now my data. Now there are 3 annual Fs (3 X 20), 3 uncertain parameters (2 X 3) + 1 global parameter
## 64 parameters!!!

## make it random walk
fish_rate <- abs(cbind(abs(cumsum(rnorm(20,0,0.075)) + 0.1),cumsum(rnorm(20,0,0.075)) + 0.6,cumsum(rnorm(20,0,0.075)) + 0.25,rnorm(20,0.0001,0.0001)))

model_run_act <- run_LeMans(toy_params, years=20, effort=fish_rate,N0=model_run@N[,,1001])
plot_biomass(toy_params,model_run_act)
catches <- get_CPG(toy_params,model_run_act,effort=fish_rate)
com_catches_t <- apply(catches[,1:3,],c(2,3),sum)

sur_catch <-  com_catches <- matrix(0,20,3)
for (i in 1:20){
  com_catches[i,] <- log(rowSums(com_catches_t[,10 * (i-1) + (1:10)])) - log(1e6)
  sur_catch[i,] <- log(rowSums(catches[,4,10 * (i-1) + (1:10)])) - log(1e6)
}

### add some noise
obs_c <- com_catches + matrix(rnorm(60,0,c(0.38,0.44,0.5)),20,3,byrow=T)
obs_s <- sur_catch + matrix(rnorm(60,0,c(0.6,0.7,0.6)),20,3,byrow=T)

sur_ef <- fish_rate[,4]
params <- toy_params
save(obs_c,obs_s,sur_ef,params,file="data.Rdata")
