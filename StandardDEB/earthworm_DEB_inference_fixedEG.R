#Bayesian inference for the scaled standard DEB model
#Supplementary code for the manuscript Boersch-Supan and Johnson 2018: Two case studies detailing Bayesian parameter inference for dynamic energy budget models

library(devtools)
#install the development version of deBInfer that is able to recalculate boundary values in the MCMC
install_github("pboesu/debinfer", ref = "recalc-inits")
#install DEB utility functions
install_github("pboesu/DEButilities")

#load additional required packages
library(deSolve)
library(deBInfer)
library(truncdist)


#data
tW <- readr::read_delim("data/Lumbricus_terrestris_mass_Butt1993.txt", delim = " ")
tN <- readr::read_delim("data/Lumbricus_terrestris_production_Butt1993.txt", delim = " ")

#construct parameter vector using AmP entry for Lumbricus terrestris
#http://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/Lumbricus_terrestris/Lumbricus_terrestris_par.html
params_Lumter <- c(L_m = 1.406,
                   p_Am = 898.555,
                   v = 0.1426,
                   k_J = 0.002,
                   kap = 0.985,
                   T_A = 5000,
                   T_ref = 293.15,
                   T_b = 293.15,
                   E_G = 4150,
                   f = 1,
                   E_Hb = 1.428,
                   E_Hp = 728.6)

#chemical parameters
w_E = 23.9 # molecular weight of reserve g mol^-1
d_v = d_E = 0.16 # specific density of structure
mu_E = 550000 # chemical potential of reserve J / mol
mu_v = 500000 # chemical potential of structuure J / mol
kap_R = 0.95 # fixed canonical value


#a function to calculate selected compound parameters
compound_params <- function(params){
  L_m = params["L_m"]
  E_G = params["E_G"]
  v = params["v"]
  kap = params["kap"]
  k_J = params["k_J"]
  E_Hp = params["E_Hp"]
  p_Am =params["p_Am"]
  f = params["f"]
  a_b = params["a_b"]
  
  E_m = p_Am/ v
  g = E_G/ kap/ E_m;
  k_M = v / (g * L_m)
  l_T = 0
  k = k_J / k_M
  V_m =(v/(k_M*g))^3
  
  w_E = 23.9 # molecular weight of reserve g mol^-1
  d_v = 0.16 # specific density of structure
  mu_E = 550000 # chemical potential of reserve J / mol 
  
  w = p_Am * w_E / (v *d_v * mu_E) #omega
  
  u_Hp = E_Hp/(g*E_m*V_m)
  
  return(list(E_m = unname(E_m), g = unname(g), k_M = unname(k_M), l_T = unname(l_T), k = unname(k), V_m = unname(V_m), u_Hp = unname(u_Hp), w = unname(w)))
}

#calculate compound parameters
parscomp_Lumter <- compound_params(params_Lumter)

#load R DEB model function
#Using the R implementation is not recommended because it slows down inference considerably, but the R model is included for reference
#source("src/scaled_std_deb.R")

#load C model
#this requires a C compiler
#details on how to specify and compile C model code is available in the corresponding deSolve and deBInfer vignettes
#they can be displayed with the commands
# vignette("deBInfer_compiled_code", package = "deBInfer")
# vignette("compiledCode", package = "deSolve")
system("R CMD SHLIB src/scaled_std_deb.c")
# on a windows machine the compiled object will have the file extension .dll rather than .so
dyn.load("src/scaled_std_deb.so") 


#calculate scaled maturity at birth and puberty 
v_Hb <- DEButilities::get_v_Hb(params_Lumter[["E_Hb"]], params_Lumter[["kap"]], parscomp_Lumter[["g"]], params_Lumter[["p_Am"]], params_Lumter[["v"]], params_Lumter[["L_m"]])
v_Hp <- DEButilities::get_v_Hb(params_Lumter[["E_Hp"]], params_Lumter[["kap"]], parscomp_Lumter[["g"]], params_Lumter[["p_Am"]], params_Lumter[["v"]], params_Lumter[["L_m"]])

u_Hb = params_Lumter[["E_Hb"]] / params_Lumter[["p_Am"]] * parscomp_Lumter[["g"]]^2 * parscomp_Lumter[["k_M"]]^3 / params_Lumter[["v"]]^2
u_Hp = params_Lumter[["E_Hp"]] / params_Lumter[["p_Am"]] * parscomp_Lumter[["g"]]^2 * parscomp_Lumter[["k_M"]]^3 / params_Lumter[["v"]]^2

#calculate initial scaled reserve and length at birth
uE0_lb = DEButilities::get_ue0(p = c(parscomp_Lumter[["g"]], parscomp_Lumter[["k"]], v_Hb) )
uE0 = uE0_lb["uE0"]
l_b = uE0_lb["lb"]
L_b = l_b * params_Lumter[["L_m"]]


#assemble initial value vector to evolve DEB model from birth
states_Lumbricus_terrestris <- c(e = 1, l = unname(l_b), uH = unname(u_Hb), uR = 0)
#set time intervals at which to solve the ODE
times_Lumbricus_terrestris <- seq(0,400, by = 1)

#solve the ODE at the specified times
out <- ode(y = states_Lumbricus_terrestris, times = times_Lumbricus_terrestris, 
           func = "d_scaled_std_deb", dllname = "scaled_std_deb", initfunc = "init_scaled_std", nout = 1, outnames = "EWw",
           parms = params_Lumter, method = 'lsoda')
#plot trajectories of the state variable
plot(out)




#conversion to observable quantities
#simulated egg mass
w_E = 23.9 # molecular weight of reserve g mol^-1
d_v = d_E = 0.16 # specific density of structure
mu_E = 550000 # chemical potential of reserve J / mol
#egg weight
Ww_0_sim = uE0 * params_Lumter[["v"]]^2 / (parscomp_Lumter[["g"]]^2 * parscomp_Lumter[["k_M"]]^3) * params_Lumter[["p_Am"]] * w_E/ mu_E/ d_E; # g, egg wet weight

#simulated age and mass at birth and puberty
tau_bp_sim = DEButilities::get_tp(p  = c(parscomp_Lumter[["g"]], parscomp_Lumter[["k"]], parscomp_Lumter[["l_T"]], v_Hb, v_Hp))
tau_b_sim = tau_bp_sim["tb"]
tau_p_sim = tau_bp_sim["tp"]
age_b_sim = tau_b_sim / parscomp_Lumter[["k_M"]]

Ww_b_sim = L_b^3 *(1 + params_Lumter["f"] * parscomp_Lumter[["w"]])

#simulated length and mass at puberty
l_p = tau_bp_sim["lp"]
Ww_p_sim = (l_p*params_Lumter["L_m"])^3 *(1 + params_Lumter["f"] * parscomp_Lumter[["w"]])

#simulated age at puberty
tau_p_sim = DEButilities::get_tp(p  = c(parscomp_Lumter[["g"]], parscomp_Lumter[["k"]], parscomp_Lumter[["l_T"]], v_Hb, v_Hp))["tp"]

age_p_since_birth_sim = tau_p_sim / parscomp_Lumter[["k_M"]] - age_b_sim


#figure of observable quantities for paper
pdf("earthworm_dat1.pdf", height = 4)
par(mfrow = c(1,2))
#calculate time-weight trajectory - add-my-pet approach predict using von Bert growth rate
ir_B = 3/ parscomp_Lumter[["k_M"]] + 3 * params_Lumter[["f"]] * params_Lumter[["L_m"]]/ params_Lumter[["v"]];
rT_B = 1/ ir_B; # d, 1/von Bert growth rate
L = params_Lumter[["L_m"]] - (params_Lumter[["L_m"]] - L_b) * exp( - rT_B * times_Lumbricus_terrestris);           # cm, structural length at time
EWw = L^3 * (1 + params_Lumter[["f"]] * parscomp_Lumter[["w"]]);                                 # g, wet weight


plot(times_Lumbricus_terrestris, EWw, col='red', type = 'l', xlab = "Time since birth (days)", ylab = "Wet weight (g)")
points(tW$days_since_hatching, tW$mass_g)

#calculate time-reproduction trajectory
kap_R = 0.95
#E_0 = 205.884	
R = kap_R * out[,"uR"] / uE0

plot(times_Lumbricus_terrestris, R, col='red', type = 'l', xlab = "Time since birth (days)", ylab = "Cumulative number of cocoons")
points(tN$days_since_puberty + 91.5, tN$cumulative_eggs)
dev.off()
#abline(v = age_p_since_birth_sim, lty  =2)


#sample observations from simulation
Lumter_simulated_obs <- data.frame(time = times_Lumbricus_terrestris, EWw = EWw, R = R)[seq(1, length(times_Lumbricus_terrestris), by = 30),]

#add simulated noise
sdlog_for_Ww = 0.05
sdlog_for_R = 0.1
#set RNG seed to ensure consistency of simulated observations
set.seed(20180118)
Lumter_simulated_obs$EWw_noisy <- rlnorm(nrow(Lumter_simulated_obs), meanlog = log(Lumter_simulated_obs$EWw), sdlog = sdlog_for_Ww)
Lumter_simulated_obs$R_noisy <- rlnorm(nrow(Lumter_simulated_obs), meanlog = log(Lumter_simulated_obs$R), sdlog = sdlog_for_R)

#plot simulated noisy "observations"
plot(EWw ~ time, data = Lumter_simulated_obs, type = "l", ylim = c(0, max(Lumter_simulated_obs$EWw_noisy)))
points(EWw_noisy ~ time, data = Lumter_simulated_obs)
plot(R ~ time, data = Lumter_simulated_obs, type = "l", ylim = c(0, max(Lumter_simulated_obs$R_noisy)))
points(R_noisy ~ time, data = Lumter_simulated_obs)


# observation model for inference
Lumter_obs_model_WN<-function(data, sim.data, samp){
  ec<- 1e-4 #numerical correction, follows Johnson et al. 2013
  
  #wet weight (not using the von-Bert approach)
  w_E = 23.9 # molecular weight of reserve g mol^-1
  d_v = d_E = 0.16 # specific density of structure
  mu_E = 550000 # chemical potential of reserve J / mol
  mu_v = 500000 # chemical potential of structuure J / mol
  kap_R = 0.95 # fixed canonical value
  
  #recalculate compound parameters
  parscomp_Lumter <- compound_params(samp)
  
  #calculate maturities etc. at transitions
  v_Hb = DEButilities::get_v_Hb(E_Hb = samp[["E_Hb"]], kap = samp[["kap"]], g = parscomp_Lumter[["g"]], p_Am = samp[["p_Am"]], v = samp[["v"]], L_m = samp[["L_m"]])
  v_Hp = DEButilities::get_v_Hb(E_Hb = samp[["E_Hp"]], kap = samp[["kap"]], g = parscomp_Lumter[["g"]], p_Am = samp[["p_Am"]], v = samp[["v"]], L_m = samp[["L_m"]])
  uE0_lb_info = DEButilities::get_ue0(p = c(parscomp_Lumter[["g"]], parscomp_Lumter[["k"]], v_Hb) )
  uE0 = uE0_lb_info["uE0"]
  l_b = uE0_lb_info["lb"]
  
  #calculate omega
  omega = unname(samp['p_Am'] * w_E / (samp['v'] * d_v * mu_E))
  
  #calculate weight and reproduction trajectories at current parameter values
  w.temp <- (sim.data[,"l"]*samp["L_m"])^3 * (1 + samp["f"] * omega)
  r.temp <- kap_R * sim.data[,"uR"] / uE0
  
  #likelihood on mass observations
  llik.Ww<-sum(dlnorm(data$EWw_noisy+ec, meanlog=log(w.temp+ec), sdlog=samp["sdlog.EWw"], log=TRUE)) #observation uncertainty is hardcoded to the data estimate
  #likelihood on cocoon observations
  #subset to post-puberty times to remove zero observations which aren't well captured by the lognormal probability model
  r.temp <- r.temp[sim.data[, "time"] > 100]
  r.data <- data$R_noisy[sim.data[, "time"] > 100]
  llik.R<-sum(dlnorm(r.data, meanlog=log(r.temp+ec), sdlog=samp["sdlog.R"], log=TRUE))
  
  
  #constraint on birth and puberty maturity
  llik.Hbp <- log(as.numeric(samp[["E_Hb"]] < samp[["E_Hp"]]))
  
  #constraint on l_b (reach_birth.m)
  llik.lb <- log(((l_b <= samp[["f"]])*(parscomp_Lumter[["k"]]*v_Hb) <= (samp[["f"]]/(parscomp_Lumter[["g"]] + samp[["f"]]) * l_b^2 * (parscomp_Lumter[["g"]] + l_b) )))
  
  #constraint on reaching puberty
  llik.puberty <- log((parscomp_Lumter[["k"]]*v_Hp) <= (samp[["f"]]*(samp[["f"]] - parscomp_Lumter[["l_T"]])^2))
  
  #likelihood on egg mass
  Ww_0 = uE0 * samp[["v"]]^2 / (parscomp_Lumter[["g"]]^2 * parscomp_Lumter[["k_M"]]^3) * samp[["p_Am"]] * w_E/ mu_E/ d_E; # g, egg wet weight
  llik.Ww0 <- dtrunc(Ww_0, "norm", mean = Ww_0_sim, sd = 0.01*Ww_0_sim, a = 0, log = TRUE) #mean from Butt 1993, sd guestimated from Garcia & Fragoso 2002 and Lofs-Holmin 1982
  
  
  #age at birth and puberty
  tau_bp = DEButilities::get_tp(p  = c(parscomp_Lumter[["g"]], parscomp_Lumter[["k"]], parscomp_Lumter[["l_T"]], v_Hb, v_Hp))
  tau_b = tau_bp["tb"]
  tau_p = tau_bp["tp"]
  age_b = tau_b / parscomp_Lumter[["k_M"]]
  age_p_since_birth = tau_p / parscomp_Lumter[["k_M"]] - age_b
  
  Ww_b = (l_b*samp["L_m"])^3 *(1 + samp["f"] * parscomp_Lumter[["w"]])
  
  llik.Ww_b = dtrunc(Ww_b, "norm", mean = Ww_b_sim, sd = 0.01*Ww_b_sim, a = 0, log = TRUE)
  
  #length at puberty
  l_p = tau_bp["lp"]
  Ww_p = (l_p*samp["L_m"])^3 *(1 + samp["f"] * parscomp_Lumter[["w"]])
  
  llik.Ww_p = dtrunc(Ww_p, "norm", mean = Ww_p_sim, sd = 0.01*Ww_p_sim, a = 0, log = TRUE)
  
  #likelihood on hatching age
  llik.age_b <- dtrunc(age_b, "norm", mean = age_b_sim, sd = 0.01*age_b_sim, a = 0, log = TRUE) #mean from Butt 1993, sd guestimated
  
  #likelihood on puberty age
  llik.age_p <- dtrunc(age_p_since_birth, "norm", mean = age_p_since_birth_sim, sd = 0.01*age_p_since_birth_sim, a = 0, log = TRUE) #mean from Butt 1993, sd guestimated
  
  #constraint on E_G
  llik.EG <- log(as.numeric(samp["E_G"] > d_v*mu_v/w_E))
  
  #total likelihood
  llik<-llik.Ww + llik.R + llik.Hbp + llik.lb + llik.puberty + llik.Ww0 + llik.Ww_p +  llik.Ww_p + llik.age_b + llik.age_p + llik.EG
  
  return(llik)
}


#now set up parameters for inference procedure
kap <- debinfer_par(name = "kap", var.type = "de", fixed = FALSE,
                    value = 0.8, prior="beta", hypers=list(shape1 = 2, shape2 = 2),
                    prop.var=5e-7, samp.type="rw")

L_m <- debinfer_par(name = "L_m", var.type = "de", fixed = FALSE,
                    value = 1, prior="trunc", hypers=list(spec = "norm", a = 0, mean = 1, sd = 1),
                    prop.var=0.0001, samp.type="rw")

p_Am<- debinfer_par(name = "p_Am", var.type = "de", fixed = FALSE,
                    value = 100, prior="norm",  hypers=list(mean = 900, sd = 300),
                    prop.var=100, samp.type="rw")

v <- debinfer_par(name = "v", var.type = "de", fixed = FALSE,
                  value = 0.2, prior="trunc",  hypers=list(spec = "norm", a = 0, mean = 0.2,  sd = 0.2), #pseudodata approach
                  prop.var=2e-5, samp.type="rw")

k_J <- debinfer_par(name = "k_J", var.type = "de", fixed = TRUE,
                    value = 0.002)

T_A <- debinfer_par(name = "T_A", var.type = "de", fixed = TRUE,
                    value = 5000)

T_ref <- debinfer_par(name = "T_ref", var.type = "de", fixed = TRUE,
                      value = 293)

T_b <- debinfer_par(name = "T_b", var.type = "de", fixed = TRUE,
                    value = 293)

f <- debinfer_par(name = "f", var.type = "de", fixed = TRUE,
                  value = 1)

E_G <- debinfer_par(name = "E_G", var.type = "de", fixed = TRUE,
                    value = 4150, #the following prior and sampler specification is ignored when fixed = TRUE. When fixed = FALSE, an informative prior based on AmP data for all species is used 
                    prior="trunc", hypers=list(spec = "norm", a = d_v*mu_v/w_E, mean = 4200,  sd = 2000),
                    prop.var=200, samp.type="rw")# back of envelope calculation from add-my-pet median of 1e4.4 for E_G/d_v

E_Hb <- debinfer_par(name = "E_Hb", var.type = "de", fixed = FALSE,
                     value = 1,  prior="trunc", hypers=list(spec = "norm", a = 0, mean = 0,  sd = 100),
                     prop.var=0.02, samp.type="rw")

E_Hp <- debinfer_par(name = "E_Hp", var.type = "de", fixed = FALSE,
                     value = 1000,  prior="norm", hypers=list(mean = 1000,  sd = 1000),
                     prop.var=100, samp.type="rw")


#observation parameters for uni-variate data
sdlog.EWw <- debinfer_par(name = "sdlog.EWw", var.type = "obs", fixed = FALSE,
                          value = sdlog_for_Ww, prior="trunc", hypers=list(spec = "norm", a = 0, mean=0,sd=0.1),
                          prop.var=5e-4, samp.type="rw-ref")

sdlog.R <- debinfer_par(name = "sdlog.R", var.type = "obs", fixed = FALSE,
                        value = sdlog_for_R, prior="trunc", hypers=list(spec = "norm", a = 0, mean=0,sd=0.1),
                        prop.var=1e-3, samp.type="rw-ref")

#function that recalculate initial values l_b and uHb 
# ----deinitfunc-------------------------------------------------
get_lb_uHb <- function(inits, params){
  #recalculate compound parameters
  parscomp <- compound_params(params)
  v_Hb <- DEButilities::get_v_Hb(params["E_Hb"], params["kap"], parscomp[["g"]], params["p_Am"], params["v"], params["L_m"])
  u_Hb <- v_Hb * (1 - params["kap"])
  lb <- DEButilities::get_lb(c(parscomp[["g"]], parscomp[["k"]], v_Hb))["lb"]
  inits["l"] <- unname(lb)
  inits["uH"] <- unname(u_Hb)
  return(inits)
}

#specification of initial values (at birth)
# ----inits---------------------------------------------------------------
e_init <- debinfer_par(name = "e", var.type = "init", fixed = TRUE, value = 1)

#calculate length at birth  
l_init <- debinfer_par(name = "l", var.type = "initfunc", fixed = TRUE, deinitfunc = get_lb_uHb, value = NA)

uH_init <- debinfer_par(name = "uH", var.type = "init", fixed = TRUE, value = unname(u_Hb))

uR_init <- debinfer_par(name = "uR", var.type = "init", fixed = TRUE, value = 0)

#collate inference parameters
mcmc.pars <- setup_debinfer(L_m, p_Am, v, k_J, kap, T_A, T_ref, T_b, E_G, f, E_Hb, E_Hp, sdlog.EWw, sdlog.R, e_init, l_init, uH_init, uR_init)

#collate data (i.e. simulated observations)
Lumter_data <- data.frame(time = Lumter_simulated_obs$time, EWw = Lumter_simulated_obs$EWw_noisy, R = Lumter_simulated_obs$R_noisy)



#run the MCMC inference


iter = 1000 #(for the figures in the paper inference was conducted using three chains (i.e. three runs of this entire script) with 150,000 iterations each. This took about eight hours per chain to compute on a machine with Intel(R) Xeon(R) CPUs (Model E5-2603 v3 @ 1.60GHz).
#1000 iterations should take 2-3 minutes to compute on a modern desktop computer.
#at 1000 iterations the chain will likely be far from good parameter estimates. In our experience, this model requires a burn-in of 5000-10000 iterations before the chain arrives in a realistic region of the parameter space


#reset RNG seed by querying system time (so the MCMC chain is not identical every time we run the script)
set.seed(Sys.time())


#if the inference procedure fails on the start values it often helps to debug the observation model to understand what causes the initial likelihood to be infinite 
#debug(Lumter_obs_model_WN)

#run inference. see ?de_mcmc for documentation or refer to the deBInfer package vignettes for more details
Lumter_samples <- de_mcmc(N= iter, # run MCMC for this many iterations
                          data=Lumter_simulated_obs, #observations
                          de.model="d_scaled_std_deb", #ODE model function; see # vignette("deBInfer_compiled_code", package = "deBInfer") for details on how to call a C model for inference
                          dllname = "scaled_std_deb", initfunc = "init_scaled_std", nout = 1, outnames = "EWw", #further parameters to the C code
                          obs.model=Lumter_obs_model_WN, #specify observation model function
                          all.params=mcmc.pars, # specify collated MCMC parameters
                          Tmax = max(Lumter_simulated_obs$time), #maximum timestep to solve the ODE for
                          data.times=c(Lumter_simulated_obs$time), #time points for whih there are data
                          cnt=200, #interval at which to update 
                          plot=FALSE, #plot traces as you run the inference. This slows down the inference, but can be informative when tuning samplers
                          sizestep=0.1, 
                          solver="ode", # solver wrapper choice, here deSolve::ode
                          verbose = TRUE, # verbose output from solver
                          verbose.mcmc = TRUE, #verbose updates from the MCMC sampler
                          maxsteps = 5000, # argument to solver
                          method = "lsoda" #ODE solver method choice
                          )

#some basic plots of the inference output
#MCMC trace
plot(Lumter_samples, density = FALSE, smooth = FALSE)
#posterior densities
plot(Lumter_samples, trace = FALSE, smooth = FALSE)

#simulate posterior trajectories of state variable
post_traj <- post_sim(Lumter_samples, n=100, times=1:500, burnin = 100, output = 'all', prob = 0.95, maxsteps = 15000, dllname = "scaled_std_deb",
                      initfunc = "init_scaled_std", nout = 1, outnames = "EWw")
#plot posterior trajectories
plot(post_traj, plot.type = "ensemble", col="#D60000")

#compare the simulated "true" trajectory with the posterior trajectories
#this requires the plyr package, which can be installed from CRAN using
#install.packages("plyr")
par(mfrow = c(2,2))
plot(out, which = "e", mfrow = NULL, lwd = 3)
plyr::l_ply(post_traj$sim, function(x) lines(x[,c("time", "e")], lty = 2, col = "red"))
plot(out, which = "l", mfrow = NULL, lwd = 3)
plyr::l_ply(post_traj$sim, function(x) lines(x[,c("time", "l")], lty = 2, col = "red"))
plot(out, which = "uH", mfrow = NULL, lwd = 3)
plyr::l_ply(post_traj$sim, function(x) lines(x[,c("time", "uH")], lty = 2, col = "red"))
plot(out, which = "uR", mfrow = NULL, lwd = 3)
plyr::l_ply(post_traj$sim, function(x) lines(x[,c("time", "uR")], lty = 2, col = "red"))