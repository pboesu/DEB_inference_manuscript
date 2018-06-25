#Likelihood function for the snail DEBKiss model

DEBKiss_obs_model <- function(data, sim.data, samp){

    w.obs <- simdat <- NA
    test <- FALSE

    ## choose at which times to evaluate the likelihood
    w.obs <- which(sim.data[,"time"] %in% data$time)
   
    ## The simulated data set at only the correct tims
    simdat <- sim.data[w.obs,]

    ## this is a small correction used to keep the arguments of the
    ## log normal distribution from being exactly zero
    ec <- 1e-6

    ## The calculate of the cummulative number of eggs, calculated
    ## from the total amount of energy in the reproductive buffer,
    ## times the efficiency of the conversion of the reproductive
    ## buffer to eggs, divided by the energy per egg
    cumEggs <- (simdat[,"WR"]*as.numeric(samp[["yBA"]])/as.numeric(samp[["WB0"]])) + ec

    ## this is code to produce a plot of the data + simulation for
    ## testing purposes
    if(test){
        par(mfrow=c(1,2), bty="n")
        plot(simdat[,"time"], simdat[,"Lw"], type="l", ylab="structural length",
             xlab="time")
        points(data$time, data$L, col=2)
        plot(simdat[,"time"], cumEggs, type="l", ylab="total eggs", xlab="time")
        points(data$time, data$Egg, col=3)
    }

    ## log likelihood for the length portion of the data 
    llik.L <- sum(dlnorm(data$L, meanlog = log(simdat[,"Lw"] + ec),
                         sdlog = samp[["sdlog.L"]], log = TRUE))

    ## log likelihood for the cumulative eggs
    llik.E <- sum(dlnorm(data$Egg+ec, meanlog = log(cumEggs),
                         sdlog = samp[["sdlog.E"]], log = TRUE))

    ## sum the egg and length components together
    llik<- llik.L + llik.E

    if(test) print(c(llik.L, llik.E, llik)) ## check the likelihood

    return(llik)
}


