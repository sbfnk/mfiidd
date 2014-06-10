# Generate observation from simulations


```r
my_SIR_generateObservation <- function(model.traj, theta){

    # daily incidence needed
    
    # sample from observation process

    # return the same data.frame updated
    return(model.traj)
}
```

Check:


```r
SIR <- fitmodel(name="SIR",state.variables,list.fitparams,my_SIR_initialiseState,my_SIR_simulateDeterministic,my_SIR_generateObservation,verbose=TRUE)
```

# Load data in fitmodel and plot some fit
Load measles data:


```r
data(measles)
SIR <- fitmodel(name="SIR",state.variables,list.fitparams,myy_SIR_initialiseState,my_SIR_simulateDeterministic,my_SIR_generateObservation,data=measles,verbose=TRUE)
```

Simulate and plot some fits:


```r
plotThetaFit(SIR$theta,SIR)
```

Since the observation process is stochastic you might want to simulate several observed simulation to appreciate their variability:


```r
plotThetaFit(SIR$theta,SIR,n.replicates=100)
```

# Joint log-likelihood


```r
my_SIR_logLikelihood <- function(data, theta, model.traj){

    # daily incidence needed

    # compute log likelihood of data given theta and model.traj

    # return log likelihood value

}
```

Update fitmodel.


```r
SIR <- fitmodel(name="SIR",state.variables,list.fitparams,myy_SIR_initialiseState,my_SIR_simulateDeterministic,my_SIR_generateObservation,data=data,log.likelihood=my_SIR_logLikelihood,verbose=TRUE)
```

# Marginal likelihood function

Generic function (type 2).


```r
my_marginalLogLikelihoodDeterministic <- function(theta, fitmodel) {

    # simulate model at successive observation times of data

    # compute log-likelihood

    # return log-likelihood
    return(log.likelihood)
}
```

# Bonus
Fit some parameters using a frequentist approach.

# Navigate
Next: [Run a MCMC](mcmc.md)
Previous: [My first fitmodel](first_fitmodel.md)
