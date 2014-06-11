# Generate observation from simulations

Since the measles data represent daily count of new cases (incidence), you need to account for the fact that probably not all cases are observed. You can do that by defining an observation process, that is a probabilistic model that can link your true - unobserved - state to your noisy - observed - data.


```r
SIR_generateObservation <- function(model.traj, theta){

    # compute daily incidence 
    
    # sample from observation process

    # return the same model.traj data.frame with one additional column called  `observation`
}
```

Check:


```r
SIR <- fitmodel(name="SIR",state.variables, list.fitparams, SIR_initialiseState, SIR_simulateDeterministic, SIR_generateObservation, verbose=TRUE)
```

# Load data in fitmodel and plot some fit
Load measles data:


```r
data(measles)
SIR <- fitmodel(name="SIR", state.variables, list.fitparams, SIR_initialiseState, SIR_simulateDeterministic, SIR_generateObservation, data=measles, verbose=TRUE)
```

Simulate and plot some fits:


```r
plotThetaFit(SIR$theta,SIR)
```

Since the observation process is stochastic you might want to plot several replicates of observed simulations to appreciate their variability:


```r
plotThetaFit(SIR$theta,SIR,n.replicates=100)
```

# Joint log-likelihood

Now you have an observation process you can formulate the log-likelihood of your `data` given our model parameter `theta` and the corresponding trajectory `model.traj`.


```r
SIR_logLikelihood <- function(data, theta, model.traj){

    # compute daily incidence

    # compute log-likelihood of data given theta and model.traj

    # return log-likelihood value

}
```

Update fitmodel.


```r
SIR <- fitmodel(name="SIR",state.variables, list.fitparams, SIR_initialiseState, SIR_simulateDeterministic, SIR_generateObservation, data=measles, log.likelihood=SIR_logLikelihood, verbose=TRUE)
```

# Marginal likelihood function

Since we are interested into the marginal likelihood of `theta`, you can now write a function that will only take `theta` as input, simulate the trajectory and compute the log-likelihood. This is your first generic function (type 2) as it takes a `fitmodel` as argument.


```r
my_marginalLogLikelihoodDeterministic <- function(theta, fitmodel) {

    # simulate model at successive observation times of data

    # compute log-likelihood

    # return log-likelihood
}
```

# Bonus
You can use optimization methods (`?optim`) to find the `theta` that maximizes the likelihood, i.e. using a frequentist approach.

# Navigate
Previous: [My first fitmodel](first_fitmodel.md) Next: [Run a MCMC](mcmc.md)
