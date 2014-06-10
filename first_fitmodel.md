
# Objectives
In this section, you'll build your first `fitmodel`: an SIR model. By the end of the day you'll fit this model to a measles outbreak dataset:

![plot of chunk measles-data](knitr/figure/measles-data.png) 

Before doing any coding you should think how to model this outbreak, what assumptions you want to make and make a draw of your model showing the notation for parameters and state variables. Then you can build the corresponding `fitmodel` step-by-step.

# Parameters
Define all parameters as `fitparam` objects (only the `name` and `value` arguments are required for the moment) and group them as a list:


```r
R0 <- fitparam(name="R0",value=10)
infectious.period <- fitparam(name="IP",value=5)
list.fitparams <- list(R0,infectious.period)
```

This is of course an incomplete example, and you probably need more parameters to describe your model (think about initial conditions for instance). As a rule of thumb, a parameter can be anything you might want to estimate in your model. Of course feel free to choose different parameter `name`s, usually a short acronym is convenient because your are going to reuse it a lot later ;)

# State variables and initial conditions
Based on your choice, create a vector `state.variables` with the state variables names. To simulate your model you'll need some initial values for you state variables. However, since most of these initial conditions are unknown and will have to be estimated you'll need a function to set them, based on the parameter values `theta <- getParameterValues(list.fitparams)`:


```r
SIR_initialiseState <- function(theta) {
# set initial conditions as a function of the parameter vector theta

# return a named vector of state variables values
}
```

At this stage you can start to check your work by building your first `fitmodel`:


```r
SIR <- fitmodel(name="SIR",state.variables,list.fitparams,SIR_initialiseState,verbose=TRUE)
```

Your `SIR` object should contain 4 non `NULL` elements: `name`, `theta`, `state.variables` and `initialise.model`. Note: you can use the option `verbose=FALSE` to hide the details on the checks performed.

# Simulate the model

It might be useful to first write down the system of ordinary differential equations (ODEs). Now, write a function to simulate your model deterministically, you will use the function `ode()` of the `R` package `deSolve`.


```r
SIR_simulateDeterministic <- function(theta,state.init,times) {
 # The following function compute the derivative of the ODE system
 SIR_ode <- function(time,state,theta) {
 
 # type ?ode for guidance on what to put here
 
 }
 
 # simulate and return a data.frame
 trajectory <- data.frame(ode(y=state.init,times=times,func=SIR_ode,parms=theta))
 
 return(trajectory)
 }
```

Check your code by updating you `fitmodel`:


```r
SIR <- fitmodel(name="SIR",state.variables, list.fitparams, SIR_initialiseState, SIR_simulateDeterministic, verbose=TRUE)
```

# Plot some simulations
You should now be able to simulate your `fitmodel`:


```r
traj <- SIR$simulate.model(SIR$theta, SIR$initialise.state(SIR$theta), times=0:10)
```

and plot the output using the helpful `plotModelTraj()` function:


```r
plotModelTraj(SIR,traj)
```

Type `?plotModelTraj` for more potting options.

# Navigate
Previous: [Overview of fitcourseR](fitcourseR.md) Next: [Data and likelihood](data_likelihood.md)



