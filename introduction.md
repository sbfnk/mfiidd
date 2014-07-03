

# Introduction to model fitting in R

## Objectives
The aim of this first session is to set you up with a framework for model fitting. To do this, we will introduce you to the `fitR` package, which we have created to facilitate interaction during this course. This is not a full-fledged model fitting suite such as POMP. Instead, we use it to provide you with model code that follows a common, data and solutions to the practical questions.

In this session, you will

1. familiarise yourselves with the structure of models in the `fitR` package
2. combine the prior and likelihood to calculate the posterior for a simple SIR model
3. explore the posterior for a model that has a single parameter

## A simple SIR model

First, we will work with a simple SIR model. We will later fit this to a simple data set. The model has two parameters, the basic reproductive number $R_0$ and the infectious period $D$. The model equations are

$$
\begin{aligned}
\frac{dS}{dt} &= - \beta S \frac{I}{N}\\
\frac{dI}{dt} &= \beta S \frac{I}{N} - \nu I\\
\frac{dR}{dt} &= \nu I\\
\end{aligned}
$$

where $\beta=R_0/D$, $\nu = 1/D$ and $N = S + I + R$ is constant. To load this model into your **R** session, type


```r
example(SIR)
```

This will load a `fitmodel` object called `SIR` into your **R** session. A `fitmodel` is simply a collection of functions and variables that define a model. To see the objects that the `SIR` model contains, type


```r
names(SIR)
```

```
## [1] "name"         "state.names"  "theta.names"  "simulate"    
## [5] "genObsPoint"  "logPrior"     "pointLogLike"
```

A fitmodel object provides its description, as well as the names of its model variables and parameters. These can be accessed with


```r
SIR$name
```

```
## [1] "SIR with constant population size"
```

```r
SIR$state.names
```

```
## [1] "S" "I" "R"
```

```r
SIR$theta.names
```

```
## [1] "R0" "D"
```

Moreover, each `fitmodel` contains four functions: `simulate` to run the model, `logPrior` to calculate the logarithm of the prior, `pointLogLike` to calculate the log-likelihood of a given data point with respect to the model, and `generateObs` to generate observations from a model run. We will now look at these one after the other. If at any time you would like more information on the components of a `fitmodel`, you can do so by typing


```r
?fitmodel
```

### Simulate

To simulate a model run, we use the `simulate` contained within a `fitmodel`. To run the SIR model, we have to provide parameter values, initial conditions and the times at which we want model output.


```r
parameters <- c(R0 = 3, D = 1.5)
state.init <- c(S = 999, I = 1, R = 0)
times <- 1:100
traj <- SIR$simulate(parameters, state.init, times)
```

We can now look at the output of the model run


```r
head(traj)
```

```
##   time     S       I       R
## 1    1 999.0   1.000   0.000
## 2    2 994.8   3.774   1.393
## 3    3 979.4  13.997   6.606
## 4    4 925.8  48.804  25.348
## 5    5 775.1 140.292  84.572
## 6    6 514.9 264.171 220.923
```

To visualise a model, you can use `plotTraj`


```r
plotTraj(traj)
```

![](figure/unnamed-chunk-8.png) 

Try running the model with different values of the parameters and initial conditions, and output times.

Remember that if you want more detail on what the `simulate` and `plotTraj` functions do, you can always look at their implementation by typing their function name. To look at the `simulate` function of SIR, type


```r
SIR$simulate
```

### Prior

To evaluate the (logarithm of the) prior for a certain combination of parameters, use the `logPrior` function


```r
SIR$logPrior(parameters)
```

```
## [1] -7.996
```

Have a look at the `logPrior` function by typing


```r
SIR$logPrior
```

You will see that this calculates the prior of the parameters using uniform distributions on $R_0$ (between 1 and 100) and $D$ (between 0 and 30).

### Likelihood

The `pointLogLike` function is used to evaluate the likelihood of a data point


```r
SIR$pointLogLike(data.point = c(obs = 13), model.point = c(I = 24), parameters)
```

```
## [1] -5.237
```

Have a look at the `pointLogLike` function by typing


```r
SIR$pointLogLike
```

You will see that this calculates the likelihood of a data point taking its "obs" member (an observation) and evaluating it with respect to a Poisson distribution centred around the "I" member of the model point. In other words, it assumes that the observations follows a Poisson distribution centred around the current prevalence.

Let's load a test data set using


```r
data(epi)
```

This is contains several epidemic data sets. The first one, called `epi1`, has been created using this same SIR model, with an infectious period of 2, in a population of 1000, with 1 initial infected and the remaining population susceptible. Later, we will try to estimate the value of $R_0$. For now, you can look at this data set using


```r
epi1
```

```
##    time obs
## 1     1   3
## 2     2   1
## 3     3   1
## 4     4   2
## 5     5   5
## 6     6   8
## 7     7   9
## 8     8  18
## 9     9  26
## 10   10  27
## 11   11  38
## 12   12  60
## 13   13  81
## 14   14  90
## 15   15 111
## 16   16 122
## 17   17 123
## 18   18 120
## 19   19  92
## 20   20 105
## 21   21  84
## 22   22  65
## 23   23  66
## 24   24  41
## 25   25  52
## 26   26  43
## 27   27  29
## 28   28  15
## 29   29  12
## 30   30  14
## 31   31   7
## 32   32   5
## 33   33   9
## 34   34   2
## 35   35   4
## 36   36   1
## 37   37   1
## 38   38   0
```

The observations ("obs") are based on observed prevalence ("I"). You can plot the data using.


```r
plotTraj(epi1)
```

![](figure/unnamed-chunk-16.png) 

To calculate the log-likelihood of a set of parameters and initial conditions, we can use `trajLogLike`. This takes a model, a parameter vector, an initial state and a data set and simulates the model using the `simulate` function of the model with the given parameters and initial state. It then evaluates the log-likelihood at every data point using the `pointLogLike` function of the model. The result returned is the likelihood $p(\mathrm{Data}|X_0, \theta)$ of the chosen set of parameters and initial conditions.


```r
trajLogLike(SIR, parameters, state.init, epi1)
```

```
## [1] -6792
```

### Generate observations

The function `genObsPoint` generates single observation point from a model point. In that sense, it can be seen as the inverse of `pointLogLike`. Whereas `loglikePoint` evaluates the log-likelihood at a data point with respect to the model, `genObsPoint` takes a (randomly sampled) data point from the model.


```r
SIR$genObsPoint(model.point = c(I = 24), parameters)
```

```
## [1] 18
```

Of course, you might see a different number, as the result of this command is the outcome of a random draw.

To generate a whole trajectory, we can use `genObsTraj`. This uses the `simulate` function to simulate the model, and then the `genObsPoint` function at every time point to generate a trajectory of observations. It returns the trajectory, with an added column for the observations.


```r
obs.traj <- genObsTraj(SIR, parameters, state.init, epi1$time)
head(obs.traj)
```

```
##   time     S       I       R obs
## 1    1 999.0   1.000   0.000   1
## 2    2 994.8   3.774   1.393   6
## 3    3 979.4  13.997   6.606  15
## 4    4 925.8  48.804  25.348  43
## 5    5 775.1 140.292  84.572 146
## 6    6 514.9 264.171 220.923 254
```

If you run this multiple times, you will find that the outcome is different every time. This is because the observations are outcomes of random draws from the (deterministic) model trajectory. If we changed `genObsPoint` to be deterministic (instead of a random draw from a Poisson distribution), the outcome of `genObsTraj` would be the same every time.

## Calculate the posterior

We are now ready to calculate the value of the posterior $p(\theta, X_0|\mathrm{Data})$ (up to normalisation) for a given set of parameters and initial conditions, with respect to a data set. Let us write a function that does that. Below you find the skeleton of such a function. We have inserted comments for every line that you should insert. If you are struggling with at any point, click on the link below the code for some hints.


```r
# This is a function that takes 4 arguments:
# - fitmodel, a fitmodel object that defines the model dynamics,
#   prior and likelihoods.
# - theta, a named vector of parameters
# - state.init,  a named vector of initial conditions
# - data, the data set we are fitting the model to
# It should return the posterior for the given model, parameters,
# initial conditions and data.
my_posterior <- function(fitmodel, theta, state.init, data) {

    # calculate the fitmodel prior for parameter theta

    # calculate the fitmodel likelihood for parameter theta and
    # initial conditions state.init, with respect to the data set data.

    # return the logged posterior probability
}
```



If you have trouble filling any of the empty bits in, have a look at our [more guided example](posterior_example.md).

Check that your function returns a sensible value.


```r
my_posterior(SIR, parameters, state.init, epi1)
```

```
## [1] -6800
```

## Explore the posterior

Now we are ready to do parameter estimation by exploring the posterior at different values of the parameter. In the next practical we will see how to automate this step using MCMC, but for now let us simply evaluate the posterior at different values of the single unknown parameter, $R_0$. As stated above, the infectious period of the `epi1` data set was 2. You can evaluate the posterior at different values of $R_0$ using the `my_posterior` function you wrote above (or the one provided by clicking through to our [solution](posterior_example_solution.html)). In which range of $R_0$ is the posterior maximised?

You can also change the prior and likelihood definitions of the `SIR` model. Remember that you can see the definition of `SIR$logPrior` and `SIR$pointLogLike` by typing


```r
SIR$logPrior
SIR$pointLogLike
```

To change the prior, copy and paste the functions, change them, and reassign to their variables. For example, to change the point log-likelihood to follow a normal distribution, you could type


```r
SIR$pointLogLike <- function(data.point, model.point, theta) {
        # the prevalence is observed through a normally distributed process
        return(dnorm(x = data.point[["obs"]], mean = model.point[["I"]], log = TRUE))
}
```

Try different distributions for the prior and likelihood distributions. Do they change the shape of the posterior?

## Going further

If you knew that on average only 10% of cases were reported, how would you change the point likelihood? You can test this with a second data set, `epi2`, which you can have a look at with


```r
epi2
```

This was created with infectious period 2 and a reporting rate of 10%. Can you estimate $R_0$?

<div>
# Navigate
Previous: [Installation](README.md) Next: [MCMC](mcmc.md)
</div>
