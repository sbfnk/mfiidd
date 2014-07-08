% Modelling interlude: Tristan da Cunha outbreak & the SEITL model.





# Objectives

This session will be focus on modelling rather than fitting (no worry, we will come back to fitting on the next session). The aim of this session is to familiarise yourselves with a data set: the Tristan da Cunha outbreak and a model called SEITL. We felt this session was necessary as we will use this dataset and this model as a toy example to illustrate all the concepts and methods from now until the end of the short-course. 

As you will read below, the SEITL model has been proposed as a mechanistic explanation for a two-wave influenza A/H3N2 epidemic that occurred on the remote island of Tristan da Cunha in 1971. Given the small population size of this island (284 inhabitants), random effects at the individual level may have had important consequences at the population level. This phenomenon is called *demographic stochasticity* and it could have played a significant role in the dynamics of the epidemic on Tristan da Cunha. 

Unsurpisingly, to account for demographic stochasticity we need a stochastic model. However, as you will see in this session and the following, simulating and fitting for a stochastic model is much more costly than for a deterministic model. This is why it is important to understand when and why you need to use stochastic or deterministic models. Actually, we hope that by the end of the day you'll be convinced that both approaches are in fact complementary.

In addition, because this session is devoted to modelling rather than fitting, we sought it was useful to show you how you can make your model biologically more realistic by using other distributions than the exponential distribution for the time spent in each compartment. 

In brief, in this session you will:

1. Familiarise yourself with the structure of the SEITL model and try to guess its parameter values from the litterature.
2. Compare deterministic and stochastic simulations in order to explore the role of demographic stochasticity in the dynamics of the SEITL model.
3. Use another distribution for the time spent in a compartment and assess its effect on the shape of the epidemic.

But first of all, let's have a look at the data.

# Tristan da Cunha outbreak

![](figure/import/TdC.png)

[Tristan da Cunha](http://en.wikipedia.org/wiki/Tristan_da_Cunha) is a volcanic island in the South Atlantic Ocean. It has been inhabited since the $19^{th}$ century and in 1971, the 284 islanders were living in the single village of the island: Edinburgh of the Seven Seas. Whereas the internal contacts were typical of close-knit village communities, contacts with the outside world were infrequent and mostly due to fishing vessels that occasionally took passengers to or from the island. These ships were often the cause of introduction of new diseases into the population. As for influenza, no epidemic had been reported since an epidemic of A/H1N1 in 1954. In this context of a small population with limited immunity against influenza, an unusual epidemic occurred in 1971, 3 years after the global emergence of the new subtype A/H3N2.

On August 13, a ship returning from Cape Town landed five islanders on Tristan da Cunha. Three of them developed acute respiratory disease during the 8-day voyage and the other two presented similar symptoms immediately after landing. Various family gatherings welcomed their disembarkation and in the ensuing days an epidemic started to spread rapidly throughout the whole island population. After three weeks of propagation, while the epidemic was declining, some islanders developed second attacks and a second peak of cases was recorded. The epidemic faded out after this second wave and lasted a total of 59 days.

Among the 284 islanders, 273 (96%) experienced at least one attack and 92 (32%) experienced two attacks, which is remarkable for influenza. Unfortunately, only 312 of the 365 attacks (85%) are known to within a single day of accuracy and constitute the dataset reported by [Mantle \& Tyrrell in 1973](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2130434/).

The dataset of daily incidence can be loaded and plotted as follows:


```r
data(FluTdC1971)
head(FluTdC1971)
##         date time obs
## 1 1971-08-13    1   0
## 2 1971-08-14    2   1
## 3 1971-08-15    3   0
## 4 1971-08-16    4  10
## 5 1971-08-17    5   6
## 6 1971-08-18    6  32

# plot daily observed incidence
plotTraj(data = FluTdC1971)
```

<img src="figure/SEITL_model/data.png" title="" alt="" style="display: block; margin: auto;" />

# SEITL model

One possible explanation for the rapid reinfections reported during this two-wave outbreak is that following recovery from a first infection, some islanders did not develop a long-term protective immunity and remained fully susceptible to reinfection by the same influenza strain that was still circulating. This can be modelled as follows:

![The SEITL model extends the classical SEIR model to account for the dynamics and host heterogeneity of the immune response among the islanders. Following recovery, hosts remain temporarily protected against reinfection thanks to the cellular response. Accordingly, they enter the T stage (temporary protection). Then, following down-regulation of the cellular response, the humoral response has a probability $\alpha$ to reach a level sufficient to protect against reinfection. In this case, recovered hosts enter the L stage (long-term protection) but otherwise they remain unprotected and re-enter the susceptible pool (S).
](figure/import/SEITL.png)

The SEITL model can be described with five states (S, E, I, T and L) and five parameters:

1. basic reproductive number ($R_0$)
2. latent period ($D_\mathrm{lat}$)
3. infectious period ($D_\mathrm{inf}$)
4. temporary-immune period ($D_\mathrm{imm}$)
5. probability of developing a long-term protection ($\alpha$).

and the following deterministic equations:

$$
\begin{cases}
\begin{aligned}
\frac{dS}{dt} &= - \beta S \frac{I}{N} + (1-\alpha) \tau T\\
\frac{dE}{dt} &= \beta S \frac{I}{N} - \epsilon E\\
\frac{dI}{dt} &= \epsilon E - \nu I\\
\frac{dT}{dt} &= \nu I - \tau T\\
\frac{dL}{dt} &= \alpha \tau T\\
\end{aligned}
\end{cases}
$$

where $\beta=R_0/D_\mathrm{inf}$, $\epsilon=1/D_\mathrm{lat}$, $\nu = 1/D_\mathrm{inf}$, $\tau=1/D_\mathrm{imm}$ and $N = S + I + R$ is constant. 

As you saw in the lecture, there is an analogy between the deterministic equations and the algorithm that performs stochastic simulations of the model, so that usually only the deterministic equations are given. 

Actually, in order to fit the SEITL model to the Tristan da Cunha dataset we need to add one state variable and one parameter to the model:

* The dataset represents daily incidence counts: we need to create a $6^\mathrm{th}$ state variable - called $\mathrm{Inc}$ - to track the (cumulative) number of daily new cases. Assuming that new cases are reported when they become symptomatic and infectious, we have the following equation for the new state variable:
$$
\frac{d\mathrm{Inc}}{dt} = \epsilon E
$$
* The dataset is incomplete: only 85% of the cases were reported. In addition, we need to account for potential under-reporting of asymptomatic cases. We assume that the data were reported according to a Poisson process with reporting rate $\rho$. Since this reporting rate is unknown (we can only presume that it should be below 85% due to reporting errors) we will include it as an additional parameter.

The deterministic and stochastic SEITL models are already implemented as two `fitmodel` objects called `SEITL_deter` and `SEITL_stoch` respectively. These can be loaded into your **R** session by typing:


```r
data(SEITL_deter)
data(SEITL_stoch)
```

__Take 5 min__ to look at the different elements of both models. In particular, you might be interested in how the daily incidence is computed in the function `SEITL_deter$simulate` and how the likelihood function `SEITL_deter$pointLogLike` accounts for under-reporting. 

You should also note that `SEITL_stoch` mainly differs from `SEITL_deter` through the `simulate` function, which replaces the deterministic equations solver by a stochastic simulation algorithm. More precisely, this algorithm takes a list of transitions (`SEITL_transitions`) and a function to compute the transition rate (`SEITL_rateFunc`). Make sure you understand all the transitions and rates and see the analogy with the deterministic equations.

If you are curious about how the stochastic model is simulated, you can have a look at the code of the function `simulateModelStochastic`. You will note that it calls the function `ssa.adaptivetau` of the **R** package `adaptivetau`, and then process the returned data frame in order to extract the state of the model at the desired observation times given by `times`.

However, before doing any simulation, we need some parameter values and initial state for the model.

## Parameter values and initial state

Based on the description of the outbreak above and the informations below found in the literature, can you think of one or more set(s) of values for the parameters (`theta`) and initial state (`init.state`) of the model? 

1. The $R_0$ of influenza is commonly estimated around 2. However, it can be significantly larger ($R_0>10$) in close-knit communities with exceptional contact rates.
2. Both the average latent and infectious periods for influenza have been estimated around 2 days each.
3. The down-regulation of the cellular response is completed 15 days after symptom onset on average.
4. Serological surveys have shown that the seroconversion rate to influenza (probability to develop specific antibodies) is around 80%. However, it is likely that __not__ all seroconvereted individuals acquire a long-term protective humoral immunity.
5. Between 20 and 30% of the infection with influenza are asymptomatic.
6. There is very limited cross-immunity between influenza viruses A/H1N1 and A/H3N2.

# Deterministic vs Stochastic simulations

Now, let's assess whether the SEITL model can reproduce the two-wave outbreak of Tristan da Cunha with your guess values for the initial state and the parameters. Here we propose to start with the deterministic model and then move on to the stochastic model but feel free to proceed in a different way.

## Deterministic simulations

You can use the `plotFit` function you encountered in the [first practical](introduction.md#assess-the-model-fit) to generate an observation trajectory and display it against the data:


```r
theta.bad.guess <- c(R0 = 2, D.lat = 2, D.inf = 2, alpha = 0.9, D.imm = 13, 
    rho = 0.85)
init.state.bad.guess <- c(S = 250, E = 0, I = 4, T = 0, L = 30, Inc = 0)
plotFit(SEITL_deter, theta.bad.guess, init.state.bad.guess, data = FluTdC1971)
```

<img src="figure/SEITL_model/oneTraj.png" title="" alt="" style="display: block; margin: auto;" />

You can display all states variables (not only the sampled observations) by passing `all.vars=TRUE`


```r
plotFit(SEITL_deter, theta.bad.guess, init.state.bad.guess, data = FluTdC1971, 
    all.vars = TRUE)
```

<img src="figure/SEITL_model/allVars.png" title="" alt="" style="display: block; margin: auto;" />

Although the simulation of the trajectory is deterministic, the observation process is stochastic, hence the noisy `obs` time-series.
You can appreciate the variability of the observation process by plotting several replicates (use the argument `n.replicates`):


```r
plotFit(SEITL_deter, theta.bad.guess, init.state.bad.guess, data = FluTdC1971, 
    n.replicates = 100)
```

<img src="figure/SEITL_model/manyObs.png" title="" alt="" style="display: block; margin: auto;" />

By default, this function plots the mean, median as well as the $95^\mathrm{th}$ and $50^\mathrm{th}$ percentiles of the replicated simulations. Alternatively, you can visualise all the simulated trajectories by passing `summary=FALSE`:


```r
plotFit(SEITL_deter, theta.bad.guess, init.state.bad.guess, data = FluTdC1971, 
    n.replicates = 100, summary = FALSE)
```

<img src="figure/SEITL_model/allTraj.png" title="" alt="" style="display: block; margin: auto;" />

Now, __take 10 min__ to explore the dynamics of the model for different parameter and initial state values. In particular, try different values for $R_0\in[2-15]$ and $\alpha\in[0.3-1]$. For which values of $R_0$ and $\alpha$ do you get a decent fit?

If you didn't manage to get the two waves, you can have look at [our solution](TODO).

You should find that it is hard to reproduce the low incidence in the data with the model. We will now test if this can be better reproduced with a stochastic model which models explicitly the discrete nature of individuals.

## Stochastic simulations


Then __take 10 min__ to explore the dynamics of the stochastic SEITL model with the function `plotFit`. Note that `SEITL_stoch` has the same `theta.names` and `state.names` as `SEITL_deter`, which means that you can use the same parameters and initial state vectors as in the previous section. In addition, you can plot the time-series of the proportion of faded out simulations by passing `p.extinction=TRUE` to `plotFit`.

What differences do you notice between stochastic and deterministic simulations? Conclude on the role of demographic stochasticity on the dynamics of the SEITL model.

# Exponential vs Erlang distributions

So far, we have assumed that the time spent in each compartment was following an [exponential distribution](http://en.wikipedia.org/wiki/Exponential_distribution). This distribution has the well known property of being memoryless, which means that the probability of leaving a compartment doesn't depend on the time already spent in this compartment. Although mathematically convenient, this property is not realistic for many biological processes such as the contraction of the cellular response. In order to include a memory effect, it is common to replace the exponential distribution by an [Erlang distribution](http://en.wikipedia.org/wiki/Erlang_distribution). This distribution is parametrised by its mean $m$ and shape $k$ and can be modelled by $k$ consecutive sub-stages, each being exponentially distributed with mean $m/k$. As illustrated below the flexibility of the Erlang distribution ranges from the exponential ($k=1$) to Gaussian-like ($k>>1$) distributions.

<img src="figure/SEITL_model/erlang.png" title="" alt="" style="display: block; margin: auto;" />

We can extend the SEITL as follows:

![The SEIT2L model extends the SEITL model to account for memory effect in the contraction of the cellular response. The time spend in the $T$ compartment (temporary protection) follow an Erlang distribution with mean $D_\mathrm{imm}=1/\tau$ and shape equal to 2.
](figure/import/SEIT2L.png)

The deterministic and stochastic SEIT2L models are already implemented as `fitmodel` objects, which can be loaded into your **R** session by typing:


```r
data(SEIT2L_deter)
data(SEIT2L_stoch)
```
__Take 5 min__ to have a look at the function `simulate` of these SEIT2L models and make sure you understand how the Erlang distribution for the $T$ compartment is coded.

Now, __take 10 min__ to compare the dynamics of the SEITL and SEIT2L model using `plotFit`. Use your best guess from the previous exercise as well as stochastic and deterministic simulations. Note that although SEITL and SEIT2L share the same parameters, their state variables differ so you need to modify the `init.state` of SEIT2L accordingly.

Can you notice any differences? If so, which model seems to provide the best fit? Do you understand why the shape of the epidemic changes as you change the distribution of the $T$ compartment?


# To go further

* The Poisson process: see bottom of page 3 of this [reference](http://data.princeton.edu/wws509/notes/c4.pdf) for more details.
* Modify the SEIT2L model

<div>
# Navigate
Previous: [MCMC diagnostics](mcmc_diagnostics.md) Next: [MCMC and model comparison](mcmc_and_model_comparison.md)
</div>
