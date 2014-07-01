# @knitr generate_samples
trace <- my_mcmcMH(target.dist = my_posterior_R0, # target distribution
                   init.theta = 1, # intiial parameter (R0) guess
                   proposal.sd = 0.1, # standard deviation of Gaussian proposals
                   n.iterations = 1000) # number of interations
