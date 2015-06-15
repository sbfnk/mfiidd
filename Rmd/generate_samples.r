# @knitr generate_samples
trace <- my_mcmcMH(target = my_dLogPosterior_R0_epi1, # target distribution
                   init.theta = 1, # intial parameter guess
                   proposal.sd = 0.1, # standard deviation of
                       				  # Gaussian proposal: 0.1
                   n.iterations = 1000) # number of iterations
