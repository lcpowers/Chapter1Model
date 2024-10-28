## ParameterVals
sizes = c(seed=1,juv=8,flrg=10)
P.flr = 0.85 # Probability of flowering
P.germ = 0.01 # Probability of seeds germinating

surv_params = list(sizes=sizes,
                   max.survival.rates = c(seed=0.34,juv=0.75,flrg=0.98),
                   min.survival.rates = c(seed=0.01,juv=0.25,flrg=0.5),
                   S_2 = c(-0.001,-0.001,-0.001),
                   S_3 = c(-0.0075,-0.015,-0.015),
                   asp.mag = 0,
                   asp.effect = 0,
                   clim = 0)

growth_params = list(sizes=sizes,
                     max.growth.rates = c(seed=0.3,juv=0.49),
                     min.growth.rates = c(seed=0.01,juv=0.025),
                     G_2 = c(0.001,0.001),
                     G_3 = c(-0.01,-0.009),
                     asp.mag = 0,
                     asp.effect = 0,
                     clim = 0)

# list2env(surv_params,envir = .GlobalEnv)
