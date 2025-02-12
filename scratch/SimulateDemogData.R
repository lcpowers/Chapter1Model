# Simulating demography data #
rm(list = ls())

asp.effects <- c(-1,1)
asp.mag <- 2
clim.sd <- 2
n.subpops <- length(asp.effects)
yrs <- 100
seeds <- seeds_df$avg.seeds[seeds_df$asp.effects== paste(asp.effects,collapse = ",")&seeds_df$asp.mag==asp.mag]
rho <- 0.1
move.rate = 0.05

starting.pop.vec <- rep(1000,3)

