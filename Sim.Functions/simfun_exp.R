simfun_exponential = function(n.subpops, burn.in.yrs, sim.yrs, clim.sd, move.mx, asp.mag, asp.effects){
    
    # total number of years to run the loop
    years = burn.in.yrs + sim.yrs
    
    # Full pop mx start cell positions
    start.cells = seq(1,n.subpops*3,by=3)
    
    # Initiate output data.frame
    output_df = data.frame(year = 1:years,
                           clim.sd = clim.sd,
                           asp.mag = asp.mag,
                           clim = NA,
                           pop.lam = NA)

    # Define model params
    surv_params.n = surv_params
    growth_params.n = growth_params
    
    # Add the strength of the aspect effect 
    surv_params.n$asp.mag = growth_params.n$asp.mag = asp.mag
    
    # Set the starting population vectors
    pop.size = 1000
    pop.vec.t0 = rep(pop.size,3*n.subpops)
    
    for(yr in 1:years){
      
      # yr = 1
      
      # Climate in this year
      clim.yr = rnorm(n=1, mean=0, sd=clim.sd)
      
      # add climate value to output df
      output_df$clim[yr] = clim.yr
      
      # Add the climate value to survival and growth function parameter lists
      surv_params.n$clim = growth_params.n$clim = clim.yr
      
      # Initiate a full population MX for current year
      pop.mx = matrix(data=0,nrow=n.subpops*3,ncol=n.subpops*3)
      
      # Build the full population mx
      for(n in 1:n.subpops){
        
        # Add subpop asp.effect to parameters
        surv_params.n$asp.effect = growth_params.n$asp.effect = asp.effects[n]
        
        seeds = seeds_df$seeds[seeds_df$asp.mag==asp.mag & seeds_df$aspect.effect == asp.effects[n]]
        # seeds = 40
        # Create sub-population matrix
        subpop.mx = subpop_mx_fun(surv_params = surv_params.n,
                                  growth_params = growth_params.n,
                                  seeds = seeds)
        
        Re(eigen(subpop.mx)$values[1])
        start.cell = start.cells[n]
        pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx
        rm(subpop.mx)
      } # End subpop loop
      
      # Multiply movement and pop matrices
      m.pop.mx = move.mx%*%pop.mx
      
      # Calc t1 population vector
      pop.vec.t1 = m.pop.mx%*%pop.vec.t0
      
      # Get annual lambda for full population and store in output dataframe
      ann.lam = sum(pop.vec.t1)/sum(pop.vec.t0)
      output_df$pop.lam[yr]=ann.lam
      
      # Rescale the population vector. 
      pop.vec.t1 = pop.vec.t1/sum(pop.vec.t1)*pop.size
      
      pop.vec.t0 = pop.vec.t1
      rm(pop.vec.t1)
      
    } # End year loop
    return(output_df)
  }
