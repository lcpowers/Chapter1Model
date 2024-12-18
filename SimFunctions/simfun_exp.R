simfun_exponential = function(n.subpops, burn.in.yrs, sim.yrs, clim.sd, move.mx, asp.mag, asp.effects,seeds,rho){
    
    # total number of years to run the loop
    years = burn.in.yrs + sim.yrs
    
    # Full pop mx start cell positions
    start.cells = seq(1,n.subpops*3,by=3)
    
    # Initiate output data.frame
    output_df = data.frame(year = 1:years,
                           clim.sd = clim.sd,
                           asp.mag = asp.mag,
                           rho = rho,
                           clim = NA,
                           pop.lam = NA)
    output_df[paste0("s",1:n.subpops)] <- NA
    
    # Define model params
    surv_params.n = surv_params
    growth_params.n = growth_params
    
    # Add the strength of the aspect effect 
    surv_params.n$asp.mag = growth_params.n$asp.mag = asp.mag
    
    # Set the starting population vectors
    pop.size = 1000
    pop.vec.t0 = rep(pop.size,3*n.subpops)
    
    # e_old <- rnorm(1, mean=0, sd=clim.sd)
    
    if(rho==0){
      clims = rnorm(years,mean=0,sd=clim.sd)
    } else if(rho>0){
      clims <- arima.sim(n = years, list(ar = rho), sd = clim.sd)
    }
    
    for(yr in 1:years){
      
      # # yr = 1
      # 
      # z_t <- rnorm(1, mean = 0, sd = clim.sd)
      # 
      # e_new <- rho * e_old + z_t
      # output_df$clim[yr] <- e_new
      # 
      ### Old way to add climate autocorrelation ###
      # Climate in this year
      # Random climate draw for this year
      # z_t <- rnorm(1,mean=0,sd=clim.sd)
      # 
      # M&D equation 4.16 (pg. 139) 
      # e_new <- rho*e_old + clim.sd*((1 - rho^2)^(1/2))*z_t
      ### End old way to add climate autocorrelation ###
      # 
      # add climate value to output df
      output_df$clim[yr] = clims[yr]

      # Add the climate value to survival and growth function parameter lists
      surv_params.n$clim = growth_params.n$clim = clims[yr]
      
      # Make this years adjusted clim last years & Remove this years
      # e_old <- e_new;rm(e_new)
      
      # Initiate a full population MX for current year
      pop.mx = matrix(data=0,nrow=n.subpops*3,ncol=n.subpops*3)
      
      # Build the full population mx
      for(n in 1:n.subpops){
        
        # Add subpop asp.effect to parameters
        surv_params.n$asp.effect = growth_params.n$asp.effect = asp.effects[n]
        
        # Create sub-population matrix
        subpop.mx = subpop_mx_fun(surv_params = surv_params.n,
                                  growth_params = growth_params.n,
                                  seeds = seeds)
        
        # Re(eigen(subpop.mx)$values[1])
        start.cell = start.cells[n]
        pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx
        output_df[yr,paste0("s",n)] <- Re(eigen(subpop.mx)$values[1])
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
