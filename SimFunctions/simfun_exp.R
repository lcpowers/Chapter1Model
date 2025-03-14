simfun_exponential = function(burn.in.yrs, sim.yrs, clim.sd, move.mx, asp.mag, asp.effects, seeds, rho){
    
    # total number of years to run the loop
    years = burn.in.yrs + sim.yrs
    
    n.subpops = length(asp.effects)
    
    # Full pop mx start cell positions
    start.cells = seq(1,n.subpops*3,by=3)
    
    # Initiate output data.frame
    output_df = data.frame(year = 1:years,
                           clim.sd = clim.sd,
                           asp.mag = asp.mag,
                           rho = rho,
                           clim = NA,
                           pop.lam = NA)
    output_df[paste0("s_",asp.effects)] <- NA
    output_df[paste0("size_",asp.effects)] <- NA
    
    # Define model params
    surv_params.n = surv_params
    growth_params.n = growth_params
    
    # Add the strength of the aspect effect 
    surv_params.n$asp.mag = growth_params.n$asp.mag = asp.mag
    
    # Set the starting population vectors
    start.size = 1000
    pop.vec.t0 = rep(start.size/(3*n.subpops),3*n.subpops)

    e_old <- rnorm(1, mean=0, sd=clim.sd)
    beta <- sqrt(1-rho^2) # rho is a function parameter
    
    # years = 10
    for(yr in 1:years){
      
      # Generate new value based on previous value and random noise
      z_t <- rnorm(1)

      e_new <- rho*e_old + clim.sd*beta*z_t
      
      # add climate value to output df
      output_df$clim[yr] = e_new

      # Add the climate value to survival and growth function parameter lists
      surv_params.n$clim = growth_params.n$clim = e_new
      e_old <- e_new;rm(e_new)
      
      # Initiate a full population MX for current year
      pop.mx = matrix(data=0,nrow=n.subpops*3,ncol=n.subpops*3)
    
      for(n in 1:n.subpops){
        
        # Add subpop asp.effect to parameters
        surv_params.n$asp.effect = growth_params.n$asp.effect = asp.effects[n]
        
        # Create sub-population matrix
        subpop.mx = subpop_mx_fun(surv_params = surv_params.n,
                                  growth_params = growth_params.n,
                                  seeds.in = seeds)
        
        # Add to main matrix
        start.cell = start.cells[n]
        pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx
        
        output_df[yr,paste0("s_",asp.effects[n])] <- Re(eigen(subpop.mx)$values[1])
        rm(subpop.mx)
        
      } # End subpop loop

      # Multiply movement and pop matrices
      m.pop.mx = move.mx%*%pop.mx
      
      # Calc t1 population vector
      pop.vec.t1 = m.pop.mx%*%pop.vec.t0
      
      # Get annual lambda for full population and store in output dataframe
      ann.lam = sum(pop.vec.t1)/sum(pop.vec.t0)
      output_df$pop.lam[yr]=ann.lam
      
      # Now do that for subpops
      for(n in 1:n.subpops){
        
        start.cell = start.cells[n]
        
        t1sum <- sum(pop.vec.t1[start.cell:(2+start.cell)])
        t0sum <- sum(pop.vec.t0[start.cell:(2+start.cell)])
        
        if(t1sum==0){
          
          subpop.lam <- 0
          
        } else{
          subpop.lam <- t1sum/t0sum
          
        }
        
        # output_df[yr,paste0("s_",asp.effects[n])] <- subpop.lam
        output_df[yr,paste0("size_",asp.effects[n])] <- t1sum
        
      }
      
      # Rescale the population vector
      pop.vec.t1 = pop.vec.t1/sum(pop.vec.t1)*start.size
      
      pop.vec.t0 = pop.vec.t1
      rm(pop.vec.t1)
      
    } # End year loop
    return(output_df)
  }
