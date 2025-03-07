simulate_pop_data <- function(surv_params, growth_params, seeds,
                              timesteps, clim_vals, asp.effects, asp.mag, move.rate) {
  
  # Set the microhabitat effect parameters
  growth_params$asp.mag <- surv_params$asp.mag <- asp.mag # Shift distance
  n.subpops <- length(asp.effects)
  
  # New pop df for each subpop
  seed_counts <- rep(0, n.subpops) # ADDED
  subpop.list.t0 <- list()
  
  for(n in 1:n.subpops){
    
    growth_params$asp.effect <- surv_params$asp.effect <- asp.effects[n]
    
    ## Finding the Subpop specific SSD for when asp.effects asymmetrical
    mean.mx <- subpop_mx_fun(surv_params = surv_params, growth_params = growth_params, seeds.in = seeds)
    ssd.props <- Re(eigen(mean.mx)$vectors[,1])/Re(sum(eigen(mean.mx)$vectors[,1]))
    start_pop <- round(ssd.props * (50/ssd.props[3]) )

    juv_flrg_pop <- sum(start_pop[2:3])
    dt.n <- data.table(
      individual_id = seq_len(juv_flrg_pop),
      start.stage = NA,
      end.stage = rep(c("juv", "flrg"), start_pop[2:3]),
      survival_status = 1,
      growth_status = 0,
      repro_status = 0,
      timestep = 0,
      seeds = 0,
      clim = -9999,
      asp.mag = asp.mag,
      asp.effect = asp.effects[n]
    ) 

    dt.n$seeds[dt.n$end.stage=="flrg"][1] <- start_pop[1]
    subpop.list.t0[[n]] <- dt.n
    rm(dt.n)
    seed_counts[n] <- start_pop[1]
    
  }
  
  all_population_data <- list()
  all_population_data[[1]] <- data.table::rbindlist(subpop.list.t0)

  # Start to loop through timesteps
  for (t in 1:timesteps) {
    # t = 1
    ## initiate subpop list
    subpop.list.t1 <- list()
    new_seed_counts <- rep(0,n.subpops)

    # Put clim value into the survival and growth parameter lists
    clim.t <- clim_vals[t]
    growth_params$clim <- surv_params$clim  <- clim.t
    
    for(n in 1:n.subpops){
      
      # Add aspect effect to growth and surv params
      growth_params$asp.effect <- surv_params$asp.effect <- asp.effects[n] # Shift direction
      
      # Get last timestep dt for subpop n
      subpop.n <- subpop.list.t0[[n]]
      
      # Filter for living individuals
      subpop.n <- subpop.n[survival_status == 1]
      if (nrow(subpop.n) == 0) subpop.n <- data.table()
      
      # Process prior year's seeds into juveniles
      seed_count <- seed_counts[n] # ADDED
      if (seed_count > 0) {

        # Seed survival and growth rates 
        P_surv <- do.call(survival_fun, surv_params)["seed"]
        P_grow <- do.call(growth_fun, growth_params)["seed"]
        
        # Number of new juveniles this year from seeds last year
        n_new_juvs <- rbinom(1, seed_count, P_surv*P_grow)
        
        # If there are new juveniles, add them to the data.table
        if (n_new_juvs > 0) {
          new_juvs_dt <- data.table(
            individual_id = if (nrow(subpop.n) > 0) max(subpop.n$individual_id, na.rm = TRUE) + 1:n_new_juvs else 1:n_new_juvs,
            start.stage = "seed",
            end.stage = "juv",
            survival_status = 1,
            growth_status = 1,
            repro_status = 0,
            timestep = t,
            seeds = 0,
            clim = clim.t,
            asp.mag = asp.mag,
            asp.effect = asp.effects[n]
          )
          # subpop.n <- rbind(subpop.n, new_juvs_dt)
        } # End adding new juveniles
        
      } # End processing last years seeds

      # If there are still living individuals to track
      if (nrow(subpop.n) > 0) {
        
        # Set last year's end.stage to this years start.stage
        subpop.n[, `:=` (start.stage = end.stage, end.stage = NA, timestep = t,
                         survival_status=0, growth_status=0, repro_status=0,
                         seeds = 0, clim = clim.t)]
        
        # Update survival and growth rates based on the new clim
        surv.rates <- do.call(survival_fun, surv_params)
        growth.rates <- do.call(growth_fun, growth_params)
        
        # Used if tracking every individual seed
        # subpop.n[start.stage == "seed", survival_status := rbinom(.N, 1, surv.rates["seed"])]
        # subpop.n[start.stage == "seed" & survival_status == 1, growth_status := rbinom(.N, 1, growth.rates["seed"])]
        # subpop.n[start.stage == "seed" & growth_status == 1, end.stage := "juv"]
        
        subpop.n[start.stage == "juv", survival_status := rbinom(.N, 1, surv.rates["juv"])]
        subpop.n[start.stage == "juv" & survival_status == 1, growth_status := rbinom(.N, 1, growth.rates["juv"])]
        subpop.n[start.stage == "juv" & growth_status == 1, end.stage := "flrg"]
        subpop.n[start.stage == "juv" & survival_status == 1 & growth_status == 0, end.stage := "juv"]
        subpop.n[start.stage == "juv" & survival_status == 0, end.stage := "dead"]
        
        subpop.n[start.stage == "flrg", survival_status := rbinom(.N, 1, surv.rates["flrg"])]
        subpop.n[start.stage == "flrg" & survival_status == 1, repro_status := rbinom(.N, 1, P.flr)]
        subpop.n$seeds[subpop.n$start.stage == "flrg" & subpop.n$repro_status == 1] = seeds 
        subpop.n[start.stage == "flrg" & survival_status == 1, end.stage := "flrg"]
        subpop.n[start.stage == "flrg" & survival_status == 0, end.stage := "dead"]
        
        # Fill in empty end.stage cells
        # subpop.n[is.na(end.stage) & survival_status == 1 & growth_status == 0, end.stage := start.stage]
        # subpop.n[survival_status == 0, end.stage := "dead"]
        
        new_seed_counts[n] <- sum(subpop.n$seeds)
         
      }
      
    if( exists("new_juvs_dt") ) {subpop.n <- rbind(subpop.n,new_juvs_dt)}
      
    subpop.list.t1[[n]] <- subpop.n
      
    } # End subpopulation fates dataframe
    
    ### Dealing with seed movement
    
    # Add new seeds to the population dataframe
    if (sum(new_seed_counts) > 0) {
      
      # movement matrix to adjust seed vector
      move.mx <- matrix(data = NA,
                        nrow = n.subpops,
                        ncol = n.subpops)
      diag(move.mx) <- 1 - move.rate # Along diagonal is what stays
      move.mx[is.na(move.mx)] <- move.rate/(n.subpops-1) # Other elements are what moves
      new_seed_counts_move <- as.numeric(move.mx%*%new_seed_counts)
      new_seed_counts <- new_seed_counts_move;rm(new_seed_counts)
    }
    
    # Turn subpop list in dataframe and add to output list
    all_subpops <- rbindlist(subpop.list.t1, use.names = TRUE, fill = F)
    all_population_data[[t + 1]] <- all_subpops
    
    # Set this years subpops equal to t0 and remove t1
    subpop.list.t0 <- subpop.list.t1;rm(subpop.list.t1)
    
    # *** CHECK FOR ZERO TOTAL POPULATION HERE ***
    # If theres no one left, break the loop/stop the sim
    if (nrow(all_subpops) == 0 && sum(seed_counts) == 0) break
    rm(all_subpops)
    
  } # End timestep loop
  
  all_population_output <- rbindlist(all_population_data, use.names = TRUE, fill = TRUE)
  
  if (max(all_population_output$timestep) < timesteps) {
    extinct_df <- data.table(
      individual_id = -9999,
      start.stage = "extinct",
      end.stage = NA,
      survival_status = 0,
      growth_status = 0,
      repro_status = 0,
      timestep = seq(max(all_population_output$timestep) + 1, timesteps),
      seeds = 0,
      clim = NA,
      asp.mag = asp.mag,
      asp.effect = asp.effects
    )
    
    all_population_output <- rbind(all_population_output, extinct_df, fill = TRUE)
  
    }
  return(all_population_output)
} # End function
