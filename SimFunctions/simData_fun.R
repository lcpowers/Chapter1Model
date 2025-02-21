simulate_pop_data <- function(surv_params, growth_params, seeds,
                              timesteps, clim_vals, asp.effects, asp.mag, move.rate) {
  
  # Set the microhabitat effect parameters
  growth_params$asp.mag <- surv_params$asp.mag <- asp.mag # Shift distance

  # New pop df for each subpop
  n.subpops <- length(asp.effects)
  subpop.list.t0 <- list()
  
  for(n in 1:n.subpops){
    
    growth_params$asp.effect <- surv_params$asp.effect <- asp.effects[n]

    ## Finding the Subpop specific SSD for when asp.effects asymmetrical
    mean.mx <- subpop_mx_fun(surv_params = surv_params, growth_params = growth_params, seeds.in = seeds)

    ssd.props <- eigen(mean.mx)$vectors[,1]/sum(eigen(mean.mx)$vectors[,1])
    ssd.scalar <- 50/ssd.props[3]
    start_pop <- round(ssd.props*ssd.scalar)
    print(start_pop)
    total_pop <- sum(start_pop)
    
    # Initialize population 
    pop <- data.frame(
      individual_id = 1:total_pop,
      start.stage = NA,
      end.stage = c(rep("seed",start_pop[1]),rep("juv",start_pop[2]),rep("flrg",start_pop[3])),
      survival_status = rep(1, total_pop),  # Track if they survive (1 = survive, 0 = don't survive)
      growth_status = rep(NA, total_pop),   # Track growth status (NA = not alive, 1 = grew, 0 = didn't grow)
      repro_status = rep(NA, total_pop),    # Track reproduction events to make sure the align with probability
      timestep = rep(0, total_pop),         # Start at timestep 0
      stringsAsFactors = FALSE
    ) %>% 
      mutate(seeds = NA,
             clim = NA,
             asp.mag = asp.mag,
             asp.effect = asp.effects[n])
    
    assign(paste0("subpop.",n),pop)
    subpop.list.t0[[n]] <- eval(as.name(paste0("subpop.",n)))
    rm(pop)
    
  }

  all_subpops <- data.table::rbindlist(subpop.list.t0)
  # Store population data for each timestep in a list
  all_population_data <- list()
  all_population_data[[1]] <- all_subpops;rm(all_subpops)  # Store initial state

  # Start to loop through timesteps
  for (t in 1:timesteps) {

    # print(t)
    # t = 2
    subpop.list.t1 <- list()
    new_seeds_vec <- rep(0,n.subpops)

    for(n in 1:n.subpops){

      # Get subpop.n dataframe
      subpop.n <- subpop.list.t0[[n]] %>%
        mutate(start.stage = end.stage,
               end.stage = NA)

      # *** CHECK FOR ZERO SUBPOPULATION HERE ***
      if (is.null(subpop.n)) {

        # If the output from the. list is NULL go to next loo
        next  # Uncomment to stop the simulation entirely

        } else {

          # Otherwise, filter for survived individuals
        subpop.n <- filter(subpop.n, survival_status == 1)

      }

      # And if noone survived last year
      if (nrow(subpop.n)==0) {

        next  # Carry on to next loop

        }

      # Drop rows of individuals that died in the last time step
      # subpop.n <- filter(subpop.n,survival_status==1)

      # Add aspect effect to growth and surv params
      growth_params$asp.effect <- surv_params$asp.effect <- asp.effects[n] # Shift direction

      # Get clim value from input vector of clim values. Input vector rather than selected here to make sure that the same values are used for different microhabitats/subpops
      clim.t <- clim_vals[t]

      # Put clim value into the survival and growth parameter lists
      growth_params$clim <- surv_params$clim  <- clim.t
      subpop.n$clim <- clim.t # Store this timestep's climate value

      # Update survival and growth rates based on the new clim
      surv.rates <- do.call(survival_fun, surv_params)
      growth.rates <- do.call(growth_fun, growth_params)

      # Process each individual in the subpopulation
      for (i in 1:nrow(subpop.n)) {

        # Update timestep for each individual
        subpop.n$timestep[i] <- t

        if (subpop.n$start.stage[i] == "seed") {

          # This is a very stepwise way to do this...

          # Survival and growth probabilities for seeds
          survival_prob <- surv.rates["seed"]
          growth_prob <- growth.rates["seed"]

          # print(paste("Timestep:", t, "Subpopulation:", n, "Row:", i))
          # print(paste("survival_prob:", survival_prob))

          # Survival status
          subpop.n$survival_status[i] <- rbinom(1, 1, survival_prob)

          # If survived, grow?
          if (subpop.n$survival_status[i] == 1) {

            # Growth (transition to juvenile) status
            subpop.n$growth_status[i] <- rbinom(1, 1, growth_prob)

            # If growth occurs, move to juvenile stage
            if (subpop.n$growth_status[i] == 1) {

              subpop.n$end.stage[i] <- "juv"
              # seed_to_juv = seed_to_juv + 1

            } else if(subpop.n$growth_status[i] == 0){

              # This else if wasn't here before, but without it, seeds can survive without transitioning to juv stage.
              # i.e. there was a seed bank
              subpop.n$end.stage[i] <- NA
              subpop.n$survival_status[i] <- 0

            }

          }

        } else if (subpop.n$start.stage[i] == "juv") {

          # Survival and growth probabilities for juveniles
          survival_prob <- surv.rates["juv"]
          growth_prob <- growth.rates["juv"]

          # Determine survival
          subpop.n$survival_status[i] <- rbinom(1, 1, survival_prob)

          # If survived, grow
          if (subpop.n$survival_status[i] == 1) {
            subpop.n$growth_status[i] <- rbinom(1, 1, growth_prob)

            # If growth occurs, move to repro/flowering stage
            if (subpop.n$growth_status[i] == 1) {
              subpop.n$end.stage[i] <- "flrg"
            } else {

              subpop.n$end.stage[i] <- "juv"

            }
          }

        } else if (subpop.n$start.stage[i] == "flrg") { ## Flowering stage
          # Survival probs flowering plants (no growth possible)
          survival_prob <- surv.rates["flrg"]

          # Determine survival
          subpop.n$survival_status[i] <- rbinom(1, 1, survival_prob)

          # If survived, reproduce and produce new seeds
          if (subpop.n$survival_status[i] == 1) {

            subpop.n$end.stage[i] <- "flrg"

            # Prob reproduction
            subpop.n$repro_status[i] <- rbinom(1, 1, P.flr)

            # If the reproductive-age flower did in fact produce flowers, how many seeds?
            if(subpop.n$repro_status[i] == 1){

              new_seeds <- seeds #rpois(n = 1,lambda = seed_vec)
              subpop.n$seeds[i] <- new_seeds # Store seed number to check later
              # total_new_seeds <- total_new_seeds + new_seeds

            } # End seeds

          } # End repro

        } # End flrg stage

      } # End loop going through each row of subpop dataframe

      subpop.n$end.stage[subpop.n$survival_status==0] = "dead"

      new_seeds_vec[n] = sum(subpop.n$seeds,na.rm=T)
      assign(paste0("subpop.",n),subpop.n)
      rm(subpop.n)

    } # End subpopulation fates dataframe

    ### Dealing with new seeds and seed movement

    # Add new seeds to the population dataframe
    if (sum(new_seeds_vec) > 0) {

      # movement matrix to adjust seed vector
      move.mx <- matrix(data = NA,
                        nrow = n.subpops,
                        ncol = n.subpops)
      diag(move.mx) <- 1 - move.rate # Along diagonal is what stays
      move.mx[is.na(move.mx)] <- move.rate/(n.subpops-1) # Other elements are what moves
      new_seeds_move_vec <- c(move.mx%*%new_seeds_vec)

      for(n in 1:n.subpops){

        new_seeds <- round(sum(new_seeds_move_vec[n]))

        # Individual IDs for new seeds
        subpop.n <- eval(as.name(paste0("subpop.",n)))
        new_individual_ids <-(max(subpop.n$individual_id)+1):(max(subpop.n$individual_id)+new_seeds)

        new_seeds_df <- data.frame(
          individual_id = new_individual_ids,
          start.stage = NA,
          end.stage = rep("seed", new_seeds),
          survival_status = rep(1, new_seeds),
          growth_status = rep(NA, new_seeds),
          repro_status = rep(NA, new_seeds),
          timestep = rep(t, new_seeds),
          clim = rep(clim.t,new_seeds),
          stringsAsFactors = FALSE
        ) %>%
          mutate(seeds = NA,
                 asp.mag = asp.mag,
                 asp.effect = asp.effects[n])

        subpop.n <- rbind(subpop.n,new_seeds_df)
        subpop.list.t1[[n]] <- subpop.n
      }

    } # End adding and moving new seeds

    # Turn subpop list of dataframes into one dataframe
    all_subpops <- data.table::rbindlist(subpop.list.t1)

    # add that to the all population list
    all_population_data[[t+1]] <- all_subpops

    # Set this years subpops equal to t0 and remove t1
    subpop.list.t0 <- subpop.list.t1;rm(subpop.list.t1)

    # *** CHECK FOR ZERO TOTAL POPULATION HERE ***
    if (nrow(all_subpops)==0) {

      # If theres no one left, break the loop/stop the sim
      break
    }

    rm(all_subpops)

  } # End timestep loop

  all_population_output <- data.table::rbindlist(all_population_data, use.names = TRUE, fill = TRUE)

  if(max(all_population_output$timestep)<timesteps){
    extinct_df <- expand.grid(individual_id = -9999,
                              start.stage = "extinct",
                              end.stage = NA,
                              survival_status=0,
                              growth_status=0,
                              repro_status=0,
                              timestep = seq((max(all_population_output$timestep)+1),timesteps,1),
                              seeds = NA,
                              clim = NA,
                              asp.mag = NA,
                              asp.effect = NA)
    all_population_output <- rbind(all_population_output,extinct_df)
  }
  return(all_population_output)
} # End function
