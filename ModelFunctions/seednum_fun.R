
basic_seeds_fun <- function(seeds){

  sg.mx[1,3] <- as.numeric(seeds*P.flr*surv.rates["flrg"])
  lam <- Re(eigen(sg.mx)$values[1])
  x = (1.00005-lam)^2
  return(x)

  }


#### This finds the appropriate seed number for the MEGA MATRIX based on a given climate standard deviation #####
complex_seeds_fun <- function(seeds, clim.sd, asp.mag, years){

  # # Full pop mx start cell positions
  # start.cells = F.rows = seq(1,n.subpops*3,by=3)
  # 
  # annual.df <- data.frame(yr=1:years,
  #                         lam=NA)
  # 
  # growth_params.n <- growth_params
  # surv_params.n <- surv_params
  #
  # growth_params.n$asp.mag <- surv_params.n$asp.mag <- asp.mag
  
  for(yr in 1:years){
  
    # clim.yr = rnorm(n = 1,mean = 0,sd = clim.sd)
    
    growth_params.n$clim <- surv_params.n$clim <- clim.yr
    
    # Initiate a full population MX for current year
    pop.mx <- matrix(data=0,nrow=n.subpops*3,ncol=n.subpops*3)
    
    # Build the full population mx
    for(n in 1:n.subpops){
  
      # Add subpop asp.effect to parameters
      surv_params.n$asp.effect <- growth_params.n$asp.effect <- asp.effects[n]
  
      # Create sub-population matrix
      subpop.mx = subpop_mx_fun(surv_params = surv_params.n,
                                growth_params = growth_params.n,
                                seeds.in = 0)
      
      subpop.mx["seed","flrg"] <- surv.rates["flrg"]*P.flr*1002
      
      start.cell = start.cells[n]
      pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx
      rm(subpop.mx)
    }
    lam <- Re(eigen(pop.mx)$values[1])
    annual.df$lam[yr] <- lam
  }

  annual.df2 <- annual.df %>% filter(yr>burn.in)
  lt.lam <- geoMean(annual.df2$lam)
  x = (1.00005-lt.lam)^2
  
  return(x)
}
