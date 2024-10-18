
basic_seeds_fun <- function(seeds){

  sg.mx[1,3] <- as.numeric(seeds*P.flr*surv.rates["flrg"])
  lam <- Re(eigen(sg.mx)$values[1])
  x = (1.00005-lam)^2
  return(x)

  }

complex_seeds_fun <- function(seeds){
  
  # Full pop mx start cell positions
  start.cells = F.rows = seq(1,n.subpops*3,by=3)
  
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
    
    ## Adjust fecundity cells based on movement rates
    # Get original F cell value
    f.cell.val = subpop.mx[1,3]
    
    # Reduce F cell based on movement rate out of this subpop
    subpop.mx[1,3] = f.cell.val*(1-move.rate)
    
    # put subpopulation MX into full population MX
    start.cell = start.cells[n]
    pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx
    
    ## Now add values to cells for movement from this subpop to others
    # Get rows where seeds are "moving to"
    move.to.rows = setdiff(start.cells,start.cells[n])
    
    # Fec. cell for
    F.col.n = start.cells[n]+2
    
    # Original fecundity cell value * fraction that is moving / the number of places moving to. 
    pop.mx[move.to.rows,F.col.n] = f.cell.val*move.rate/(n.subpops-1)
  }
  
  lam <- Re(eigen(pop.mx)$values[1])
  x = (1.00005-lam)^2
  return(x)
}
