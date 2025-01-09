moveMx_fun = function(n.subpops,move.rate,symmetrical=TRUE){
  
  ## Creating all combinations of movement rates
  subpop.names = LETTERS[1:n.subpops] # Name subpops with letters
  
  move.dirs = expand.grid(from=subpop.names,to=subpop.names) %>% # get all combos
    # Remove rows that imply no movement (A to A)
    filter(from!=to) %>% 
    # Concatenate
    mutate(col.names = paste0(from,"to",to)) %>% 
    # save last column as a vector
    arrange(col.names) %>% pull() 
  
  # Create movement data.frame to loop over
  move.rates.df = matrix(data=rep(move.rate,length(move.dirs)),ncol=length(move.dirs)) %>% # generate matrix with movement rates
    as.data.frame() %>% 
    # Make all combos of movement rates
    expand.grid()
  
  if(symmetrical==T){
    
    ## This filters for symmetrical movement rates
    move.rates.df <- move.rates.df %>% filter(if_all(starts_with("V"),~ .==V1))
    
  }
  
  names(move.rates.df)=move.dirs # rename columns to make movement direction clear
  rm(subpop.names,move.dirs,move.rate) 
  
  # Tell me if move rates != number of subpops
  # stopifnot("Number of movement rates not equal to number of subpops"=length(move.rates)==n.subpops)
  
  # Initiate a zero matrix
  move.mx = matrix(data = 0,nrow=n.subpops*3, ncol=n.subpops*3)
  
  # get the upper L cell of each subpop within the larger pop mx
  start.cells = seq(1,(1+2*n.subpops),by=3)
  
  # Each iteration of this loop corresponds to the column being filled in
  # First I fill in the dispersal cells (from i to others) then find that non-dispersal cell
  for(j in 1:n.subpops){
    
    # 'j' is the current start cell. This is the column that will be filled in in this loop
    
    # Get the indexes of other subpops. These will be the "FROM" subpops in the inner loop
    other.subpops = setdiff(1:n.subpops,j)
    
    for(k in other.subpops){
      
      # From the list of movement rates, get the from I to J movement value
      move.rate.k = move.rates.df[[paste0(LETTERS[j],"to",LETTERS[k])]] %>% unlist()
      
      # Row k is the pop moving to
      # Col j is the pop moving from
      
      move.mx[start.cells[k],start.cells[j]] = move.rate.k
      
    }
    move.mx[start.cells[j],start.cells[j]] = 1-sum(move.mx[,start.cells[j]])
  } 
  
  ## add ones to the position along the diagonal that are still = 0
  diag(move.mx)[diag(move.mx)==0]=1 
  return(move.mx)
}