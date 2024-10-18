moveMx_fun = function(n.subpops,move.rates){
  
  # Tell me if move rates != number of subpops
  # stopifnot("Number of movement rates not equal to number of subpops"=length(move.rates)==n.subpops)
  
  # Initiate a zero matrix
  move.mx = matrix(data = 0,nrow=n.subpops*3, ncol=n.subpops*3)
  
  # get the upper L cell of each subpop within the larger pop mx
  start.cells = seq(1,(1+2*n.subpops),by=3)
  
  # Each iteration of this loop corresponds to the column being filled in
  # First I fill in the dispersal cells (from i to others) then find that non-dispersal cell
  for(i in 1:n.subpops){
    
    # 'i' is the current start cell. This is the column that will be filled in in this loop
    
    # Get the indexes of other subpops. These will be the "FROM" subpops in the inner loop
    other.subpops = setdiff(1:n.subpops,i)
    
    for(j in other.subpops){
      
      # From the list of movement rates, get the from I to J movement value
      move.rate.j = move.rates[[paste0(LETTERS[i],"to",LETTERS[j])]] %>% unlist()
      # Row j is the pop moving to
      # Col i is the pop moving from
      
      move.mx[start.cells[j],start.cells[i]] = move.rate.j
      
    }
    move.mx[start.cells[i],start.cells[i]] = 1-sum(move.mx[,start.cells[i]])
  } 
  
  ## add ones to the position along the diagonal that are still = 0
  diag(move.mx)[diag(move.mx)==0]=1 
  return(move.mx)
}