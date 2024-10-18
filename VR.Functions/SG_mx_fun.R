#' Build survival/growth matrix that doens't yet include fecundity
#' 
#' 
#' @param surv.rates description
#' @param growth.rates description
#' 

sg_mx_fun = function(surv.rates,growth.rates){
  
  # Initiate empty matrix
  mx = matrix(data=0,nrow=3,ncol=3,dimnames = list(c("seed","juv","flrg"),c("seed","juv","flrg")))
  
  ##### Seeds #####
  # Seedlings survive and don't transition to juveniles
  # mx["seed","seed"] = as.numeric(surv.rates["seed"]*(1-growth.rates["seed"]))
    
  ## survive and grow to juv
  mx["juv","seed"] = as.numeric(surv.rates["seed"]*growth.rates["seed"])
  #####
  
  ##### Juvenile #####
  # survive but don't grow
  mx["juv","juv"] = as.numeric(surv.rates["juv"]*(1-growth.rates["juv"]))
  
  # survive and grow
  mx["flrg","juv"] = as.numeric(surv.rates["juv"]*(growth.rates["juv"]))
  #####
  
  ##### Adults ######
  # survive 
  mx["flrg","flrg"] = as.numeric(surv.rates["flrg"])
  return(mx)
  
}