#' Wrapper function to create matrix
#' @param surv_fun_params Input values need for surv_fun
#' @param growth_fun_params Input values need for growth_fun
#' 

subpop_mx_fun = function(surv_params,growth_params,seeds){
  
  source("VR.Functions/logit2prob.R")
  source("VR.Functions/survival_fun.R")
  source("VR.Functions/growth_fun.R")
  source("VR.Functions/subpop_mx_fun.R")
  
  surv.rates <- do.call(survival_fun,surv_params)
  growth.rates <- do.call(growth_fun,growth_params)
  
  # list2env(surv_params,envir = environment());list2env(growth_params,envir = environment())
  
  # Initiate empty matrix
  mx = matrix(data=0,nrow=3,ncol=3,dimnames = list(c("seed","juv","flrg"),c("seed","juv","flrg")))
  
  ##### seeds #####
  # survive but don't grow
  # mx["seed","seed"] = surv.rates["seed"]*(1-growth.rates["seed"])
  
  ## survive and grow
  mx["juv","seed"] = surv.rates["seed"]*growth.rates["seed"]

  ##### Juvenile #####
  # survive but don't grow
  mx["juv","juv"] = surv.rates["juv"]*(1-growth.rates["juv"])
  
  # survive and grow
  mx["flrg","juv"] = surv.rates["juv"]*growth.rates["juv"]
  #####
  
  ##### Adults ######
  # survive 
  mx["flrg","flrg"] = surv.rates["flrg"]
  
  ##### Fecundity #####
  mx["seed","flrg"] <- surv.rates["flrg"]*P.flr*seeds
  
  # out.list = list(seeds=seeds,surv.rates = surv.rates, growth.rates = growth.rates,subpop_mx=mx)
  #return(out.list)
  return(mx)
  ######
  
}

