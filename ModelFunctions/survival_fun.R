#' Survival function
#' Creates vectors to but used to build 3x3 matrices for sub populations
#' 
#' @param sizes A vector of names elements indicating the average sizes of y, j, and a size classes
#' @param survival.rates Average survival rates for y, j, and a size classes 
#' @param S_c Not sure about description yet
#' @param S_d not sure about description yet
#' @param clim average macroclimate?. Default = 1
#' @param asp.effect Parameter that shifts macroclimate value to be cooler or warmer
#' @param asp.mag Magnitude of the aspect.effect
#' 

survival_fun = function(sizes, max.survival.rates, min.survival.rates,S_2,S_3,clim=0,asp.effect=0,asp.mag=0){

  # list2env(surv_params,globalenv())
  survs_lm = suppressWarnings(glm(max.survival.rates~sizes,family = binomial( link = "logit" )))
  
  S_0 = survs_lm$coefficients[1]
  S_1 = survs_lm$coefficients[2]
  
  s.rates = logit2prob(S_0 + S_1*sizes + S_2*(clim + asp.mag*asp.effect) + S_3*(clim + asp.mag*asp.effect)^2)
  
  # New to make logistic outputs floor > 0
  s.rates = min.survival.rates + s.rates*(max.survival.rates-min.survival.rates)
  
  return(s.rates)
  }
