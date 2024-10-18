#' Growth function
#' Creates vectors to but used to build 3x3 matrices for sub populations
#' 
#' @param sizes A vector of names elements indicating the average sizes of y, j, and a size classes
#' @param growth.rates Average survival rates for y, j, and a size classes 
#' @param G_c Not sure about description yet
#' @param G_d not sure about description yet
#' @param clim average regional climate
#' @param asp.effect Parameter that shifts macroclimate value to be cooler or warmer
#' @param asp.mag Magnitude of the aspect.effect
#' 
# Logistic function to modify growth

growth_fun = function(sizes,min.growth.rates,max.growth.rates,G_2,G_3,clim=0,asp.effect=0,asp.mag=0){
  
  # list2env(growth_params,globalenv())
  grow_lm = suppressWarnings(glm(max.growth.rates~sizes[c("seed","juv")],family = binomial( link = "logit" )))
  G_0 = grow_lm$coefficients[1]
  G_1 = grow_lm$coefficients[2]

  g.rates = logit2prob(G_0 + G_1*sizes[c("seed","juv")] + G_2*(clim + asp.mag*asp.effect) + G_3*(clim + asp.mag*asp.effect)^2)
  
  g.rates = min.growth.rates + g.rates*(max.growth.rates-min.growth.rates)

  return(g.rates)
}

