clim_vals_fun <- function(clim.sd, rho, timesteps) {
  
  beta <- sqrt(1 - rho^2)
  e_old <- rnorm(1, mean = 0, sd = clim.sd)
  
  # Use sapply to iteratively generate values
  out_vec <- sapply(1:timesteps, function(i) {
    z_t <- rnorm(1)
    e_new <- rho * e_old + clim.sd * beta * z_t
    e_old <<- e_new  # Use <<- to update e_old in the outer environment
    return(e_new)
  })
  
  return(out_vec)
}
