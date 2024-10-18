## Building movement matrices ##
# source("VRfunctions/growth_fun.R")
# source("VRfunctions/survival_fun.R")
# source("VRfunctions/SG_mx_fun.R")
# source("VRfunctions/logit2prob.R")

rm(list = ls())
source("VRfunctions/ParamVals.R")
source("VRfunctions/subpop_mx_fun.R")


## Movement matrix ##
n.subpops = 3
move.mx = matrix(data = 0,nrow=n.subpops*3, ncol=n.subpops*3)
start.cells = seq(1,(1+2*n.subpops),by=3)

for(i in start.cells){
  
  # get current start cell. This is the row that will be filled in in this loop
  
  # Get the other start cell values. These columns will be filled in in this loop
  other = setdiff(start.cells,i)
  
  for(j in other){
    
    move.mx[i,j] = round(runif(1,min=0,max=0.1),2)
    
  }
  move.mx[i,i] = 1-sum(move.mx[i,])
} 
diag(move.mx)[diag(move.mx)==0]=1
###

### Population matrix ###
# Build the full population mx
asp.effect = 1
pop.mx = matrix(data = 0,nrow=n.subpops*3, ncol=n.subpops*3)
for(n in 1:n.subpops){
  
  # Add subpop asp.effect to parameters
  surv_params$asp.effect = growth_params$asp.effect = asp.effect
  
  # Create sub-population matrix
  subpop.mx = subpop_mx_fun(surv_params = surv_params,
                            growth_params = growth_params,
                            seeds = 1000)
  
  
  
  
  # put subpopulation MX into full population MX
  start.cell = start.cells[n]
  pop.mx[start.cell:(2+start.cell),start.cell:(2+start.cell)]=subpop.mx

  rm(subpop.mx)
} # End subpop loop

pop.mx = round(pop.mx,2);pop.mx

pop.vec = rep(100,9)
move.mx
pop.mx
move.mx%*%pop.mx

