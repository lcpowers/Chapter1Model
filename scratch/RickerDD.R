rm(list=ls())
library(tidyverse)

sapply(list.files("VR.Functions",pattern=".R$",full.names =T), source)
fec_df = read_csv("VR.Functions/fecundity.csv")

surv.rates <- do.call(survival_fun,surv_params)
grow.rates <- do.call(growth_fun,growth_params)

seeds <- fec_df$seeds[fec_df$aspect.effect==1&fec_df$asp.mag==0][1]

mx <- subpop_mx_fun(surv_params = surv_params, growth_params = growth_params, seeds=seeds)
mx
Re(eigen(mx)$values[1])

param_df <- expand.grid(B0 = seq(0,0.1,by=0.025),
                        B1 = seq(0,0.001,by=0.00025))
max_F <- mx[1,3]
years <- 1000





out_df <- NULL

for(i in 1:nrow(param_df)){
  
  pop.vec = rep(500,3)
  tmp.mx <- mx
  tmp.df <- expand.grid(B0 = param_df$B0[i],
                        B1 = param_df$B1[i],
                        year = 0:years,
                        s1 = NA, s2 = NA, s3 = NA, total = NA)
  
  tmp.df[1,c("s1","s2","s3")] <- pop.vec
  tmp.df$total[1] <- sum(pop.vec)
  
  tmp.B0 <- param_df$B0[i]
  tmp.B1 <-  param_df$B1[i]
  
  for(yr in 1:years){
    
    # set F value based on current adult class size
    tmp.mx[1,3] <- max_F*exp(tmp.B0 - tmp.B1*pop.vec[3])
    
    pop.vec <- tmp.mx%*%pop.vec
    
    tmp.df[yr+1,c("s1","s2","s3")] <- pop.vec
    tmp.df$total[yr+1] <- sum(pop.vec)
  }
  
  out_df <- rbind(out_df,tmp.df)
  
}

ggplot(filter(out_df,B1>0&year>10),aes(x=year,y=total,color=as.factor(B1)))+
  geom_line()+
  facet_wrap(~B0)

adult_pop = 100

param_df <- expand.grid(B0 = seq(0,3,by=0.025),
                        B1 = seq(0,0.05,by=0.025)) %>% 
  mutate(F_100 = max_F*exp(-(B0 + B1*100)),
         F_1000 = max_F*exp(-(B0 + B1*1000)))

ggplot(param_df,aes(x=B0,color=factor(B1)))+
  geom_line(aes(y=F_100))+
  ggtitle("100")

ggplot(param_df,aes(x=B0,color=factor(B1)))+
  geom_line(aes(y=F_1000))+
  ggtitle("1000")




# For a set of values, the slope of the line is steeper when base adult population size is smaller





#### OLD #######
# Before DD
mx <- matrix(data=0,nrow=3,ncol=3)

mx[1,1] <- 0.001 
mx[2,2] <- 0.4391621 # J surv
mx[3,3] <- 0.9515129 # A surv
mx[3,2] <- 0.2866531 # J to A
mx[1,3] <- 2.55
mx[2,1] <- 0.03734508

mx
Re(eigen(mx)$values[1])

years = 1000
pop.vec = rep(1000,3)

exp_out_df = data.frame(year = 0:years,
                    s1 = NA,
                    s2 = NA,
                    s3 = NA,
                    total = NA,
                    lam = NA)
exp_out_df[1,2:4] = pop.vec
exp_out_df$total[1] = sum(pop.vec)

for(yr in 1:years){
  
  pop.vec = mx%*%pop.vec
  row.i = which(exp_out_df$year==yr)
  
  exp_out_df[row.i,2:4] = pop.vec
  exp_out_df$total[row.i] = sum(pop.vec) 
  exp_out_df$lam[row.i] = exp_out_df$total[row.i]/exp_out_df$total[row.i-1]
  
}


## With DD
mx2 <- mx

Re(eigen(mx2)$values[1])

max.g = 0.03734508
B0 = 0.35
pop.vec = rep(1000,3)

dd_out_df = data.frame(year = 0:years,
                        s1 = NA,
                        s2 = NA,
                        s3 = NA,
                        total = NA,
                        lam = NA)
dd_out_df[1,2:4] = pop.vec
dd_out_df$total[1] = sum(pop.vec)

for(yr in 1:years){
  
  adult.prop = (pop.vec[2]+pop.vec[3])/sum(pop.vec)
  
  mx2[2,1] = max.g*exp(B0-B1*adult.prop)
  
  pop.vec = mx2%*%pop.vec
  row.i = which(dd_out_df$year==yr)
  
  dd_out_df[row.i,2:4] = pop.vec
  dd_out_df$total[row.i] = sum(pop.vec) 
  dd_out_df$lam[row.i] = dd_out_df$total[row.i]/dd_out_df$total[row.i-1]
  
}

ggplot(filter(dd_out_df,year>10),aes(x=year,y=total))+
  geom_point()







