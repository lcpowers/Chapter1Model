# why can't I remember how matrices work
rm(list=ls())
pop = matrix(data = c(0.1, 0, 100,      0,0,0,
                      0.01, 0.5, 0,   0,0,0,
                      0, 0.1, 0.8,    0,0,0,
                      
                      0,0,0,    0.1, 0, 100,
                      0,0,0,    0.01, 0.5, 0,
                      0,0,0,    0, 0.1, 0.8),
             byrow=T,ncol=6)
pop.vec = rep(100,6)

move.mx


pop%*%pop.vec

pop

Re(eigen(pop)$values[1])

a = matrix(data=c(0,0,1,
                  0,0,0),
           byrow=T,
           nrow=2)
b = matrix(data=rep(0,6),nrow=2)

a = matrix(data=round(runif(n=9),1),
           nrow=3)
a

b = matrix(data=rep(1,3))
a%*%a

# 
# fec.mx = matrix(nrow=6,byrow=T,
#                        data=c(0,0,100,0,0,0,
#                               0,0,0,0,0,0,
#                               0,0,0,0,0,0,
#                               0,0,0,0,0,100,
#                               0,0,0,0,0,0,
#                               0,0,0,0,0,0))
# 
# move.mx = matrix(nrow=6,byrow=T,
#                  data=c(0,0,0.5, 0,0,0.5,
#                         0,0,0,   0,0,0,
#                         0,0,0,   0,0,0,
#                         
#                         0,0,0.5, 0,0,0.5,
#                         0,0,0,   0,0,0,
#                         0,0,0,   0,0,0))
# move.mx
# fec.mx%*%move.mx

  
move.mx = matrix(nrow=6,byrow=T,
                 data=c(0.5,0,0, 0.5,0,0,
                        0,  0,0,   0,0,0,
                        0,  0,0,   0,0,0,
                        
                        0.5,0,0, 0.5,0,0,
                        0,  0,0,   0,0,0,
                        0,  0,0,   0,0,0))


fec.mx = matrix(nrow=6,byrow=T,
                data=c(0,0,100,  0,0,0,
                       0,0,0,    0,0,0,
                       0,0,0,    0,0,0,
                       
                       0,0,0,    0,0,100,
                       0,0,0,    0,0,0,
                       0,0,0,    0,0,0))
move.mx%*%fec.mx

#####

n.pops = 2
move.mx = matrix(data = 0,nrow=n.pops*3, ncol=n.pops*3)
start.cells = seq(1,(1+2*n.pops),by=3)

for(i in start.cells){
  
  # get current start cell. This is the row that will be filled in in this loop
  
  # Get the other start cell values. These columns will be filled in in this loop
  other = setdiff(start.cells,i)
  
  for(j in other){
    
    move.mx[i,j] = round(runif(1,min=0,max=0.1),2)
    
  }
  move.mx[i,i] = 1-sum(move.mx[i,])
} 
move.mx

fec.mx = matrix(data=0,nrow=9,ncol=9);fec.mx[1,3] = fec.mx[4,6] = fec.mx[7,9] = 100

move.mx%*%fec.mx

rm(list=ls())
pop.mx = matrix(data = c(0.1, 0,   100,   0,0,0,
                         0.1, 0.5, 0,     0,0,0,
                         0,   0.1, 0.8,   0,0,0,
                         
                         0,0,0,    0.1, 0,   100,
                         0,0,0,    0.1, 0.5, 0,
                         0,0,0,    0,   0.1, 0.8),
                byrow=T,ncol=6)

move.mx = matrix(data=0,nrow = 6,ncol = 6)
move = 0.25
move.mx[4,1]=move
move.mx[1,4]=move
diag(move.mx)=1
move.mx[1,1]=1-move
move.mx[4,4]=1-move
move.mx%*%pop.mx
