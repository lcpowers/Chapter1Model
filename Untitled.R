rm(list=ls())

pop = matrix(byrow=T,nrow=3,
           data=c(1,0,100,
                  1,1,0,
                  0,1,1))


move = matrix(byrow=T,nrow=3,
           data=c(0.5,0,0,
                  0,1,0,
                  0,0,1))

move%*%pop

pop = matrix(byrow=T,ncol=3,
             data=c(1,0,100,
                    1,1,0,
                    0,1,1))

move = matrix(byrow=T,nrow=3,
              data=c(1,1,0.5,
                     1,1,1,
                     1,1,1))

pop*move

# 1 matrix col == 2 matrix rows

pop%*%move
# number of columns of the first matrix must be the same as the number of rows of the second matrix. If the first matrix has dimensions m × n, and is multiplied by a second matrix of dimensions n × p, then the dimensions of the product matrix will be m × p.


### Kronecker product
move = matrix(data = c(0.1,0.25),nrow=1)
move

pop = matrix(byrow=T,ncol=3,
             data=c(1,0,100,
                    1,1,0,
                    0,1,1))


big.pop = matrix(data=0,nrow=6,ncol=6)

big.pop[1:3,1:3] <- big.pop[4:6,4:6] <- pop
big.pop

kronecker(Y=big.pop,X=move)


