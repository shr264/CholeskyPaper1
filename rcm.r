#### reverse cuthill mckee ordering

install.packages("spam")
library(spam)

G = matrix(1,nrow=4,ncol=4)
G[3,1] = G[1,3] = G[4,2] = G[2,4] = 0
G
Q =  R = rep(0,0) ### prepare empty queue
colsums = apply(G,2,sum)
P = which.min(colsums)
R = append(R,P)
Q = append(Q,which(G[,R]>0))
nodesleft = setdiff(1:4,Q[1:colsums[P]])
while(length(nodesleft)>0){
    newG = G[nodesleft,nodesleft]
    if(length(newG)[1]==1){
        Q = append(Q,nodesleft)
        break
    }
    colsums = apply(newG,2,sum)
    P = which.min(colsums)
    R = append(R,P)
    Q = append(Q,which(G[,R]>0))
    nodesleft = setdiff(1:4,Q[1:colsums[P]])
}
G[Q,Q]
invPerm(Q)

out <- chol(as.spam(G), pivot ="RCM")
ordering(out)
