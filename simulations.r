library(Matrix)
library(MASS)
library(spam)
library(igraph)
library(microbenchmark)
library(ggplot2)
library(Rcpp)
sourceCpp("fillgraph.cpp")
source('concentrationgraph.r')


name = "cscs"

p = AA   # predictors
n = BB  ############# this needs to change to BB 
z = 0.0065  # fraction of non-zeros
##### for p = 500, z = 0.0065 gives 2.8% density
s = 0.5  # fraction of negative coefficients
a = 0.3  # minimum magnitude of non-zero coefficients
b = 0.7  # maximum magnitude of non-zero coefficients

plower = p*(p-1)/2

set.seed(12345) #seed for generating L
## diagonals
D = runif(p,2,5)

## off-diagonals
T = diag(p)
T[upper.tri(T)] = 0
T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) * 
		   ifelse(runif(plower)<z,  1, 0) * 
		   runif(plower, a, b))

L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
omega = t(L) %*% L            # omega
sigma = solve(omega)
sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero

sum(abs(omega)>0)/choose(p,2) #sparsity in omega

set.seed(23456 + CC) #seed for generating data
X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations

X = scale(X, center = TRUE, scale = FALSE) # centered obs

algo1time = proc.time()[3]
algo1 = concentrationgraph(X,abs(omega)>0)$Shat
algo1time = proc.time()[3] - algo1time
algo1norm = norm(algo1-omega,type="F")/norm(omega,type="F")


ipftime = proc.time()[3]
ipfshat = ipf(X,abs(omega)>0,10^(-5))$Shat
ipftime = proc.time()[3] - ipftime
ipfnorm = norm(ipfshat-omega,type="F")/norm(omega,type="F")


write.table(data.frame(algo1time = algo1time,algo1norm = algo1norm, ipftime = ipftime,ipfnorm = ipfnorm),file=paste("Concentrationp",toString(p),"n",toString(n),toString(BB),".txt",sep = ""),row.names = TRUE)



