library(Matrix)
library(MASS)
library(spam)
library(igraph)
library(microbenchmark)
library(ggplot2)
library(Rcpp)
library(ggm)
sourceCpp("fillgraph.cpp")



covariancegraph <- function(Y,G){
                                        #Step1 of Concgraph.calculate sample covariance S
    ###cat('Calculating Covariance Matrix')
    n = dim(Y)[1]
    p = dim(Y)[2]
    S = cov(Y)

    out = chol(as.spam(S*(G) + p*diag(x = max(abs(S)), p)), pivot ="RCM")
    Q = ordering(out)
    S = S[Q,Q]
    G = G[Q,Q]

    Gfill = gNM1(G,p)
    
    L = t(chol(S))
                                        #Step8: implement algorithm 1.                                    #(L = algo1(Ld,adjacencyMatrix))
    (E = matrix(0,0,2))
    for(i in 2:p){
        (zero = which(Gfill[i,1:i]==0))
        if(length(zero)>0){
            (flintrack = cbind(rep(i,length(zero)),zero))
            (E = rbind(E,flintrack))}
    }

    (a <- length(E[,1]))


    cat('...Algorithm 1')
    for(i in 1:a){
        if(E[i,2]>1){
            (L[E[i,1],E[i,2]] = -sum(L[E[i,1],1:(E[i,2]-1)]*L[E[i,2],1:(E[i,2]-1)]))
            (L[E[i,1],E[i,2]] = (L[E[i,1],E[i,2]])/(L[E[i,2],E[i,2]]))
        }
        else{       
            L[E[i,1],E[i,2]] = 0}
    }
    
        
    (omegahat = L%*%t(L))
    omegahat = omegahat[invPerm(Q),invPerm(Q)]
    G = G[invPerm(Q),invPerm(Q)]
    Gfill = Gfill[invPerm(Q),invPerm(Q)]
    return(list(Shat=omegahat))
}


p = 1000 # predictors
n = 3000 ############# this needs to change to 250
z = 0.006  # fraction of non-zeros in cholesky paramter of L
##### for p = 1000, z = 0.007 returns an omega with sparsity of 4.8%
##### for p = 1000, z = 0.005 returns an omega with sparsity of 2.8%
##### for p = 1000, z = 0.004 returns an omega with sparsity of 2.0%
##### for p = 2000, z= 0.005 returns an omega with sparsity of 4.3%
##### for p = 2000, z= 0.004 returns an omega with sparsity of 2.97%
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
sigma = t(L) %*% L            # omega
omega = solve(sigma)
omega[abs(omega)<1e-10] = 0   # set numerical error to zero


sum(abs(sigma)>0)/choose(p,2) #sparsity in omega

set.seed(23456 + 1) #seed for generating data
X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations

X = scale(X, center = TRUE, scale = TRUE) # centered obs

#g <- graph_from_adjacency_matrix(abs(omega)>0)
#cliques(g)
algo1time = proc.time()
algo1 = covariancegraph(X,abs(sigma)>0)$Shat
algo1time = proc.time() - algo1time
algo1time
norm(algo1-sigma,type="F")/norm(sigma,type="F")


icftime = proc.time()
S = (1/n)*t(X)%*%X
amat = (abs(sigma)>0) - diag(p)
dimnames(amat)[[1]] = as.list(seq(1:p))
dimnames(amat)[[2]] = as.list(seq(1:p))
dimnames(S)[[1]] = as.list(seq(1:p))
dimnames(S)[[2]] = as.list(seq(1:p))
icf = fitCovGraph(amat,S = S, n = n)$Shat
icftime = proc.time() - icftime
icftime
norm(icf-sigma,type="F")/norm(sigma,type="F")


