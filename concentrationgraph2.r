library(Matrix)
library(MASS)
library(spam)
library(igraph)
library(microbenchmark)
library(ggplot2)
library(Rcpp)
sourceCpp("fillgraph.cpp")

ipf <- function(Y,amat,tol){   # Hastie Friedman Tibshirani p. 791
    cat('...IPF')
    p = dim(Y)[2]
    n = dim(Y)[1]
    S = (1/n)*t(Y)%*%Y
    k = ncol(S)
    W0 = W = S
    it = 0
    converge = FALSE
    while( !converge ) {
        it = it+1
        for (j in 1:k){   
            W11 = W[-j,-j,drop=FALSE]     
            w12 = W[-j,j]     
            s12 = S[-j,j, drop=FALSE]
            paj = amat[j,] == 1; # neighbors
            paj = paj[-j]
	         beta = rep(0, k-1)
            if (all(!paj)){
                w = rep(0, k-1)  
            }
            else{
                beta[paj] <- solve(W11[paj, paj], s12[paj, ])
                w = W11 %*% beta
            }
            W[-j, j] = w
            W[j, -j] = w
        }
        di <- norm(W0-W)      
        if (di < tol){
            converge = TRUE
        }
        else {
            W0 <- W 
          }
    }   
    return(list(Shat = W, it=it))
}

concentrationgraph <- function(Y,G){
                                        #Step1 of Concgraph.calculate sample covariance S
    ###cat('Calculating Covariance Matrix')
    n = dim(Y)[1]
    p = dim(Y)[2]
    S = cov(Y)


                                        #using the reverse cuthill mckee ordering algorithm
    ###cat('...Reordering columns and rows')
    out = chol(as.spam(S*(G) + p*diag(x = max(abs(S)), p)), pivot ="RCM")
    Q = ordering(out)
    S = S[Q,Q]
    G = G[Q,Q]
                                        #Step6: Fill the graph - this is taking a lot of time. perhaps this can be changed
    ###cat('...Filling the graph') 
    Gfill = gNM1(G,p)
                                        #step7: calculate lower triangular matrix L
    ###cat('...Calculating Cholesky for filled graph')
    L = diag(p)
    sqD = rep(0,p)
    for (i in 1:(p-1)){
        Snew = S[i:p,i:p]
        newp = dim(Snew)[1]
        temp = which(Gfill[(i+1):p,i]>0)
        if(length(temp)>0){
            Stempi = S[(i+1):p,i]
            Sdoti = Stempi[temp]
            Stempgp = S[(i+1):p,(i+1):p]
            if(length(Stempgp)>1){
                (Sgp = Stempgp[temp,temp])
            } else {
                Sgp = Stempgp
            }
            L[(temp+i),i] = -solve(Sgp,Sdoti)
            sqD[i] = sqrt(1/(Snew[1,1] + sum(t(Sdoti)*L[(temp+i),i])))
        } else {
            sqD[i] = sqrt(1/Snew[1,1])
            L[(temp+i),i] = 0
        } 
    }
    sqD[i] = sqrt(1/S[p,p])

    L = L%*%diag(sqD)
                                        #Step8: implement algorithm 1.
                                        #(L = algo1(Ld,adjacencyMatrix))
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


p = 1000  # predictors
n = 500  ############# this needs to change to 250
z = 0.005  # fraction of non-zeros in cholesky paramter of L
##### for p = 1000, z = 0.007 returns an omega with sparsity of 4.8%
##### for p = 1000, z = 0.005 returns an omega with sparsity of 2.8%
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
omega = t(L) %*% L            # omega
sigma = solve(omega)
sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero


sum(abs(omega)>0)/choose(p,2) #sparsity in omega

set.seed(23456 + 1) #seed for generating data
X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations

X = scale(X, center = TRUE, scale = TRUE) # centered obs

#g <- graph_from_adjacency_matrix(abs(omega)>0)
#cliques(g)
algo1time = proc.time()
algo1 = concentrationgraph(X,abs(omega)>0)$Shat
algo1time = proc.time() - algo1time
algo1time
norm(algo1-omega,type="F")/norm(omega,type="F")


ipftime = proc.time()
ipf = ipf(X,abs(omega)>0,10^(-5))$Shat
ipftime = proc.time() - ipftime
ipftime
norm(ipf-omega,type="F")/norm(omega,type="F")


