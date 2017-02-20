ipf <- function(Y,amat){   # Hastie Friedman Tibshirani p. 791
    p = dim(Y)[2]
    n = dim(Y)[1]
    S <- (1/n)*t(Y)%*%Y
    W0 <- S ; W <- S
    it <- 0
    converge = FALSE
    while( !converge ) {
        it <- it+1
        for (j in 1:k){   
            W11 <- W[-j,-j,drop=FALSE]     
            w12 <- W[-j,j]     
            s12 <- S[-j,j, drop=FALSE]
            paj <- amat[j,] == 1; # neighbors
            paj <- paj[-j]
	         beta <- rep(0, k-1)
            if (all(!paj)){
                w <- rep(0, k-1)  
            }
            else{
                beta[paj] <- solve(W11[paj, paj], s12[paj, ])
                w <- W11 %*% beta
            }
            W[-j, j] <- w
            W[j, -j] <- w
        }
        di <- norm(W0-W)      
        if(pri) {
            cat(di, "\n")
        }
        if (di < tol){
            converge <- TRUE
        }
        else {
            W0 <- W 
          }
    }   
    return(list(Shat = W, it=it))
}
