cmpGraph <- function(amat){
    g <- 1*!amat
    diag(g) <- 0
    g
  }

icf <- function(Y,amat,tol, start=NULL){
    p = dim(Y)[2]
    n = dim(Y)[1]
    S = cov(Y)

    pa.each.node <-function(amat){
        ## List of the parents of each node.
        ## If amat is symmetric it returns the boundaries.
        p <- nrow(amat)
      b <- vector(p, mode="list")
      ip <- 1:p
      for(i in 1:p)
          b[[i]] <- ip[amat[,i]==1]
      b 
    }

    spo <- pa.each.node(amat)
    nsp <- pa.each.node(cmpGraph(amat)) #### what is cmpGraph
    number.spouses <- unlist(lapply(spo, length))
    no.spouses <- (1:p)[number.spouses==0]
    all.spouses <- (1:p)[number.spouses==(p-1)]
    nontrivial.vertices <- setdiff((1:p), no.spouses)

    if(length(nontrivial.vertices)==0){
      if(p==1){
        Sigma <- S
      }
      else{
        Sigma <- diag(diag(S))
        dimnames(Sigma) <- dimnames(S)
      }
      return(list(Sigmahat=Sigma, iterations=1))
    }

    if(is.null(start)){
      Sigma <- as.matrix(diag(diag(S))) # starting value
    }
    else{
      temp <- diag(start)
      start[amat==0] <- 0
      diag(start) <- temp
      diag(start)[no.spouses] <- diag(S)[no.spouses]
      Sigma <- as.matrix(start)
      if(min(eigen(Sigma)$values) <= 0){
        stop("Starting value is not feasible!")
      }
    }
    
    repeat{
      i <- i+1
      Sigma.old <- Sigma
      for(v in nontrivial.vertices){
        if(is.element(v, all.spouses)){
          B <- S[v,-v]%*%solve(S[-v,-v])
          lambda <- S[v,v]-B%*%S[-v,v]
        }
        else{
          B.spo.nsp <-
            Sigma[spo[[v]],nsp[[v]]]%*%solve(Sigma[nsp[[v]],nsp[[v]]])
          YZ <- S[v,spo[[v]]]-S[v,nsp[[v]]]%*%t(B.spo.nsp) 
          B.spo <- 
            YZ %*%
              solve( S[spo[[v]],spo[[v]]]
                    -S[spo[[v]],nsp[[v]]]%*%t(B.spo.nsp)
                    -B.spo.nsp%*%S[nsp[[v]],spo[[v]]]
                    +B.spo.nsp%*%S[nsp[[v]],nsp[[v]]]%*%t(B.spo.nsp) )
          lambda <- S[v,v]-B.spo%*%t(YZ)
          B.nsp <- -B.spo%*%B.spo.nsp
          B <- rep(0, p)
          B[spo[[v]]] <- B.spo
          B[nsp[[v]]] <- B.nsp
          B <- B[-v]
        }
        Sigma[v,-v] <- B%*%Sigma[-v,-v]
        Sigma[v,nsp[[v]]] <- 0
        Sigma[-v,v] <- t(Sigma[v,-v])
        Sigma[v,v] <- lambda + B%*%Sigma[-v,v]
      }
      if(sum(abs(Sigma.old-Sigma)) < tol) break
    }
    dimnames(Sigma) <- dimnames(S)
    return(list(Sigmahat=Sigma, iterations=i))
  }
