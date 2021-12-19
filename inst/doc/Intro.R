## ----eval=TRUE----------------------------------------------------------------
sirs<-function(X, y, a = 0.1) {
  n <- length(X[ ,1])
  p <- length(X[1, ])
  b0 <- rep(1, p)
  S <- min(floor(n/2), p)
  eps <- 1/p
  thresh <- 1e-6
  L <- 50
  M = min(L, S)
  tol = 1e-6
  k = 1
  diff = 1
  b.old = numeric(p)
  
  sirscore <- function(X, y, a = 0.1, b0 = rep(1, length(X[1, ]))){
    delta <- 1e-6
    dj <- function(t) abs(t) * (a + abs(t))/(a+1)
    D <- diag(sapply(b0, dj), p, p)
    b <- b0
    diff <- 1
    l <- 1
    dj <- function(t) abs(t) * (a + abs(t))/(a+1)
    D <- diag(sapply(b0, dj), p, p)
    while (diff > tol & l <= L){
      l = l + 1
      bold = b
      b = D%*%t(X)%*%(solve(diag(delta, n, n) + X%*%D%*%t(X)))%*%y
      diff = sqrt(sum((b - bold)^2))
      D = diag(sapply(b, dj), p, p)
    }
    b = b * (abs(b) > thresh)
    
    if (sum(b != 0) <= S){
      b_non0 = which(b != 0)
      X_non0 = X[, b_non0]
      b[b_non0] = solve(t(X_non0)%*%X_non0)%*%t(X_non0)%*%y
      b = b * (abs(b) > thresh)
    }
    return (b)
  }
  
  while (diff > tol & k <= M){
    b = sirscore(X, y, a, b0)
    l0 = sum(b != 0)
    if (l0 <= S){
      break
    }
    else {
      b_non0 = which(b != 0)
      index = order(abs(b[b_non0]), decreasing = TRUE)
      
      b_ini = rep(eps, p)
      b_ini[b_non0[index[1:k]]] = 1
      
      diff = sqrt(sum((b - b.old)^2))
      b.old = b
      k = k + 1
    }
  }
  return (b)
}

## ---- eval = TRUE-------------------------------------------------------------
cv.sirs <- function(X, y, k = 10, a = seq(0,1,0.1)){
  n <- length(X[, 1])
  cvgroup <- function(k, datasize){
    cvlist <- list()
    n <- rep(1:k, ceiling(datasize/k))[1:datasize] 
    temp <- sample(n, datasize)
    x <- 1:k
    dataseq <- 1:datasize
    cvlist <- sapply(x, function(x) {dataseq[temp==x]})
    return(cvlist)
  }
  
  loss_rho <- function(a, beta){
    rho_a <- function(t) abs(t) * (a + abs(t))/(a+1)
    loss <- rho_a(beta[beta!=0])
    return (sum(loss))
  }
  
  loss <- numeric(length(a))
  cvlist <- cvgroup(k, n)
  
  loss <-sapply((1:length(a)), FUN = function(i){
    tmp<-sapply(1:k, function(j){
      Xj = X[-cvlist[[j]],]
      yj = y[-cvlist[[j]]]
      beta = sirs(Xj, yj, a[i])
      return (beta)
    })
    return (sum(apply(tmp, MARGIN = 2, FUN = function(beta){loss_rho(a[i], beta)})))
  })
  a.best <- a[which.min(loss)]
  return (a.best)
}

## ----echo = TRUE--------------------------------------------------------------
library(MASS)
set.seed(12345)
n <- 35
p <- 500
beta_0 <- numeric(p)
beta_0[1:7] <- c(1, -0.5, 0.7,-1.2, -0.9, 0.3, 0.55)
r <- 0.1
gamma_r <- matrix(r, p, p) + diag(1 - r, p, p)
X <- mvrnorm(n, numeric(p), gamma_r)
y <- X%*%beta_0
k <- 2
a <- c(0.2, 0.3)
a.best <- (cv.sirs(X, y, k, a))
print(a.best)
beta_hat <- sirs(X, y, a.best)
print(sum(beta_hat != 0))
print(which(beta_hat != 0))

