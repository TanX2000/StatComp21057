## ----echo = TRUE--------------------------------------------------------------
x <- matrix(seq(1, 9, 1), 3, 3)
knitr::kable(x)

## ----echo = TRUE--------------------------------------------------------------
x <- seq(-1, 1, 0.1)
y <- x^2
plot(x, y)

## ----echo = TRUE--------------------------------------------------------------
ymax <- max(y)

## ----echo = TRUE--------------------------------------------------------------
n<-1e3 #每次生成1000个数据
sigma<-c(0.2,1,5,20,50,100)
par(mfrow=c(2,3)) 
for (i in 1:length(sigma)){
  set.seed(1234)
  u<-runif(n) #生成U~U(0,1)
  x<-sigma[i]*sqrt(-2*log(1-u)) #按sigma从小到大的顺序依次利用F(x)的逆生成X
  upper<-floor(max(x))+1 #找到生成的X的上界，为方便f(x)实际密度函数作图
  t<-seq(from=0,to=upper,by=upper/n) 
  hist(x,prob=TRUE,breaks=15)
  lines(t,t*exp(-t*t/2/sigma[i]/sigma[i])/sigma[i]/sigma[i],col='red',lwd=2)#在直方图中作出f(x)实际的密度函数曲线作为对比
}

## ----echo = TRUE--------------------------------------------------------------
set.seed(1234)
n <- 1000
X1 <- rnorm(n,0,1)
X2 <- rnorm(n,3,1)
r<-sample(c(0,1),n,replace=TRUE,prob=c(0.25,0.75)) #按0.25与0.75的概率比生成0-1分布
Z <- r*X1+(1-r)*X2 #r有0.75的可能取1，于是Z有0.75的可能取X1
hist(Z,prob=TRUE)

## ----echo = TRUE--------------------------------------------------------------
p=seq(0.3,0.7,0.05) #从0.3至0.7每间隔0.05作一次p的取值
#par(mfrow=c(3,3)) #生成3x3直方图,从左至右、从上到下p值依次增大0.5
for (i in 1:length(p)){
  set.seed(1234)
  r<-sample(c(0,1),n,replace=TRUE,prob=c(p[i],1-p[i]))
  Z <- r*X1+(1-r)*X2
  hist(Z,prob=TRUE)
}

## ---- echo = TRUE-------------------------------------------------------------
set.seed(1234)
lambda <- c(1,1,5,5)
a<-c(1,2,1,2)
b<-c(1,5,1,5)
t0 <- 10
estimate_mean<-numeric(4)
estimate_var<-numeric(4)
theory_mean<-numeric(4)
theory_var<-numeric(4)
for (i in 1:4){
  set.seed(1234)
  X<-replicate(10000, expr = {
    T <- rexp(200, lambda[i])  #interarrival times
    S <- cumsum(T)  #arrival times：cumsum计算累积和S[n]=T[1]+T[2]+...+T[n]
    n <- min(which(S > t0)) #arrivals+1 in [0, t0]： N(t0)为落在[0,t0]区间的次数，故先找到使得S[n]最先超过t0的n，将其-1即为所求
    Y<-rgamma(n-1,a[i],b[i]) #生成N(t)个服从Gamma(a,b)的Y
    sum(Y) #X(t)的取值为这N(t)个Y的总和
    })
  estimate_mean[i]<-mean(X) #计算均值的估计值
  estimate_var[i]<-var(X) #计算方差的估计值
  theory_mean[i]<-lambda[i]*t0*a[i]/b[i] #计算均值的理论值
  theory_var[i]<-lambda[i]*t0*((a[i]+a[i]*a[i])/(b[i]*b[i])) #计算方差的理论值
}
parameters=c('(1,1,1)','(1,2,5)','(5,1,1)','(5,2,5)') #按lambda,a,b的顺序
table1<-cbind(parameters,estimate_mean,theory_mean)
knitr::kable(table1,align='l') #生成期望估计值与理论值的比较表格
table2<-cbind(parameters,estimate_var,theory_var)
knitr::kable(table2,align='l') #生成期望估计值与理论值的比较表格


## ----echo = TRUE--------------------------------------------------------------
n<-1e4
p<-seq(0.1,0.9,0.1)
theta.hat<-numeric(length(p))
theta<-numeric(length(p))
for (i in 1:length(p)){
  set.seed(504)
  u<-runif(n,0,p[i])
  theta.hat[i] <- round(p[i]*mean(u*u*(1-u)*(1-u)/beta(3,3)),4)
  theta[i]<-round(pbeta(p[i],3,3),4)
}
table5_4=cbind(p,theta,theta.hat)
knitr::kable(table5_4,align='c')

## ----echo = TRUE--------------------------------------------------------------
sigma=1
set.seed(123)
u1<-runif(1000) #生成U(0,1)上随机数
x1<-sigma*sqrt(-2*log(1-u1))
x1_new<-sigma*sqrt(-2*log(u1)) #用X1生成其对偶变量X1'
set.seed(1234)
u2<-runif(1000)
x2<-sigma*sqrt(-2*log(1-u2))  #生成新变量X2
c(mean((x1+x1_new)/2),mean((x1+x2)/2))
var_anti<-round(var((x1+x1_new)/2),4)
var_ordinary<-round(var((x1+x2)/2),4)
var_reduction<-round((var_ordinary-var_anti)/var_ordinary,4)
table5_9<-cbind(var_anti,var_ordinary,var_reduction)
knitr::kable(table5_9,align='c')

## ----echo = TRUE--------------------------------------------------------------
x <- seq(1,10, 0.01)
w <- 2
g <- x^2*exp(-x^2/2)/sqrt(2*pi)
f1 <- x*exp(-x^2/2)
f2 <- exp(-x^2/2)/sqrt(2*pi)

par(mfrow=c(1,1))
gs <- c(expression(g(x)==x^2*e^{-x^2/2}/sqrt(2*pi)),
  expression(f[1](x)==x*e^{-x^2/2}),
  expression(f[2](x)==e^{-x^2/2}/sqrt(2*pi)))
##生成g,f1,f2的图线
plot(x,g,type="l",ylim=c(0,1), lwd = w, lty = 2,col=1,main='(A)')
lines(x,f1,type="l",ylim=c(0,1), lwd = w, lty = 2,col=2,main='(A)')
lines(x,f2,type="l",ylim=c(0,1), lwd = w, lty = 2,col=3,main='(A)')
legend("topright", legend = gs,
       lty = 1:3, lwd = w, inset = 0.02,col=1:3)
##生成g/f1,g/f2的图线
plot(x, g/f1, type = "l", ylab = "",
     ylim = c(0,130), lwd = w, lty = 2,col=2,main='(B)')
lines(x, g/f2, lty = 3, lwd = w,col=3)
legend("topright", legend = gs[-1],
       lty = 1:3, lwd = w, inset = 0.02,col=2:3)

## ----echo = TRUE--------------------------------------------------------------
m<-1e4
g = function (x) {
  x ^ 2 / sqrt(2*pi) * exp(-x^2/2) *(x>1) #只取>1部分
}
est<-numeric(2)
v<-numeric(2)
##f1计算积分
set.seed(1)
f1 = function (x) {
  x*exp(-x^2/2)
}
u<-runif(m)
x<-sqrt(-2*log(1-u))#逆变换法生成服从Rayleigh分布的随机数X
i <- c(which(x<=1))
x[i] <- -1 #只取>1部分
fg<-g(x)/f1(x)
est[1]<-mean(fg)
v[1]<-var(fg)

##f2计算积分
set.seed(1)
x <- rnorm(1000) #using f2
i <- c(which(x<=1))
x[i] <- -1 #只取>1部分
fg <- g(x) / dnorm(x)
est[2]<-mean(fg)
v[2]<-var(fg)
res <- rbind(estimate=round(est,4), variance=round(v,4))
colnames(res) <- paste0('f',1:2)
knitr::kable(res, align='c')

## ----echo = TRUE--------------------------------------------------------------
n<-20
set.seed(1)
n <- 20
alpha <- .05
x_mean<-numeric(1e3)
x_sd<-numeric(1e3)
for (i in 1:1e3){
  x<-rchisq(n,df=2)
  x_mean[i]<-mean(x) #样本均值
  x_sd[i]<-sd(x) #样本标准差
}
t_alpha<-qt(p=1-alpha/2,df=n-1)
coverage<-mean(x_mean-t_alpha*x_sd/sqrt(n)<=2&x_mean+t_alpha*x_sd/sqrt(n)>=2) #所有置信区间中能够覆盖真值2的所占的比例
round(coverage,4)

##例6.4
n <- 20
set.seed(1)
alpha <- .05
UCL <- replicate(1000, expr = {
  x <- rchisq(n, df = 2)
  (n-1) * var(x) / qchisq(alpha, df = n-1)
} )
coverage2<-mean(UCL > 4)
round(coverage2,4)

## ----echo = TRUE--------------------------------------------------------------
set.seed(1234)
n <- 1e3
alpha <- 0.1
m <- 5e3 #重复次数
p_x<-p_y<-p_z <- numeric(m) #存储p值
for (j in 1:m) {
  x <- rchisq(n,1) #X^2(1)
  y <- runif(n,0,2) #U(0,2)
  z <- rexp(n,1) #Exp(1)
  ttest_x <- t.test(x,mu=1)
  ttest_y <- t.test(y,mu=1)
  ttest_z <- t.test(z,mu=1)
  p_x[j] <- ttest_x$p.value
  p_y[j] <- ttest_y$p.value
  p_z[j] <- ttest_z$p.value
}

p_chisq.hat <- mean(p_x <= alpha)
p_unif.hat <- mean(p_y <= alpha)
p_exp.hat <- mean(p_z <= alpha)
table1<-rbind(p_chisq.hat,p_unif.hat,p_exp.hat) #3种分布的t1e
knitr::kable(table1,align='l')

se_chisq.hat <- sqrt(p_chisq.hat * (1 - p_chisq.hat) / m)
se_unif.hat <- sqrt(p_unif.hat * (1 - p_unif.hat) / m)
se_exp.hat <- sqrt(p_exp.hat * (1 - p_exp.hat) / m)
table2<-rbind(se_chisq.hat,se_unif.hat,se_exp.hat) #3种分布t1e标准差的估计值
knitr::kable(table2,align='l')

## ----echo = TRUE--------------------------------------------------------------
library(MASS)
sk <- function(data) {
  n<-nrow(data) 
  d<-ncol(data)
  mu<-numeric(d)
  mat_cent<-diag(1,n)-matrix(1/n,n,n) #中心化X时的右乘矩阵,可以看作单位阵I-所有元素全为1/n的矩阵
  #cov_hat<-t(data)%*%mat_cent%*%data/n
  cov_hat<-cov(data)*(n-1)/n #r中cov函数计算的是修正后的协方差阵,除数是n-1
  data_cent<-mat_cent%*%data #中心化data
  mat_bd<-data_cent%*%solve(cov_hat)%*%t(data_cent) #solve函数求逆矩阵
  b1d<-sum(mat_bd^3)/n^2
  return (b1d)
}

set.seed(123)
n <- c(10,20,30,50,100,300) #每次的样本量
p.reject <- numeric(length(n)) 
m <- 1000 #每次模拟的重复次数
d<-2 #多元分布的维数
dg<-d*(d+1)*(d+2)/6 #X^2分布的自由度
cv <- qchisq(.95,dg)*6/n #检验临界值
mean<-c(0,0) 
sigma<-matrix(c(1,0,0,1),nrow=d,ncol=d)
for (i in 1:length(n)) {
  cnts_reject <- numeric(m) #test decisions
  for (j in 1:m) {
    x <- mvrnorm(n[i],mean,sigma)
    cnts_reject[j] <- as.integer(sk(x) >= cv[i])
  }
  p.reject[i] <- mean(cnts_reject) #拒绝比例
}
table<-rbind(n,p.reject)
knitr::kable(table)

## ----echo = TRUE--------------------------------------------------------------
set.seed(123)
alpha <- 0.05
n <- 100 #样本量
m <- 1000 #每次重复的次数
epsilon <- c(seq(0, 0.2, 0.02), seq(0.2, 1, 0.05)) 
N <- length(epsilon)
power <- numeric(N)
cv <- qchisq(.95,dg)*6/n #检验临界值
Sigma_1<-diag(1,d,d) #MN1的协方差阵
Sigma_10<-diag(10,d,d) #MN2的协方差阵
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  cnts_reject <- numeric(m)
  for (i in 1:m) { 
    sigma <- sample(c(1, 10), replace = TRUE,
                    size = n, prob = c(1-e, e)) #以概率1-e由N1生成,概率2由多元正态分布2生成
    n_1<-length(which(sigma==1)) #由MN1生成的数量
    n_10<-length(which(sigma==10)) #由MN2生成的数量
    if (n_1>0 & n_10>0){
      x1<-mvrnorm(n_1,mean,Sigma_1)
      x2<-mvrnorm(n_10,mean,Sigma_10)
      x<-rbind(x1,x2)
    }
    else if (n_1>0){
      x<-mvrnorm(n_1,mean,Sigma_1)
    }
    else {
      x<-mvrnorm(n_10,mean,Sigma_10)
    }
    cnts_reject[i] <- as.integer(sk(x) >= cv)
  }
  power[j] <- mean(cnts_reject)
}
plot(epsilon, power, type = "p",
     xlab = bquote(epsilon), ylim = c(0,1),main="empirical power curve")
lines(epsilon,power,lwd=2,col="red")
abline(h =0.05, lty = 3,lwd=2) #直线:y=0.05
se <- sqrt(power * (1-power) / m) #标准差
lines(epsilon, power+se, lty = 3) 
lines(epsilon, power-se, lty = 3)
abline(v=0.16,lty=3,lwd=2) #最高点在0.16附近取得
table2<-rbind(epsilon,power)
knitr::kable(table2)

## ----echo = TRUE--------------------------------------------------------------
library(bootstrap)
attach(scor)
lambda_hat<-eigen(cov(scor))$values #eigen求特征值
theta_hat<-max(lambda_hat)/sum(lambda_hat)
n<-nrow(scor)
B<-1000
set.seed(1234)
theta_B<-numeric(B)
for (i in 1:B){
  ind<-sample(1:n,size=n,replace=TRUE)
  scor_i<-scor[ind,]
  lambda_i<-eigen(cov(scor_i))$values
  theta_B[i]<-max(lambda_i)/sum(lambda_i)
}
theta_B_hat<-mean(theta_B)
bias_B<-theta_B_hat-theta_hat
se_B<-sd(theta_B)
knitr::kable(cbind(theta_B_hat,bias_B,se_B))

## ----echo = TRUE--------------------------------------------------------------
theta_j<-numeric(n)
for (i in 1:n){
  scor_i<-scor[-i,]
  lambda_i<-eigen(cov(scor_i))$values
  theta_j[i]<-max(lambda_i)/sum(lambda_i)
}
theta_j_hat<-mean(theta_j)
bias_j<-(n-1)*(theta_j_hat-theta_hat)
se_j<-sqrt((n-1)/n*sum((theta_j-theta_j_hat)^2))
knitr::kable(cbind(theta_j_hat,bias_j,se_j),algn='c')

## ----echo = TRUE--------------------------------------------------------------
library(boot)
set.seed(1234)
boot.theta<- function(x,i){
  lambda_i<-eigen(cov(x[i,]))$values
  theta_B<-max(lambda_i)/sum(lambda_i)
  theta_B
}
boot.obj <- boot(scor, statistic = boot.theta, R=500)
ci<-boot.ci(boot.obj, type=c("perc","bca"))
ci.perc<-ci$percent[4:5] ##percentile CI
ci.bca<-ci$bca[4:5] ##Bca CI
knitr::kable(rbind(ci.perc,ci.bca),align='c')

## ----echo = TRUE--------------------------------------------------------------
library(boot)
m<-1e3
n<-20
sk <- function(x,i) {
  #computes the sample skewness coeff.
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}
sk_chi_5 <- sqrt(8/5)
set.seed(0)
ci1.norm<-ci1.basic<-ci1.perc<-ci2.norm<-ci2.basic<-ci2.perc<-matrix(0,m,2)
for(i in 1:m){
  x1<-rnorm(n);x2<-rchisq(n,5)
  de1 <- boot(data=x1,statistic=sk, R = 1e3);de2<-boot(data=x2,statistic=sk,R=2e2)
  ci1 <- boot.ci(de1,type=c("norm","basic","perc"));ci2<-boot.ci(de2,type=c("norm","basic","perc"))
  ci1.norm[i,]<-ci1$norm[2:3];ci2.norm[i,]<-ci2$norm[2:3]
  ci1.basic[i,]<-ci1$basic[4:5];ci2.basic[i,]<-ci2$basic[4:5]
  ci1.perc[i,]<-ci1$percent[4:5];ci2.perc[i,]<-ci2$percent[4:5]
}
#coverage
normal_norm<-mean(ci1.norm[,1]<=0 & ci1.norm[,2]>=0);chisq_norm<-mean(ci2.norm[,1]<=sk_chi_5 & ci2.norm[,2]>= sk_chi_5)
normal_basic<-mean(ci1.basic[,1]<=0 & ci1.basic[,2]>=0);chisq_basic<-mean(ci2.basic[,1]<= sk_chi_5 & ci2.basic[,2]>=sk_chi_5)
normal_perc<-mean(ci1.perc[,1]<=0 & ci1.perc[,2]>=0);chisq_perc<-mean(ci2.perc[,1]<= sk_chi_5 & ci2.perc[,2]>= sk_chi_5)
#the proportion of times that the confidence intervals miss on the left
normal_norm_left<-mean(ci1.norm[,1]>0);chisq_norm_left<-mean(ci2.norm[,1]> sk_chi_5)
normal_basic_left<-mean(ci1.basic[,1]>0);chisq_basic_left<-mean(ci2.basic[,1]> sk_chi_5)
normal_perc_left<-mean(ci1.perc[,1]>0);chisq_perc_left<-mean(ci2.perc[,1]> sk_chi_5)
#the proportion of times that the confidence intervals miss on the right
normal_norm_right<-mean(ci1.norm[,2]<0);chisq_norm_right<-mean(ci2.norm[,2]< sk_chi_5)
normal_basic_right<-mean(ci1.basic[,2]<0);chisq_basic_right<-mean(ci2.basic[,2]< sk_chi_5)
normal_perc_right<-mean(ci1.perc[,2]<0);chisq_perc_right<-mean(ci2.perc[,2]< sk_chi_5)
ci_type <- c("norm", "basic", "perc")
#standard normal
normal_cover <- c(normal_norm, normal_basic, normal_perc)
normal_left<-c(normal_norm_left,normal_basic_left,normal_perc_left)
normal_right<-c(normal_norm_right,normal_basic_right,normal_perc_right)
table1 <- data.frame(ci_type,normal_cover,normal_left,normal_right)
knitr::kable(table1)
#chisq
chisq_cover <- c(chisq_norm, chisq_basic, chisq_perc)
chisq_left<-c(chisq_norm_left,chisq_basic_left,chisq_perc_left)
chisq_right<-c(chisq_norm_right,chisq_basic_right,chisq_perc_right)
table2 <- data.frame(ci_type,chisq_cover,chisq_left,chisq_right)
knitr::kable(table2)

## ----echo = TRUE--------------------------------------------------------------
set.seed(123)
x <- rnorm(20)
y <- rnorm(20)
R <- 399;z <- c(x, y);K <- 1:length(z);n<-length(x)
reps <- numeric(R);
t0 <- cor(x, y,method="spearman")
for (i in 1:R) {
  ind <- sample(K, size = n) # 置换后x对应的下标
  x1 <- z[ind];y1 <- z[-ind]
  reps[i] <- cor(x1, y1,method="spearman")
}
p <- round(mean(abs(c(t0, reps)) >= t0),3)
p0 <- round(cor.test(x,y,method="spearm")$p.value,3)
round(c(p,p0),3)

## ----echo = TRUE--------------------------------------------------------------
library("boot")
library("RANN")
library("energy")
library("Ball")

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

m <- 3e2; k<-3; p<-2; mu <- 0.3; set.seed(123)
n1 <- n2 <- 50; R<-99; n <- n1+n2; N = c(n1,n2)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1),ncol=p);
  y <- matrix(rnorm(n2*p,0,1.4),ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
}
alpha <- 0.1;
name<- c('knn','energy','ball')
power <- colMeans(p.values<alpha)
knitr::kable(rbind(name,power),align="c")

## ----echo = TRUE--------------------------------------------------------------
set.seed(123)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,0,1),ncol=p);
  y <- matrix(rnorm(n2*p,0.3,1.2),ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
}
alpha <- 0.1;
power <- colMeans(p.values<alpha)
knitr::kable(rbind(name,power),align="c")

## ----echo = TRUE--------------------------------------------------------------
set.seed(1234)
for(i in 1:m){
  x<-matrix(rt(n1*p,df=1),ncol=p)
  y <-rnorm(n2 * p,mean=c(0,0),sd = sample(c(1,1.3), prob = c(0.5,0.5),replace = TRUE))
  y<- matrix(y,ncol=p)
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=99,seed=i)$p.value
}
alpha <- 0.1;
power <- colMeans(p.values<alpha)
knitr::kable(rbind(name,power),align="c")

## ----echo = TRUE--------------------------------------------------------------
n1 <- 50;n2 <- 5;n <- n1+n2;N = c(n1,n2)
for(i in 1:m) {
  x <- matrix(rnorm(n1*p,0,1), ncol = p)#50 samples
  y <- matrix(rnorm(n2*p,0.9,1), ncol = p)#5 samples
  z <- rbind(x, y)
  p.values[i, 1] <- eqdist.nn(z, N, k)$p.value
  p.values[i, 2] <- eqdist.etest(z, sizes = N, R = R)$p.value
  p.values[i, 3] <- bd.test(x = x, y = y, R = 99, seed = i)$p.value
}
alpha <- 0.1
power <- colMeans(p.values<alpha)
knitr::kable(rbind(name,power),align="c")

## ----echo = TRUE--------------------------------------------------------------
cauchy.metropolis<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for (i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if (u[i]<=dcauchy(y)/dcauchy(x[i-1]))
      x[i]<-y
    else {
      x[i]<-x[i-1]
      k<-k+1
    }
  }
  x
}
set.seed(1234)
N<-5000 #链长为5000
sigma<-2.5 #sigma值选为2.5
x0<-0 #初值选为0
burn<-1e3 #预烧期为1000
chain<-cauchy.metropolis(sigma,x0,N)
#plot(chain,type="l")
a<-ppoints(25)
QR<-qcauchy(a)
Q<-quantile(chain[(burn+1):N],a)
qqplot(QR,Q,xlab="theorical",ylab="practical")
abline(c(0,0),c(1,1),lwd=2)
hist(chain[(burn+1):N],breaks="scott",freq=FALSE,ylim=c(0,0.1+1/pi),main="Hist",xlab=expression(X[t]))
lines(seq(-20,20,0.01),dcauchy(seq(-20,20,0.01)),lwd=2,col="red")

## ----echo = TRUE--------------------------------------------------------------
gibbs<-function(n,a,b,x0,N=5000){
  X <- matrix(0, N, 2)
  X[1, ] <- x0
  for (i in 2:N){
    x2 <- X[i-1,2]
    X[i,1]<-rbinom(1,n,x2)
    x1<-X[i,1]
    X[i,2]<-rbeta(1,x1+a,n-x1+b)
  }
  X
}
x0<-c(10,0.5) #x,y初值分别选为10,0.5
a<-b<-5;n<-20 #a,b,n,分别选为5,5,20
X<-gibbs(n,a,b,x0,N)
burn<-1000
x<-X[(burn+1):N,]
plot(x,xlab="x",ylab="y")
#knitr::kable(t(colMeans(x))) #x,y均值

## ----echo = TRUE--------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

## 9.3
set.seed(1234)
k <- 4 #生成链的条数
N <- 5000 #每条链的长度
burn <- 2000 #预烧期
x0 <- c(-5, -2, 2, 5) #4个不同初值
X <- matrix(0, nrow=k, ncol=N)
sigma<-2.5
for (i in 1:k)
  X[i, ] <- cauchy.metropolis(sigma, x0[i],N)
psi <- t(apply(X, 1, cumsum)) #cumsum计算累加和
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi)) #psi储存累积均值; 1:ncol(psi)对应每个时刻的链长
print(Gelman.Rubin(psi))
# par(mfrow=c(2,2))
# for (i in 1:k)
#   plot(psi[i, (b+1):N], type="l",xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1))
r.hat <- rep(0, (N-burn))
for (j in (burn+1):N)
  r.hat[j-burn] <- Gelman.Rubin(psi[,1:j])#生成Gelman-Rubin统计量r.hat
plot(r.hat, type="l", xlab="t", ylab=expression(hat(R)),ylim=c(1,max(r.hat)+0.1))
abline(h=1.2, lty=2) #临界值1.2

##9.8 
set.seed(123)
k <- 4 #生成链的条数
N <- 10000 #每条链的长度
burn <- 1000 #预烧期
X0<-matrix(c(7,0.6,10,0.55,10,0.45,13,0.4),ncol=2,nrow=4,byrow=TRUE)
X1 <- matrix(0, nrow=k, ncol=N)
X2 <- matrix(0, nrow=k, ncol=N)
for (i in 1:k){
  X1[i, ] <- gibbs(n,a,b,X0[i,],N)[,1]
  X2[i, ] <- gibbs(n,a,b,X0[i,],N)[,2]
}

psi1 <- t(apply(X1, 1, cumsum)) #cumsum计算累加和
psi2 <- t(apply(X2, 1, cumsum)) 
for (i in 1:nrow(psi1)){
  psi1[i,] <- psi1[i,] / (1:ncol(psi1)) #psi储存累积均值; 1:ncol(psi)对应每个时刻的链长
  psi2[i,] <- psi2[i,] / (1:ncol(psi2)) 
}
print(c(Gelman.Rubin(psi1),Gelman.Rubin(psi2)))
par(mfrow=c(1,1))
r1.hat <- r2.hat <- rep(0, (N-burn))
for (j in (burn+1):N){
  r1.hat[j-burn] <- Gelman.Rubin(psi1[,1:j])#生成Gelman-Rubin统计量r.hat
  r2.hat[j-burn] <- Gelman.Rubin(psi2[,1:j])
}
plot(r1.hat, type="l", xlab="t", ylab=expression(hat(R[1])),ylim=c(1,max(r1.hat)+0.1))
abline(h=1.2, lty=2) #临界值1.2
plot(r2.hat, type="l", xlab="t", ylab=expression(hat(R[2])),ylim=c(1,max(r2.hat)+0.1))
abline(h=1.2, lty=2) #临界值1.2

## ---- echo = TRUE-------------------------------------------------------------
###(1)
f_k<-function(a,k){
  d<-length(a)
  a_mod<-sqrt(as.numeric(a%*%a)) #计算a的2-范数
  f_k0<-a_mod^2*gamma((d+1)/2)*(gamma(3/2))/gamma(d/2+1)/2 #k=0时单独计算
  if (k==0){
    return (f_k0)
  }
  else {return ( (-1)^k * exp(log(gamma((d+1)/2)) + log(gamma(k+3/2)) - 
                                log(gamma(k+d/2+1)) - sum(log(seq(1,k,1))) - k*log(2) +
                                (2*k+2)*log(a_mod) - log(2*k+1) - log(2*k+2)))} #防止超界先取对数
}

###(2)
f_sum<-function(a){
  k<-0
  s<-0
  while (abs(f_k(a,k))>=1e-6){
    s<-s+f_k(a,k)
    k<-k+1
  } #只计算>10^{-6}的项
  return (s)
}
###(3)
a<-c(1,2)
print(round(f_sum(a),4))

## ----echo = TRUE--------------------------------------------------------------
# g <- function(u,k){
#   return ((1+u^2/k)^(-(k+1)/2))
# }
# I <- function(a,k){
#   ck<-sqrt(a^2*k/(k+1-a^2))
#   I1<-integrate(g,lower=0,upper=ck,k=k)$value
#   return (2*exp(log(gamma((k+1)/2))-log(gamma(k/2)))/sqrt(pi*k)*I1)
# }
I<-function(a,k){
  ck<-sqrt(a^2*k/(k+1-a^2))
  return (2*integrate(dt,lower=0,upper=ck,df=k)$value)
}
diff<-function(a,k){
  I(a,k)-I(a,k-1)
}
#为确定根的范围, 作出下列a-diff图示
k0<-c(4,8,15,25,100,500,1000)
a <- seq(0, sqrt(4)-1e-2, by=0.01)
res<-matrix(0,ncol=length(a),nrow=length(k0))
for (i in 1:length(k0)){
  for (j in 1:length(a)){
    res[i,j]=diff(a[j],k0[i])
  }
}
#分四组分别做出图示
#par(mfrow=c(2,2))
plot(a,res[1,],col=2,type="l",lwd=2,ylab="",xlab="")
lines(a,res[2,],col=3,lwd=2)
abline(h=0,col=4,lty=2)
legend("topright",col=2:3,cex=0.5,legend=c("k=4","k=8"),lwd=rep(2,2))
plot(a,res[3,],col=2,lwd=2,type="l",main="",ylab="",xlab="")
lines(a,res[4,],col=3,lwd=2)
abline(h=0,col=4,lty=2)
legend("topright",col=2:3,cex=0.5,legend=c("k=15","k=25"),lwd=rep(2,2))
plot(a,res[5,],col=2,type="l",lwd=2,main="",ylab="",xlab="")
abline(h=0,col=4,lty=2)
legend("topright",col=2:3,cex=0.5,legend="k=100",lwd=2)
plot(a,res[6,],col=2,type="l",lwd=2,main="",ylab="",xlab="")
lines(a,res[7,],col=3,lwd=2)
abline(h=0,col=4,lty=2)
legend("topright",col=2:3,cex=0.5,legend=c("k=500","k=1000"),lwd=rep(2,2))

#从上述图线中可以发现(1,2)之间均存在零点, 故可以猜想uniroot的搜索区间定为[1,2].
#为此给出k不同取值时在a=1,2处的值.

k<-c(4:25,100,500,1000)
for (i in 1:length(k)){
  print(round(c(k[i],diff(1,k[i]),diff(2,k[i])),7))
}
#可以看到, 在k不同取值时1,2处符号均相反. 故[1,2]中一定存在根.

roots<-numeric(length(k))
for (i in 1:length(k)){
  roots[i]<-uniroot(function(a){diff(a,k[i])},lower=1,upper=2)$root
}
knitr::kable((cbind(k,round(roots,3))))

##题11.4
S<-function(a,k){
  tk<-sqrt(a^2*k/(k+1-a^2))
  return (1-pt(tk,df=k))
}
A<-function(a,k){
  return (S(a,k)-S(a,k-1))
}
k<-c(4:25,100,500,1000)
roots2<-numeric(length(k))
for (i in 1:length(k)){
  roots2[i]<-uniroot(function(a){A(a,k[i])},lower=1,upper=2)$root
}
knitr::kable((cbind(k,round(roots2,3))))

## ----echo = TRUE--------------------------------------------------------------
y<-c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau<-1
x<-y[y!=1.00]
n<-length(y)
n1<-length(x)
n2<-n-n1
s<-sum(x)
lambda<-seq(1,100,1)
for (i in 2:100){
  lambda[i]<-(s+(n2*(lambda[i-1]+tau)))/n
  if (abs(lambda[i]-lambda[i-1])<1e-6){
    break
  }
}
name<-c("times","lambda")
knitr::kable(rbind(name,round(c(i,lambda[i]),4)))

## ----echo=FALSE---------------------------------------------------------------
lambda_star<-round(sum(y)/n1,4)

## ----echo=FALSE---------------------------------------------------------------
lambda_mle<-(s+n2*tau)/n1

## ----eval = TRUE--------------------------------------------------------------
set.seed(12345)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

# for循环
lm_for = vector("list", 10)
for (i in 1:10){
  x<-bootstraps[[i]]
  lm_for[[i]]<-lm(mpg ~ wt, data = x)
}
lm_for

# lapply方法
lm_lapply = lapply(bootstraps, function(x) lm(formula = mpg ~ wt, data = x))
lm_lapply

## ----echo = TRUE--------------------------------------------------------------
## (1)
set.seed(1234)
a <- matrix(sample(1:50,25), nrow=5, ncol=5)
a <- as.data.frame(a)
a_sd <- vapply(a, sd, FUN.VALUE=c(sde=0))
round(a_sd,3)
## (2)
a$V1[1]<-'a' #将a的v1的第一个分量改为字符
isNum <- vapply(a, is.numeric, logical(1)) #第一次vapply判断是否为数值型
a_sd <- vapply(a[isNum], sd, FUN.VALUE=c(0)) #第二次vapply计算sd
round(a_sd,3)

## ----echo = TRUE--------------------------------------------------------------
library(parallel)
mcsapply<-function(x,f,num=4){
  cl <- makePSOCKcluster(num) #num为核的数量
  g <-  parSapply(cl,x,f)
  stopCluster(cl)
  return(g)
}
mcsapply(seq(1,30,2), sqrt, 4)

mcvapply <- function(x, f, f.value, num=4) {
    cl <- makePSOCKcluster(num) #num为核的数量
    out_list <-  parSapply(cl,x,f)
    stopCluster(cl)
    out <- matrix(rep(f.value, length(x)), nrow = length(x))
    for (i in seq_along(x)) {
      res <- out_list[[i]]
      stopifnot(
      length(res) == length(f.value),
      typeof(res) == typeof(f.value)
      )
      out[i, ] <- res
      }
      out
    }
        
mcvapply(seq(1,30,2), sqrt,numeric(1), 4)

## ----echo = TRUE--------------------------------------------------------------
library(Rcpp) 

cppFunction('NumericMatrix gibbsC(int N, int thin, int n, int a, int b) {
  NumericMatrix mat(N, 2);
  double x = 10, y = 0.5; 
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x + a, n - x + b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}')

gibbsR <- function(N, thin, n, a, b) {
  mat <- matrix(nrow = N, ncol = 2)
  x <-10; y <- 0.5 #初值选为x = 10, y = 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, x + a, n - x + b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}

library(microbenchmark)
set.seed(1234)
a <- b <- 5;n <- 20;N <- 1000 #a,b,n,分别选为5,5,20; 生成数据总数为1000对
time<-microbenchmark(
  gibC <- gibbsC(N, 10, n, a, b),
  gibR <- gibbsR(N, 10, n, a, b)
)
qqplot(gibC[,1], gibR[,1])
abline(0,1)
qqplot(gibC[,2],gibR[,2])
abline(0,1)
summary(time)

