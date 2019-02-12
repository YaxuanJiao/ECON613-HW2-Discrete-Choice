#Exercise 1 Data Creation
set.seed(100)
#construct X1
X1 <- as.vector(runif(10000,1,3)) 
#construct X2
X2 <- as.vector(rgamma(10000,3,2))
#construct X3
X3 <- as.vector(rbinom(10000,1,0.3))
#construct eps
eps <- as.vector(rnorm(10000,2,1))
#Y=0.5+1.2X1-0.9X2+0.1X3+esp 
Y <- 0.5 + 1.2*X1 - 0.9*X2 + 0.1*X3 + eps
ydum <- ifelse(Y>mean(Y),1,0)

#Exercise 2 OLS
#calculate the correlation between Y and X1
corr_Y_X1 <- cor(Y,X1)
print(corr_Y_X1)

#regression of Y on X
a <- data.frame(1,X1,X2,X3)
X <- data.matrix(a)
colnames(X) <- c("1","X1","X2","X3")
y <- data.matrix(Y)
eps <- data.matrix(eps)
#calculate the coefficients on the regression
b <- (t(X)%*%X)
det(b)
beta = (solve(t(X)%*%X)%*%t(X)%*%y)
t(beta)

#calculate the standard errors using OLS
#estimate of sigma^2
sigma2_hat <- (t(eps)%*%eps) / (nrow(X)-ncol(X))
sigma2_hat
#estimate of variance_beta_X_hat
var_beta_X_hat <- c(sigma2_hat) * solve(t(X)%*%X)
var_beta_X_hat
#ols estimate of standard error
sqrt(diag(var_beta_X_hat))

#calculate the standard errors using bootstrap
#49 replications
boot_49 <- matrix(nrow=10000,ncol=4)
colnames(boot_49)<-c("X1","X2","X3","eps")
aa <- data.frame(X1,X2,X3,eps)
bb <- data.matrix(aa)
sdboot_49 <- matrix(nrow = 49,ncol=4)
for (i in 1:49) {
  for (j in 1:10000) {
    n <- sample(1:nrow(bb), size=1,replace = TRUE)
    cc <- matrix(bb[n,])
    boot_49[j,] <- t(cc)
  }
  xmatrix <- as.matrix(cbind(1,boot_49[,1],boot_49[,2],boot_49[,3]))
  epsmatrix <- as.matrix(boot_49[,4])
  sigma2_hat1 <- t(epsmatrix)%*%epsmatrix / (nrow(xmatrix)-ncol(xmatrix))
  var_beta_X_hat1 <- c(sigma2_hat1) *  solve(t(xmatrix)%*%xmatrix)
  cal_sd <- sqrt(diag(var_beta_X_hat1))
  sdboot_49[i,] <- cal_sd
}

#499 replication
boot_499 <- matrix(nrow=10000,ncol=4)
colnames(boot_499)<-c("X1","X2","X3","eps")
aa <- data.frame(X1,X2,X3,eps)
bb <- data.matrix(aa)
sdboot_499 <- matrix(nrow = 499,ncol=4)
for (i in 1:499) {
  for (j in 1:10000) {
    n <- sample(1:nrow(bb), size=1,replace = TRUE)
    cc <- matrix(bb[n,])
    boot_499[j,] <- t(cc)
  }
  xmatrix <- as.matrix(cbind(1,boot_499[,1],boot_499[,2],boot_499[,3]))
  epsmatrix <- as.matrix(boot_499[,4])
  sigma2_hat1 <- t(epsmatrix)%*%epsmatrix / (nrow(xmatrix)-ncol(xmatrix))
  var_beta_X_hat1 <- c(sigma2_hat1) *  solve(t(xmatrix)%*%xmatrix)
  cal_sd <- sqrt(diag(var_beta_X_hat1))
  sdboot_499[i,] <- cal_sd
}

#Exercise 3
#Implement the steepest ascent optimization algorithm to maximize that likelihood
#probit likelihood function
probit_like <- function(beta,X.=X,ydum.=ydum){
  beta<-as.matrix(beta)
  cdf <- pnorm(X.%*%beta)
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log1 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(-log1)
}
#build up steepest function
#define error term
e <- 0.0001
#set initial scaling parameter
alpha = 1
#set initial difference between old and new parameters
d = 1
#set initial old parameter
bo=c(0,0,0,0)
steepest_asent <- function(beta,probit_like){
  #define tolarance level
  while(d > 0.000001){
    bo. <- matrix(bo, 4,4)
    b1 <- bo. + diag(e,4)
    #fbo. <- probit_like(beta=bo.)
    #fb1. <- probit_like(beta=b1.)
    #plug in old and new parameters into probit likelihood function
    fbo. <- apply(bo.,2,probit_like)
    fb1 <- apply(b1,2,probit_like)
    #calculate the direction 
    direc <- (fb1-fbo.)/e
    #define optimized scaling parameter with backtracking armijo line-search
    p1 <- probit_like(beta = bo - alpha*direc)
    p2 <- probit_like(beta = bo) - 0.1*alpha*t(direc)%*%direc
    if (p1>p2){
      alpha <- 0.5*alpha
      next
    }
    bn <- bo - alpha*direc
    fbn <- probit_like(beta = bn)
    fbo <- probit_like(beta = bo)
    #fbn <- apply(bn,2,probit_like)
    #fbo <- apply(bo,2,probit_like)
    d = fbn - fbo
    bo <- bn
  }
  return(bo)
}
steepest_asent(bo,probit_like)

#Exercise 4 Discrete Choice
#probit likelihood function
beta <- as.matrix(beta)
probit_likelihood <- function(beta,X.=X,ydum.=ydum){
  cdf <- pnorm(X.%*%beta)
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log2 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(-log2)
}
#probit model optimization
probit_optim <- optim(par = beta,probit_likelihood)

#logit likelihood function
logit_likelihood <- function(beta,X.=X,ydum.=ydum){
  cdf <- exp(X.%*%beta)/(1+exp(X.%*%beta))
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log3 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(log3)
}
#logit model optimization
logit_optim <- optim(par = beta,logit_likelihood)

#linear probability model optimization
linear_model <- function(belta,X.=X,ydum.=ydum){
  ydum_hat <- X%*%beta
  SSR <- sum((ydum_hat-ydum)^2)
  return(SSR)
}
#linear probability optimization
linear_optim <- optim(par = beta,linear_model,hessian = TRUE) 

#test significance
dataset <- data.frame(ydum,X1,X2,X3)
coef_probit <- glm(ydum~X1+X2+X3, family = binomial(link = "probit"),data = dataset)
summary(coef_probit)
coef_logit <- glm(ydum~X1+X2+X3, family = binomial(link = "logit"),data = dataset)
summary(coef_logit)
coef_linear <- lm(ydum~X1+X2+X3,data = dataset )
summary(coef_linear)

#Exercise 5 
install.packages("numDeriv")
library(numDeriv)
#establish a function to calculate the marginal effect by probit model
beta <- as.matrix(beta)
X <- as.matrix(cbind(1,X1,X2,X3))
marginal_probit <- function(beta.,X.){
  pmulti <- X.%*%beta.
  ppdf <- dnorm(pmulti)
  me_probit <-ppdf%*%t(beta.)
  return(me_probit)
}
ME_probit <- marginal_probit(probit_optim$par,X)

#establish a function to calculate the marginal effect by logit model
marginal_logit <- function(beta.,X.){
  lmulti <- X.%*%beta.
  lpdf <- (exp(lmulti)/((1+exp(lmulti))^2))
  me_lrobit <-lpdf%*%t(beta.)
  return(me_lrobit)
}
ME_logit <- marginal_logit(logit_optim$par,X)

write.csv(ME_probit,file = "~/Desktop/Probit Marginal Effect.csv", row.names = FALSE)
write.csv(ME_logit,file = "~/Desktop/Logit Marginal Effect.csv", row.names = FALSE)

#Compute the standard deviations using bootstrap
boot_490 <- matrix(nrow=10000,ncol=5)
margin_490 <- matrix(nrow = 490,ncol=4)
colnames(boot_490)<-c("constant","X1","X2","X3","ydum")
dd <- data.frame(1,X1,X2,X3,ydum)
ee <- data.matrix(dd)
for (m in 1:490){ 
  for (n in 1:10000){
    f <- sample(1:nrow(ee),size = 1, replace = TRUE)
    g <- matrix(ee[f,])
    boot_490[n,] <- t(g)
  }
  datax <- as.matrix(cbind(boot_490[,1],boot_490[,2],boot_490[,3],boot_490[,4]))
  dataydum <- as.matrix(boot_490[,5])
  mean_x0 <- as.matrix(mean(boot_490[,1]))
  mean_x1 <- as.matrix(mean(boot_490[,2]))
  mean_x2 <- as.matrix(mean(boot_490[,3]))
  mean_x3 <- as.matrix(mean(boot_490[,4]))
  mean_x <- cbind(mean_x0,mean_x1,mean_x2,mean_x3)
  cal_coef <- (solve(t(boot_490[,1:4])%*%boot_490[,1:4])%*%t(boot_490[,1:4])%*%boot_490[,5])
  nie <- mean_x%*%cal_coef
  cal_margin <- cal_coef%*%nie
  margin_490[m,] <- cal_margin 
}
re = as.function()
sd(margin_490)
sdboot_490 <- as.matrix(apply(margin_490, 2,sd))
