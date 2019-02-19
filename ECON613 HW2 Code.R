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

#Question 1: calculate the correlation between Y and X1
corr_Y_X1 <- cor(Y,X1)
print(corr_Y_X1)

#Question 2: regression of Y on X
a <- data.frame(1,X1,X2,X3)
X <- data.matrix(a)
colnames(X) <- c("1","X1","X2","X3")
y <- data.matrix(Y)
eps <- data.matrix(eps)
#Question 3: calculate the coefficients on the regression
b <- (t(X)%*%X)
det(b)
beta = (solve(t(X)%*%X)%*%t(X)%*%y)
t(beta)

#Question 4 (1): calculate the standard errors using OLS
#estimate of sigma^2
sigma2_hat <- (t(eps)%*%eps) / (nrow(X)-ncol(X))
sigma2_hat
#estimate of variance_beta_X_hat
var_beta_X_hat <- c(sigma2_hat) * solve(t(X)%*%X)
var_beta_X_hat
#ols estimate of standard error
sqrt(diag(var_beta_X_hat))

#Question 4 (2): calculate the standard errors using bootstrap
#49 replications
#create an empty matrix to store sample values
boot_49 <- matrix(0,nrow=10000,ncol=4)
colnames(boot_49)<-c("X1","X2","X3","eps")
aa <- data.frame(X1,X2,X3,eps)
bb <- data.matrix(aa)
sdboot_49 <- matrix(0,nrow = 49,ncol=4)
for (i in 1:49) {
  for (j in 1:10000) {
    n <- sample(1:nrow(bb), size=1,replace = TRUE)
    cc <- matrix(bb[n,])
    boot_49[j,] <- t(cc)
  }
  xmatrix <- as.matrix(cbind(1,boot_49[,1],boot_49[,2],boot_49[,3]))
  epsmatrix <- as.matrix(boot_49[,4])
  #calculate stabdard error within the bootstrape
  sigma2_hat1 <- t(epsmatrix)%*%epsmatrix / (nrow(xmatrix)-ncol(xmatrix))
  var_beta_X_hat1 <- c(sigma2_hat1) *  solve(t(xmatrix)%*%xmatrix)
  cal_sd <- sqrt(diag(var_beta_X_hat1))
  sdboot_49[i,] <- cal_sd
}

#499 replication
boot_499 <- matrix(0,nrow=10000,ncol=4)
colnames(boot_499)<-c("X1","X2","X3","eps")
sdboot_499 <- matrix(0,nrow = 499,ncol=4)
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

#Question 1: write a function that returns the likelihood of the probit
beta1 <- as.matrix(solve(t(X)%*%X)%*%t(X)%*%ydum)
probit_like <- function(beta1,X.=X,ydum.=ydum){
  cdf <- pnorm(X.%*%beta1)
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log1 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(-log1)
}

#Question 2: implement the steepest ascent optimization algorithm to maximize likelihood
#define initial error term
e <- 0.0001
#set initial scaling parameter
alpha = 1
#set initial difference between old and new parameters
d = 1
#set initial old parameter
bo=c(0,0,0,0)
steepest_asent <- function(beta1,probit_like){
  #define tolarance level
  while(d > 0.000001){
    bo. <- matrix(bo,4,4)
    b1 <- bo. + diag(e,4)
    #fbo. <- probit_like(beta=bo.)
    #fb1. <- probit_like(beta=b1.)
    #plug in old and new parameters into probit likelihood function
    fbo. <- apply(bo.,2,probit_like)
    fb1 <- apply(b1,2,probit_like)
    #calculate the direction 
    direc <- (fb1-fbo.)/e
    #define optimized scaling parameter with backtracking armijo line-search
    p1 <- probit_like(beta1  = bo - alpha*direc)
    p2 <- probit_like(beta1  = bo) - 0.1*alpha*t(direc)%*%direc
    if (p1>p2){
      alpha <- 0.5*alpha
      next
    }
    bn <- bo - alpha*direc
    fbn <- probit_like(beta1 = bn)
    fbo <- probit_like(beta1 = bo)
    #fbn <- apply(bn,2,probit_like)
    #fbo <- apply(bo,2,probit_like)
    d = fbn - fbo
    bo <- bn
  }
  return(bo)
}
steepest_asent(bo,probit_like)

#Question 3: how different are the paraameters from the true parameters? See in the output document 


#Exercise 4 Discrete Choice

#Question 1: write and optimize the probit,logit, and the linear probability model
#probit likelihood function
probit_likelihood <- function(beta1,X.=X,ydum.=ydum){
  cdf <- pnorm(X.%*%beta1)
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log2 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(-log2)
}
#probit model optimization
probit_optim <- optim(par = beta1,probit_likelihood)
#construct the probit optimization coefficient 
probit_optim_beta <- probit_optim$par
print(probit_optim_beta)

#logit likelihood function
logit_likelihood <- function(beta1,X.=X,ydum.=ydum){
  cdf <- exp(X.%*%beta1)/(1+exp(X.%*%beta1))
  #eliminate extreme values to be capable of plugging into log function
  cdf[cdf>0.99999] <- 0.99999
  cdf[cdf<0.00001]<- 0.00001
  log3 <- sum(ydum.*log(cdf))+sum((1-ydum.)*log(1-cdf))
  return(log3)
}
#logit model optimization
logit_optim <- optim(par = beta1,logit_likelihood)
#construct the logit optimization coefficient 
logit_optim_beta <- logit_optim$par
print(logit_optim_beta)

#linear probability model optimization
linear_model <- function(beta1,X.=X,ydum.=ydum){
  ydum_hat <- X%*%beta1
  SSR <- sum((ydum_hat-ydum)^2)
  return(SSR)
}
#linear probability optimization
linear_optim <- optim(par = beta1,linear_model,hessian = TRUE) 
#construct the linear optimization coefficient 
linear_optim_beta <- linear_optim$par
print(linear_optim_beta)

#Question 2: interpret and compare the estimated coefficients(see in the output document) and test on the significance
#test significance
dataset <- data.frame(ydum,X1,X2,X3)
coef_probit <- glm(ydum~X1+X2+X3, family = binomial(link = "probit"),data = dataset)
summary(coef_probit)
coef_logit <- glm(ydum~X1+X2+X3, family = binomial,data = dataset)
summary(coef_logit)
coef_linear <- lm(ydum~X1+X2+X3,data = dataset )
summary(coef_linear)

#Exercise 5 
install.packages("numDeriv")
library(numDeriv)
#establish a function to calculate the marginal effect by probit model
X <- as.matrix(cbind(1,X1,X2,X3))
marginal_probit <- function(beta1,X){
  pmulti <- X%*%beta1
  ppdf <- dnorm(pmulti)
  me_probit <-ppdf%*%t(beta1)
  return(me_probit)
}
ME_probit <- marginal_probit(probit_optim_beta,X)

#establish a function to calculate the marginal effect by logit model
marginal_logit <- function(beta1,X){
  lmulti <- X%*%beta
  lpdf <- (exp(lmulti)/((1+exp(lmulti))^2))
  me_logit <-lpdf%*%t(beta)
  return(me_logit)
}
ME_logit <- marginal_logit(logit_optim_beta,X)

write.csv(ME_probit,file = "~/Desktop/Probit Marginal Effect.csv", row.names = FALSE)
write.csv(ME_logit,file = "~/Desktop/Logit Marginal Effect.csv", row.names = FALSE)

#Compute the standard deviations of marginal effects using bootstrap
#probit
boot_490 <- matrix(0,nrow=10000,ncol=4)
pmargin_490 <- matrix(0,nrow = 490,ncol=4)
xmean <- matrix(0,nrow = 490,ncol = 4)
dd <- data.frame(X1,X2,X3,ydum)
ee <- data.matrix(dd)
for (i in 1:490){ 
  for (j in 1:10000){
    f <- sample(1:nrow(ee),size = 1, replace = TRUE)
    g <- matrix(ee[f,])
    boot_490[j,] <- t(g)
  }
  mean_x0 <- as.matrix(mean(1))
  mean_x1 <- as.matrix(mean(boot_490[,1]))
  mean_x2 <- as.matrix(mean(boot_490[,2]))
  mean_x3 <- as.matrix(mean(boot_490[,3]))
  mean_x <- as.matrix(cbind(mean_x0,mean_x1,mean_x2,mean_x3))
  xmean[i,] <- as.matrix(mean_x)
  cal_margin <- as.matrix(pnorm(xmean[i,]%*%probit_optim_beta)%*%t(probit_optim_beta))
  pmargin_490[i,] <- cal_margin
}
sdboot_probit <- as.matrix(apply(pmargin_490, 2,sd))

#logit
lmargin_490 <- matrix(0,nrow = 490,ncol=4)
for (i in 1:490){ 
  for (j in 1:10000){
    f <- sample(1:nrow(ee),size = 1, replace = TRUE)
    g <- matrix(ee[f,])
    boot_490[j,] <- t(g)
  }
  mean_x0 <- as.matrix(mean(1))
  mean_x1 <- as.matrix(mean(boot_490[,1]))
  mean_x2 <- as.matrix(mean(boot_490[,2]))
  mean_x3 <- as.matrix(mean(boot_490[,3]))
  mean_x <- as.matrix(cbind(mean_x0,mean_x1,mean_x2,mean_x3))
  xmean[i,] <- as.matrix(mean_x)
  cal_marginn <- (exp(xmean[i,]%*%logit_optim_beta)/((1+exp(xmean[i,]%*%logit_optim_beta))^2))%*%t(logit_optim_beta)
  lmargin_490[i,] <- cal_marginn
}
sdboot_logit <- as.matrix(apply(lmargin_490, 2,sd))

#compute the standard error of marginal effects using the delta method 
#calculate the variance covariance matrix for probit/logit
myprobit <- glm(ydum~X1+X2+X3,family = binomial(link = "probit"),data = dataset)
v_probit <- vcov(myprobit)
mylogit <- glm(ydum~X1+X2+X3,family = binomial,data = dataset)
v_logit <- vcov(mylogit)
