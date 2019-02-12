                             #Exercise 1 Data Creation
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
n = 10000
R = 49
result = rep(X, R)
for (i in 1:R) {
  boot_sample = sample(X, n, replace = TRUE)
  result[i] = sd(solve(t(X)%*%X)%*%t(X)%*%y)
}
result
#sqrt(diag(var((solve(t(X)%*%X)%*%t(X)%*%y))))
#R=49
#R=499
#estimated_beta=()
#for i=1:R
