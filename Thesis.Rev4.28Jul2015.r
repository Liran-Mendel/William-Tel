#### This is script for simulating one-parameters logistic regression Y = b1*X (b0=0).

# Revision 4, for meeting on 03Aug2015.
# Objective: See how the variance of the beta estimators (b0,b1) changes given
# different sigmas and p's.


# Revision 3, for meeting 28/07/2015.
# 1. Denser grid for p=0.01.
# 1. Change SD to fixed, without overlapping between ranges of different p.
# Value = delta=difference in grid (0.01).

# Revision 2, post meeting 03Feb2015.
# Goals: 
# 1. Simulate information matrix. 
# 2. Simulate information  matrix determinant. 
# 3. Evaluate 2-nd order Taylor approximation (by MSE?).

# Reset your vectors (FI=Fisher Information):
set.seed(2015)


# Set parameter values in a grid (beta's, p, n):
b0 <- 0
b1 <- 1
g.dens <- 0.01   # grid density
p <- seq(0.51, 0.99, by=g.dens) 

# test different p's, p=0.81 optimal logistic regression, non-random x
n <- 40                          # number of runs = 40 for each trial
mu.x.p <- (log(p/(1-p))-b0)/b1  # true mu.x=(p)ercentile.
sigma <- c(g.dens/5:1) # Fixed sigma at grid density to prevent overlap.

# Creat 3 vecs for: b0.hat, b1.hat, det(vcov) with dimensions (p, sigma, Mean/Sigma)
array.names <- list(paste("p=", p), paste("sigma= ",sigma), c("Mean", "Sigma"))
b0.vec <- array(NA, dim=c(length(p), length(sigma), 2 ), dimnames=array.names)
b1.vec <- array(NA, dim=c(length(p), length(sigma), 2 ), dimnames=array.names)
det.vcov.vec <- array(NA, dim=c(length(p), length(sigma), 2 ), dimnames=array.names)

#### 500 simulations for each value in the grid (p, sigma).
for (j in 1:length(p))      # use every non-noisy "p"
{
  for (k in 1:length(sigma))    # for each p, try all sigma's
  {
    
    # Zero vectors.
    b0.hat <- c()
    b1.hat <- c()
    det.vcov.hat <- c()
    for (i in 1:500)  # 500 iterations.
     {
       # Noisey x. 20 for p, 20 for 1-p.
       # Sigma=varies, mu.x.p=p percentile, 50% on p based on D-design
       x.p <- rnorm(n*0.5, mean=mu.x.p[j], sd=sigma[k])
       # Sigma=varies, mu.x.p=x.p percentile, 50% on 1-p based on D-design, 
       # Logistic is anti-symmetric (ln(p/1-p)=-ln(p/(1-p))) 1.p=1 minus p.
       x.1.p <- rnorm((n - length(x.p)), mean=-mu.x.p[j], sd=sigma[k])
       
       # Based on observed (noisy) x, calculate prob for Y=1.
       prob.x.p <- exp(b0+b1*x.p)/(1+exp(b0+b1*x.p))
       prob.x.1.p <- exp(b0+b1*x.1.p)/(1+exp(b0+b1*x.1.p))
       
       # Bernoulli success rate for each x.       
       y.x.p <- rbinom(length(prob.x.p), 1, prob.x.p)
       y.x.1.p <- rbinom(length(prob.x.1.p), 1, prob.x.1.p)
       
       # Join data for glm and estimate betas.
       logit.x <- c(x.p, x.1.p)
       logit.x
       logit.y <- c(y.x.p, y.x.1.p)
       logit.y
       mydata <- data.frame(cbind(logit.y, logit.x))
       mylogit <- glm(logit.y ~ logit.x, data = mydata, family = "binomial")
       summary(mylogit)
       # vector of all 500 simulations for each (p, sigma)
       b0.hat <- c(b0.hat, mylogit$coefficients[1])
       b1.hat <- c(b1.hat, mylogit$coefficients[2])
       det.vcov.hat <- c(det.vcov, det(vcov(mylogit)))
    }
       # create for each sigma, p, vector for b0 and b1, and vcov for 500 runs.
       b0.vec[j,k,1] <- mean(b0.hat); b0.vec[j,k,2] <- sd(b0.hat)
       b1.vec[j,k,1] <- mean(b1.hat); b1.vec[j,k,2] <- sd(b1.hat)
       det.vcov.vec[j,k,1] <- mean(det.vcov.hat); det.vcov.vec[j,k,2] <- sd(det.vcov.vec);  
  }
}



# Plot for various beta data (beta.dat) given sigma = function of p(ercentile)
p.varB.plot <- function(p, sigma, sig.col, beta.dat)
{
  lines(beta.dat ~ p, col=sig.col)
}

# bounded values, for better graphic representation
bound <- function(min, max, dat.vec) 
  return(which((dat.vec <= max) & (dat.vec >=min))) 

# For beta0.

trunc.vec <- matrix(nrow=length(sigma), ncol=length(p))
for (i in 1:length(sigma))
{
  # Upper bound by mean + 5*sd
  trunc <- bound (min=0 ,max=mean(b0.vec[,i,2])+5*sd(b0.vec[,i,2]), 
                  b0.vec[,i,2])
  trunc.vec[i] <- trunc
  # i=1 plot, then lines
  if (i==1) plot((b0.vec[trunc,i,2]) ~ p[trunc], main="VAR(b0)=f(sigma, p)",
   sub="upper trunc <= mean + 5*SD", ylab="VAR(beta0)") 
  else p.varB.plot(p[trunc], sigma[i], i, as.vector(b0.vec[trunc,i,2]))
   
}
legend("topleft", as.expression(sigma), pch=19, 
       col=1:length(sigma),bty="n", y.intersp=rep(0.5,length(sigma)))

# Print outliers
for (i in 1:5)
{
  cat("sigma=",sigma[i])
  a <- cbind(p[-trunc], b0.vec[-trunc,i,2])
  colnames(a) <- c("p", "var.comp")
  print(a)
}

# Hmm 0.97 is suspicious!

    
# For beta1.

trunc.vec <- matrix(nrow=length(sigma), ncol=length(p))
for (i in 1:length(sigma))
{
  # Upper bound by mean + 5*sd
  trunc <- bound (min=0 ,max=mean(b1.vec[,i,2])+5*sd(b1.vec[,i,2]), 
                  b1.vec[,i,2])
  trunc.vec[i] <- trunc
  # i=1 plot, then lines
  if (i==1) plot((b1.vec[trunc,i,2]) ~ p[trunc], main="VAR(b1)=f(sigma, p)",
                 sub="upper trunc <= mean + 5*SD", ylab="VAR(beta1)") 
  else p.varB.plot(p[trunc], sigma[i], i, as.vector(b1.vec[trunc,i,2]))
  
}
legend("topleft", as.expression(sigma), pch=19, 
       col=1:length(sigma),bty="n", y.intersp=rep(0.5,length(sigma)))

# Print outliers
for (i in 1:5)
{
  cat("sigma=",sigma[i])
  a <- cbind(p[-trunc], b1.vec[-trunc,i,2])
  colnames(a) <- c("p", "var.comp")
  print(a)
}

# For det.vcov

trunc.vec <- matrix(nrow=length(sigma), ncol=length(p))
for (i in 1:length(sigma))
{
  # Upper bound by mean + 5*sd
  trunc <- bound (min=0 ,max=mean(det.vcov.vec[,i,2])+5*sd(det.vcov.vec[,i,2]), 
                  det.vcov.vec[,i,2])
  trunc.vec[i] <- trunc
  # i=1 plot, then lines
  if (i==1) plot((det.vcov.vec[trunc,i,2]) ~ p[trunc], main="VAR(det.vcov)=f(sigma, p)",
                 sub="upper trunc <= mean + 5*SD", ylab="VAR(det.vcov)") 
  else p.varB.plot(p[trunc], sigma[i], i, as.vector(det.vcov.vec[trunc,i,2]))
  
}
legend("topleft", as.expression(sigma), pch=19, 
       col=1:length(sigma),bty="n", y.intersp=rep(0.5,length(sigma)))

# Print outliers
for (i in 1:5)
{
  cat("sigma=",sigma[i])
  a <- cbind(p[-trunc], det.vcov.vec[-trunc,i,2])
  colnames(a) <- c("p", "var.comp")
  print(a)
}



save(list = ls(all = TRUE), file = "WilliamTell.11Aug2015.RData")


