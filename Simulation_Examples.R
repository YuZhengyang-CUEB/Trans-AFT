# Simulation Examples
source('Simulation_Functions.R')
n0 = 150          # sample size of target cohort
nk = 200          # sample size of source cohorts
K = 8             # number of source cohorts
p = 500           # dimension
s = 20            # number of non-zero coefficients
sig.beta = 0.3    # non-zero coefficients of target cohort
sig.eta = 0.3     # exact transfer bias
exact = T         # exact configurations
h = 2            # number of exact transfer bias or degree of random transfer bias
heter = T         # heterogeneous design


Coef.list <- Coef.gen(s = s, h = h, K = K, sig.beta = sig.beta, sig.eta = sig.eta, p = p, exact = exact)
Xylist <- Data.gen(n0 = n0, nk = nk, K = K, T_c = T_c, S_c = S_c, p = p, Coef.list = Coef.list, heter = heter)
clist <- KM.Weight.Center(n0 = n0, nk = nk, K = K, Xylist = Xylist)
sse <- rep(NA, 3)
beta0 <- Coef.list$beta0
cv.las <- cv.glmnet(clist$X[1:n0, ], clist$y[1:n0], nfolds = 8, lambda = seq(1, 0.1,length.out = 10)*sqrt(2*log(p)/n0))
beta.las <- as.numeric(glmnet(clist$X[1:n0, ], clist$y[1:n0], lambda = cv.las$lambda.min)$beta)
sse[1] <- err.fun(beta.las, beta0)
beta.trans <- as.numeric(Trans.AFT(X = clist$X, y = clist$y, n0 = n0, nk = nk, K = K)$beta)
sse[2] <- err.fun(beta.trans, beta0)
cv.pool <- cv.glmnet(clist$X, clist$y, nfolds = 8, lambda = seq(1, 0.1, length.out = 10)*sqrt(2*log(p)/(n0 + nk*K)))
beta.pool <- as.numeric(glmnet(clist$X, clist$y, lambda = cv.pool$lambda.min)$beta)
sse[3] <- err.fun(beta.pool, beta0)
sse