library(randomForest)
library(caret)
library(Matrix)
# library(nnls)
library(glmnet)

# p: number of covariates
# n: sample size in each study
p = 10
n = 2e2

beta.mean = rep(1,p)

# stacking with l1 constraints
get.l1res = function( Y, X ){
  #Y is the obs values and X are columns of predicted values
  Exp.mt = cbind( diag(rep(1,ncol(X))), diag(rep(-1,ncol(X))), 0 )
  A = X%*%Exp.mt
  d = t(A)%*%Y
  Amat = matrix(c(rep(-1,2*ncol(X)),1),nrow=1)
  
  capture.output( res <- LowRankQP::LowRankQP(Vmat = t(A)%*%A, dvec = -d, Amat = Amat, bvec = 0,
                                   uvec = rep(1, 2*ncol(X)+1), method = "LU", verbose = F) )
  all.coef = res$alpha[-length(res$alpha)]
  all.coef[1:(length(all.coef)/2)] - all.coef[(length(all.coef)/2+1):length(all.coef)]
}

# stcaking with ridge penalty
ridge.w = function(x, y){
  cv.res = cv.glmnet(x=x, y=y, alpha=0, intercept=F, nfolds = 5, 
            lambda = exp(seq(log(0.001),log(5),length.out = 20)))
  glmnet.res = glmnet(x=x, y=y, alpha=0, intercept = F, lambda = cv.res$lambda.min)
  glmnet.res$beta
}
lasso.w = function(x, y){
  cv.res = cv.glmnet(x=x, y=y, alpha=1, intercept=F, nfolds = 5, 
                     lambda = exp(seq(log(0.001),log(5),length.out = 20)))
  glmnet.res = glmnet(x=x, y=y, alpha=1, intercept = F, lambda = cv.res$lambda.min)
  glmnet.res$beta
}

# get the simulation results for all scenarios
# K: number of studies
# beta.mean: mean of beta coef across studies
# beta.sd: sd of beta coef across studies (same for different beta components)
# n.rep: number of replicates
# output - average prediction MSE of regular and zero-out stacking
sim.fn = function(n, K, beta.mean, beta.sd, n.rep = 100){
  test.res <- sapply(1:n.rep, function(rep){
    if(rep%%5==0){
      print(rep)
    }
    # simulate multi-study data
    all.data = lapply(1:K, function(k){
      beta = rnorm(p, sd = beta.sd) + beta.mean
      X = matrix(rnorm(p*n), ncol = p)
      y = X%*%beta + rnorm(n, sd = 1)
      list(data=data.frame(y=y,X=X),beta=beta)
    })
    
    y.all = unlist( lapply(all.data, function(x) x$data$y) )
    # y.all.cv = unlist( lapply(all.data, function(x) x$y[(3*n/4+1):n]) )
    X.all = do.call( rbind, lapply(all.data, function(x) model.matrix(~.-y-1,x$data)) )
    
    # study-specific models
    all.mods = lapply(all.data, function(ind){
      lm(y~.-1, data = ind$data)#randomForest(y~., data = ind)
    })
    
    # all.coef = sapply(all.mods, coef)
    all.coef = sapply(all.data, function(x) x$beta)
    all.pred = X.all%*%all.coef
    
    all.pred.z = all.pred*(1-bdiag(rep(list(rep(1,n)),K)))
    
    # use l1 constraint (|w|_1 \leq 1) on stacking weights
    w.tru = get.l1res(beta.mean, all.coef)
    w.r = get.l1res(y.all, all.pred)#lasso.w(x = all.pred, y = y.all)
    w.0 = get.l1res(y.all, as.matrix(all.pred.z)*K/(K-1))#lasso.w(x = all.pred.z*K/(K-1), y = y.all)
    #w.0.new = coef(lm(y.all~as.matrix(all.pred.z)-1))*(K-1)/K
    
    c(sum((all.coef%*%w.r - beta.mean)^2), sum((all.coef%*%w.0 - beta.mean)^2),
      sum((all.coef%*%w.tru - beta.mean)^2))
    #sum((all.coef%*%w.0.new - beta.mean)^2))
    
    # y.test = unlist( lapply(test.data, function(x) x$y) )
    # X.test = do.call( rbind, lapply(test.data, function(x) model.matrix(~.-y-1,x)) )
    # # X.all.cv = do.call( rbind, lapply(all.data, function(x) model.matrix(~.-y-1,x)[(3*n/4+1):n,]) )
    # beta.merge = coef(lm(y.all~X.all-1))
  })
  # average over replicates
  rowMeans(test.res) 
}

p = 10
beta.mean = rep(1,p)
config1 = expand.grid(K = c(5, 15, 20), n = seq(20, 100, 5))
config2 = expand.grid(K = 5:50, n = c(100,200,400))
sfInit(parallel = T, cpus = 8)
sfLibrary(LowRankQP)
sfLibrary(caret)
sfLibrary(Matrix)
sfExport("p", "beta.mean", "get.l1res", "sim.fn", "config1", "config2")
config1.res = sfApply(config1, 1, function(x){
  sim.fn(n = x[2], K = x[1], beta.mean = beta.mean, beta.sd = 1, n.rep = 2e3)
})
config2.res = sfApply(config2, 1, function(x){
  sim.fn(n = x[2], K = x[1], beta.mean = beta.mean, beta.sd = 1, n.rep = 2e3)
})
sfStop()

# K constant, n increases
dr.diff = config1.res[1,] - config1.res[3,]
cs.diff = config1.res[2,] - config1.res[3,]
plt.data = config1
plt.data$dr = dr.diff
plt.data$cs = cs.diff
plt.data$var = 1/sqrt(plt.data$n)

plt.fit = lapply(split(plt.data, plt.data$K), function(x){
  data.frame(dr.fit = lm(dr~var, data = x)$fitted.values, 
             cs.fit = lm(cs~var, data = x)$fitted.values)
})

plt.data$cs.fit = c(t(sapply(plt.fit, function(x) x$cs.fit)))
plt.data$dr.fit = c(t(sapply(plt.fit, function(x) x$dr.fit)))
p.dr.n = ggplot(plt.data) + geom_point(aes(x = n, y = dr, color = as.factor(K), shape = "Approximated")) +
  geom_line(aes(x = n, y = dr.fit, color = as.factor(K), linetype = "Fitted")) + 
  theme_bw() + ylab("Edr - E0") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.7))
p.cs.n = ggplot(plt.data) + geom_point(aes(x = n, y = cs, color = as.factor(K), shape = "Approximated")) +
  geom_line(aes(x = n, y = cs.fit, color = as.factor(K), linetype = "Fitted")) + 
  theme_bw() + ylab("Ecs - E0") + xlab("K") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.7))

# n constant, K increases
dr.diff = config2.res[1,] - config2.res[3,]
cs.diff = config2.res[2,] - config2.res[3,]
plt.data = config2
plt.data$dr = dr.diff
plt.data$cs = cs.diff
plt.data$var = log(plt.data$K)/plt.data$K#1/sqrt(plt.data$n)

plt.data = plt.data[plt.data$K>=15,]

plt.fit = do.call(rbind, lapply(split(plt.data, plt.data$n), function(x){
  data.frame(dr.fit = lm(dr~var, data = x)$fitted.values, 
             cs.fit = lm(cs~var, data = x)$fitted.values)
}))

plt.data$cs.fit = plt.fit$cs.fit
plt.data$dr.fit = plt.fit$dr.fit
p.dr.K = ggplot(plt.data) + geom_point(aes(x = K, y = dr, color = as.factor(n), shape = "Approximated")) +
  geom_line(aes(x = K, y = dr.fit, color = as.factor(n), linetype = "Fitted")) + 
  theme_bw() + ylab("Edr - E0") + xlab("K") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.7))
p.cs.K = ggplot(plt.data) + geom_point(aes(x = K, y = dr, color = as.factor(n), shape = "Approximated")) +
  geom_line(aes(x = K, y = dr.fit, color = as.factor(n), linetype = "Fitted")) + 
  theme_bw() + ylab("Ecs - E0") + xlab("K") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.85,0.7))

# comparisons of cvcs to DR
# small p case
all.res = sapply(2:50, function(K){
  sapply(seq(0,4,0.2), function(beta.sd){
    res = sim.fn(K, beta.mean, beta.sd)
    # calculate the difference
    res[1] - res[2]
  })
})

# large p case
beta.mean = rep(1,1e3)
p = 1e3
all.res.highp = sapply(2:50, function(K){
  sapply(seq(0,4,0.2), function(beta.sd){
    res = sim.fn(K, beta.mean, beta.sd)
    # calculate the difference
    res[1] - res[2]
  })
})