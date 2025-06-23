# Fig 3 (b-c)
# regression scenario
U.DR = function(w, all.pred, all.Y){
  mean((all.Y - all.pred%*%w)^2)
}

U.WS = function(w, data.all, M = 5){
  # browser()
  folds = createFolds(1:nrow(data.all[[1]]$data), k = M)
  all.MSE = sapply(folds, function(fold){
    models = lapply(data.all, function(x){
      lm(Y~.-1, data = x$data[-fold,])
    })
    preds = do.call(rbind, lapply(data.all, function(x){
      cbind( x$data$Y[fold], do.call(cbind, predict(models, x$data[fold,])) )
    }))
    (preds[,1] - preds[,-1]%*%w)^2
  })
  mean(unlist(all.MSE))
}

U.CS = function(w, all.pred, all.Y){
  n = nrow(all.pred)/ncol(all.pred)
  K = ncol(all.pred)
  all.pred.z = (1-bdiag(rep(list(rep(1,n)),K)))*all.pred
  mean((all.Y - K/(K-1)*all.pred.z%*%w)^2)
}

U.tru = function(w, all.models, beta.mean, beta.sd){
  all.coef = sapply(all.models, coef)
  sum((beta.mean - all.coef%*%w)^2) + beta.sd^2*length(beta.mean) + 1
}

K = 20
p = 10

n = 100
beta.sd = 1
beta.mean = rep(1,p)

library(caret)

library(snowfall)
sfInit(parallel = T, cpus = 8)
sfLibrary(caret)
sfLibrary(Matrix)
sfExport("U.DR", "U.WS", "U.CS", "U.tru", "beta.mean", 
         "beta.sd", "K", "p", "n")
all.res = sfLapply(seq(0,0.5,0.05), function(w1){
  print(w1)
  w = c(w1, rep((1-w1)/(K-1), K-1))
  sapply(1:1e3, function(rep){
    data.all = lapply(1:K, function(k){
      X = matrix(rnorm(p*n),nrow = n)
      beta = rnorm(p, mean = beta.mean, sd = beta.sd)
      Y = X%*%beta + rnorm(n)
      list(data = data.frame(X=X, Y=Y), beta = beta)
    })
    
    all.models = lapply(data.all, function(x){
      lm(Y~.-1, data = x$data)
    })
    
    all.pred = do.call(rbind, lapply(data.all, function(x){
      do.call(cbind, predict(all.models, x$data))
    }))
    
    all.Y = unlist(lapply(data.all, function(x) x$data$Y))
    
    c(U.DR(w, all.pred, all.Y), U.WS(w, data.all), 
      U.CS(w, all.pred, all.Y), U.tru(w, all.models, beta.mean, beta.sd))
  })
})

all.sum = sapply(all.res, function(x){
  c(mean(x[1,] - x[4,]), mean(x[2,] - x[4,]), mean(x[3,] - x[4,]),
    sd(x[1,] - x[4,]), sd(x[2,] - x[4,]), sd(x[3,] - x[4,]))
})

plt.mean = data.frame(mean = c(all.sum[1,], all.sum[2,], all.sum[3,]), w1 = rep(seq(0,0.5,0.05),3),
                      method = factor(rep(c("DR", "CVWS", "CVCS"), each = 11), levels = c("DR", "CVWS", "CVCS"), ordered = T))

plt.sd = data.frame(mean = c(all.sum[4,], all.sum[5,], all.sum[6,]), w1 = rep(seq(0,0.5,0.05),3),
                    method = factor(rep(c("DR", "CVWS", "CVCS"), each = 11), levels = c("DR", "CVWS", "CVCS"), ordered = T))