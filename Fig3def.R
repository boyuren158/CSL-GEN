y1 = rnorm(100, mean = -1)
y2 = rnorm(100, mean = -5)

M = 2
all.folds = lapply(1:M, function(x) 1:(100/M) + 10*(x-1))

res = do.call(rbind, lapply(all.folds, function(fold){
  y1.mean = mean(y1[-fold])
  y2.mean = mean(y2[-fold])
  cbind(y1[fold], y1.mean, y2.mean)
}) )

t(res[,2:3])%*%res[,2:3]/100

y.new = res[,1] - res[,3]
x.new = res[,2] - res[,3]

w = seq(0,1,0.01)
res = sapply(w, function(x)sum((y.new - x.new*x)^2))
sum(y.new*x.new)/sum(x.new^2)
plot(w, res)
lm(y.new~x.new-1)

lm((res[,1] - res[,3])~(res[,2] - res[,3])+1)


y1.mean.fold = sapply(all.folds, function(x){mean(y1[x])})
y1.mean.outfold = sapply(all.folds, function(x){mean(y1[-x])})
y2.mean.fold = sapply(all.folds, function(x){mean(y2[x])})

y1.mean = mean(y1.mean.fold)
y1.2.mean = mean(y1.mean.fold^2)
y12.mean = mean(y1.mean.fold*y2.mean.fold)
y2.mean = mean(y2.mean.fold)

y1.var = sd(y1.mean.fold)^2*(M-1)/M
y2.var = sd(y2.mean.fold)^2*(M-1)/M
y12.cov = mean(y1.mean.fold*y2.mean.fold) - y1.mean*y2.mean

y1.mean*y2.mean + 1/(M-1)^2*(y12.mean - y1.mean*y2.mean)

(M*y1.mean^2 + M/(M-1)^2*(y1.2.mean - y1.mean^2))/M


(M-1)/M*y1.mean.outfold[1] + 1/M*y1.mean.fold[1]

y1.mean^2 - y1.mean*y2.mean - 1/(M-1)*y1.var + M/(M-1)^2*y2.var - 1/(M-1)^2*y12.cov

# three panel figures to show difference between CS and WS cv
p = 20
n = 1e3
beta.all = cbind(rep(1,p), rep(1,p), rnorm(p, mean = 1,sd = 4))

data.all = lapply(1:3, function(i){
  x = MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = diag(p))
  y = x%*%beta.all[,i] + rnorm(n)
  data.frame(y = y, x)
})
beta.hat.all = sapply(data.all, function(x){
  coef(lm(y~.-1, data=x))
})

# level plot
# CVWS
U.ws.hat = function(data.all, w.all){
  all.folds = createFolds(1:n, k = 10)
  all.res = do.call(rbind, lapply(all.folds, function(fold){
    lm.all.cv = lapply(data.all, function(x){
      lm(y~.-1, data = x[-fold,])
    })
    do.call(rbind, lapply(data.all, function(x){
      cbind( x$y[fold], do.call(cbind, predict(lm.all.cv, x[fold,])) )
    }))
  }))
  y.all = all.res[,1]
  x.all = all.res[,-1]
  
  apply(w.all, 1, function(w){
    -mean((y.all - x.all%*%w)^2)
  })
}

w.fine = cbind(seq(0,1,0.001), 0, 1-seq(0,1,0.001))
U.ws.fine = U.ws.hat(data.all, w.fine)
w.opt = w.fine[order(U.ws.fine, decreasing = T)[1],1]
U.opt = sort(U.ws.fine, decreasing = T)[1]

w.all = cbind(seq(0,1,0.1), 0, 1-seq(0,1,0.1))
U.ws.all = U.ws.hat(data.all, w.all)
U.ws.plt = data.frame(w1 = c(w.all[,1],w.opt), Value = c(U.ws.all,U.opt), 
                      Type = factor(c(rep("Other", nrow(w.all)), "Maximum"), 
                                    levels = c("Other", "Maximum"), ordered = T))
p.level.ws = ggplot(data = U.ws.plt) + geom_abline(aes(intercept = w1, slope = -1, color = Value, linetype = Type) ) + coord_fixed() +
  scale_x_continuous(name = expression(w[1]), limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0,0)) +
  scale_y_continuous(name = expression(w[2]), limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0,0)) +
  theme_bw() + scale_color_viridis_c() + geom_abline(intercept = 1, slope = -1) +
  ggtitle("CVWS utility function")

#CVCS
U.cs.hat = function(data.all, w.all){
  lm.all = lapply(data.all, function(x){
    lm(y~.-1, data = x)
  })
  pred.all = do.call(rbind, lapply(data.all, function(x){
    do.call(cbind, predict(lm.all, x))
  }) )
  z = bdiag(rep(list(rep(1,n)),3))
  pred.all.z = pred.all*(1-z)
  y.all = unlist(lapply(data.all, function(x) x$y))
  # test.fn = function(w){
  #   sum((y.all - pred.all.z%*%(1.5*w))^2)
  # }
  #w.opt.cs.tru = alabama::auglag(c(1/3,1/3,1/3), test.fn, hin = function(w){w})#, heq = function(w){sum(w)-1})

  apply(w.all, 1, function(w){
    w.scale = 1.5*w
    -mean((y.all - pred.all.z%*%w.scale)^2)
  })
}

w.all = expand.grid(seq(0,1,0.01), seq(0,1,0.01))
w.all = cbind(w.all, 1-rowSums(w.all))
U.cs.all = U.cs.hat(data.all, w.all)
w.cs.opt = unlist(w.all[order(U.cs.all, decreasing = T)[1],])
# w.cs.opt.nn = w.opt.cs.tru$par
U.cs.plt = data.frame(w1 = w.all[,1], w2 = w.all[,2], value = U.cs.all)
cs.plt.init = ggplot(U.cs.plt, aes(x = w1, y = w2, z = value)) + geom_contour(aes(color = after_stat(level)), bins = 30)
cs.plt.data.all = ggplot_build(cs.plt.init)$data[[1]]
cs.plt.data.use = cs.plt.data.all[(cs.plt.data.all$x + cs.plt.data.all$y)<=1,]

p.level.cs = ggplot(cs.plt.data.use, aes(x = x, y = y, color = level, group = group)) + geom_path(aes(linetype = "Other")) +
  scale_color_viridis_c() + coord_fixed() + scale_x_continuous(expand = c(0,0), breaks = seq(0,1,0.1), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.1), limits = c(0,1)) + theme_bw() + xlab(expression(w[1])) +
  ylab(expression(w[2])) + geom_abline(intercept = 1, slope = -1) + 
  geom_point(x = w.cs.opt[1], y = w.cs.opt[2], aes(shape = "Maximum")) + labs(color = "Value", shape = "Type", linetype = "") +
  ggtitle("CVCS utility function")

# projection
x.new = MASS::mvrnorm(n = 1e3, mu = rep(0,p), Sigma = diag(p))

beta.ws = beta.hat.all%*%c(w.opt, 0, 1-w.opt)
beta.cs = beta.hat.all%*%w.cs.opt.nn

beta.pca.in = cbind(beta.hat.all, beta.ws, beta.cs, rep(1,p))
y.pca.in = x.new%*%beta.pca.in
y.pca.cor = cov(t(y.pca.in))
y.pca = princomp(covmat = y.pca.cor, cor = F)
y.scores = t(y.pca.in)%*%y.pca$loadings[,1:2]
apply(y.pca.in, 2, function(x) mean(x^2))

pca.plt = data.frame(x = 0, y = 0, xend = y.scores[,1], yend = y.scores[,2], 
                     Method = factor(c("betahat1", "betahat2", "betahat3", "beta.ws", "beta.cs", "beta0"),
                                     levels = c("betahat1", "betahat2", "betahat3", "beta.ws", "beta.cs", "beta0"),
                                     ordered = T))
p.pca = ggplot(data = pca.plt) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend,
                            color = Method), arrow = arrow(length = unit(0.5, "cm"))) +
  xlab("PC1") + ylab("PC2") + theme_bw() + scale_color_hue() #+ coord_fixed()
