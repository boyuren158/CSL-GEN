# prediction of mortality rate
library(randomForest)
library(caret)

# state-level analysis
all.data = read.csv("matched_pairs110091_final.csv", as.is = T)
all.states = unique(all.data.use$STATECODE)

bst.data = split(all.data, all.data$STATECODE)

ind.models = lapply( bst.data, function(x){
  randomForest::randomForest(mortality~.-STATECODE, data = x)
})

pred.all = lapply(bst.data, function(x){
  cbind(predict(ind.models, x))
})

all.res.state = lapply(1:20, function(k){
  # generalist
  # state train
  train.state = sample(seq_along(all.states), 10)
  test.state = setdiff(seq_along(all.states), train.state)
  
  n.each = sapply(bst.data[train.state], nrow)
  weight.c = sqrt(rep(n.each, n.each))
  
  # stacking design matrix
  stack.X = do.call(rbind, lapply(pred.all[train.state], function(tmp) tmp[,train.state]))
  stack.y = unlist(lapply( bst.data[train.state]), '[[', 'mortality')
  w_dr = nnls::nnls(A = stack.X/weight.c, b = stack.y/weight.c)$x

  # cvcs
  stack.X.cvcs = do.call(rbind, lapply(train.state, function(idx){
    tmp = pred.all[[idx]]
    tmp[,idx] = 0
    tmp[,train.state]
    }))
  
  w_cvcs = nnls::nnls(A = stack.X.cvcs/weight.c, b = stack.y/weight.c)$x*(1-n.each/sum(n.each))
 
  # evaluation
  eval.X = do.call(rbind, lapply(pred.all[test.state], function(tmp) tmp[,train.state]))
  eval.y = unlist(lapply( bst.data[test.state]), '[[', 'mortality')
  c( Metrics::mse(eval.y, c(eval.X %*% w_dr)),
     Metrics::mse(eval.y, c(eval.X %*% w_cvcs)) )
})

# county level analysis
zip_to_county = read.csv("ZIP-COUNTY-FIPS_2017-06.csv")
CA_county = zip_to_county %>% filter(STATE == "CA") %>% group_by(COUNTYNAME) %>% summarise(zip = ZIP[1])

all.data.CA = all.data %>% filter( STATECODE == "CA" )
all.data.CA.use = left_join(all.data.CA, CA_county, by = join_by( STATECODE == STATE ))
bst.data.CA = split(all.data.CA.use, all.data.CA.use$COUNTYNAME)

ind.models.CA = lapply( bst.data.CA, function(x){
  randomForest::randomForest(mortality~.-STATECODE, data = x)
})

pred.all.CA = lapply(bst.data.CA, function(x){
  cbind(predict(ind.models.CA, x))
})

all.county = unique(all.data.CA.use$COUNTYNAME)

all.res.county = lapply(1:20, function(k){
  # generalist
  # state train
  train.county = sample(seq_along(all.county), 10)
  test.county = setdiff(seq_along(all.county), train.county)
  
  n.each = sapply(bst.data.CA[train.county], nrow)
  weight.c = sqrt(rep(n.each, n.each))
  
  # stacking design matrix
  stack.X = do.call(rbind, lapply(pred.all.CA[train.county], function(tmp) tmp[,train.county]))
  stack.y = unlist(lapply( bst.data.CA[train.county]), '[[', 'mortality')
  w_dr = nnls::nnls(A = stack.X/weight.c, b = stack.y/weight.c)$x
  
  # cvcs
  stack.X.cvcs = do.call(rbind, lapply(train.county, function(idx){
    tmp = pred.all.CA[[idx]]
    tmp[,idx] = 0
    tmp[,train.county]
  }))
  
  w_cvcs = nnls::nnls(A = stack.X.cvcs/weight.c, b = stack.y/weight.c)$x*(1-n.each/sum(n.each))
  
  # evaluation
  eval.X = do.call(rbind, lapply(pred.all.CA[test.county], function(tmp) tmp[,train.county]))
  eval.y = unlist(lapply( bst.data.CA[test.county]), '[[', 'mortality')
  c( Metrics::mse(eval.y, c(eval.X %*% w_dr)),
     Metrics::mse(eval.y, c(eval.X %*% w_cvcs)) )
})