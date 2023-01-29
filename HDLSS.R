pval = function(X,Y,R,t, D=NULL)
{
  #if(t==1){return(bd_test(X,Y,R))}
  if(t==2){return(bd_test_gdist(D,R))}
  if(t==3){return(FRTest(X,Y,R))}
  if(t==4){return(BFTest(X,Y,R))}
  if(t==5){return(knnTest(X,Y,3,R))}
  if(t==6){return(mmd(X,Y,R))}
  if(t==7){return(kSHP.test(X,Y))}
  if(t==8){return(BG(X,Y,R))}
}

#Example 1: Normal location Difference
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 5
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  m1 = rep(0,p)
  m2 = rep(0.15,p)
  var = diag(1,p)
  list = foreach(itr=1:ITR,.packages = c("vegan","cluster","MTSKNN"))%dopar%{
    set.seed(itr)
    X = mnormt::rmnorm(n,mean = m1,varcov = var)
    Y = mnormt::rmnorm(m,mean = m2,varcov = var)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)

#Example 2: Normal scale difference
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
R = 500
t = 5
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  var1 = diag(1,p)
  var2 = diag(1.1,p)
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = mnormt::rmnorm(n,varcov = var1)
    Y = mnormt::rmnorm(m,varcov = var2)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)

#Example 3: Normal distributional difference
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 2
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  c1 = rep(c(1,2), each = p/2)
  c2 = rep(c(2,1), each = p/2)
  var1 = diag(c1,p)
  var2 = diag(c2,p)
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = mnormt::rmnorm(n,varcov = var1)
    Y = mnormt::rmnorm(m,varcov = var2)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)


#Example 4: Cauchy location difference
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 2
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = matrix(rcauchy(n*p),ncol = p, nrow = n)
    Y = matrix(rcauchy(m*p,1,1),ncol = p, nrow = m)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)


#Example 6 : Normal equal mixture model
gen = function(n,d)
{
  v = rbinom(n,1,0.5)
  m = matrix(0,ncol = d,nrow = n)
  for(i in 1:n){
    if(v[i]==0){ m[i,] = rnorm(d,0.5,1)}else{
      m[i,] = rnorm(d,-0.5,1)
    }
  }
  return(m)
}
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 2
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  list = foreach(itr=1:ITR,.packages = c("vegan","cluster","MTSKNN"))%dopar%{
    set.seed(itr)
    X = matrix(rnorm(n*p),ncol = p, nrow = n)
    Y = gen(m,p)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)

#Example 7: Normal unequal mixture model
gen = function(n,d)
{
  v = rbinom(n,1,0.8)
  m = matrix(0,ncol = d,nrow = n)
  for(i in 1:n){
    if(v[i]==1){ m[i,] = rnorm(d,-0.25,1)}else{
      m[i,] = rnorm(d,1,1)
    }
  }
  return(m)
}
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 8
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  list = foreach(itr=1:ITR,.packages = c("vegan","cluster","MTSKNN"))%dopar%{
    set.seed(itr)
    X = matrix(rnorm(n*p),ncol = p, nrow = n)
    Y = gen(m,p)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)

#Example 8: Normal(0,2) v t_4 
d = 2^(1:10)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
n = 50
m = 50
ITR = 500
alpha = 0.05
t = 2
R = 500
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d)
{
  list = foreach(itr=1:ITR,.packages = c("vegan","cluster","MTSKNN"))%dopar%{
    set.seed(itr)
    X = matrix(rnorm(n*p,0, sd = sqrt(2)),ncol = p, nrow = n)
    Y = matrix(rt(m*p, df = 4), ncol = p, nrow = m)
    m = nrow(Y)
    n = nrow(X)
    d = ncol(X)
    N = n+m
    V = rbind(X,Y)
    D = matrix(0,N,N)
    for(i in 1:N){
      for(j in 1:N){
        l = 0
        for(k in 1:d){
          l = l+phi((V[i,k]-V[j,k])^2)
        }
        D[i,j] = l/d
      }
    }
    colnames(D) = 1:N
    rownames(D) = 1:N
    return(as.numeric(pval(X,Y,R,t, D=D)<alpha))
  }
  P = list[[1]]
  for(i in 2:ITR){P = c(P,list[[i]])}
  A = mean(P)
  print(A)
}

stopCluster(cl)

