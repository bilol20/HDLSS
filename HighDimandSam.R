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

#Example 5 and 9: Gaussian Shrinking Alternative
d = 2^(1:10)
gen1 = function(n,d,a)
{
  X = c(rnorm(d,mean = 1/(d^a),sd = 1))
  for(i in 2:n) 
    X= rbind(X,c(rnorm(d,mean = 1/(d^a),sd = 1)))
  return(X)
}

gen2 = function(n,d,a){
  X = c(rnorm(d,mean = -1/(d^a),sd = 1))
  for(i in 2:n) 
    X= rbind(X,c(rnorm(d,mean = -1/(d^a),sd = 1)))
  return(X)
}


library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)


ITR = 500
alpha = 0.05
R = 500
t = 2
delta = 0.5
gamma = 1.1
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d){
  n = 5 + as.integer(p^gamma)
  #n = 50
  m = n
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = gen1(n,p,delta)
    Y = gen2(m,p,delta)
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
  for(k in 2:ITR){P = c(P,list[[k]])}
  A = mean(P)
  print(c(delta,gamma,A))
}

stopCluster(cl)


power = matrix(rep(0), nrow = 8, ncol = 10)
colnames(power) = c("l2","l1","exp","log","FR","BF","knn","mmd","SHP","BG")

#Example 10: Normal distributional difference
d = 2^(1:8)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
ITR = 500
alpha = 0.05
t = 8
R = 500
gamma = 1
phi_x = function(x,k){
  if(k==1){return(x)}
  if(k==2){return(sqrt(x))}
  if(k==3){return(1-exp(-x/2))}
  if(k==4){return(log(1+x))}
}
for(index in 1:4){
  t=2
  q = c()
  for(p in d)
  {
    n = 5 + as.integer(p^gamma)
    m = n
    c1 = rep(c(1,5), each = p/2)
    c2 = rep(c(5,1), each = p/2)
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
            l = l+phi_x((V[i,k]-V[j,k])^2,index)
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
    q = c(q,A)
  }
  power[,index] = q
}
for(t in 3:8){
  q = c()
  for(p in d)
  {
    n = 5 + as.integer(p^gamma)
    m = n
    c1 = rep(c(1,5), each = p/2)
    c2 = rep(c(5,1), each = p/2)
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
            l = l+phi_x((V[i,k]-V[j,k])^2,1)
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
    q = c(q,A)
  }
  power[,t+2] = q
}

stopCluster(cl)


#Example 11
varf = function(d,a){
  A = matrix(1,d,d)
  for(i in 1:d){
    for(j in 1:d){
      A[i,j] = a^(abs(i-j))
    }
  }
  return(A)
}

d = 2^(1:8)
library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)

nl = length(d)
ITR = 500
alpha = 0.05
t = 2
R = 500
gamma = 1
phi_x = function(x,k){
  if(k==1){return(x)}
  if(k==2){return(sqrt(x))}
  if(k==3){return(1-exp(-x/2))}
  if(k==4){return(log(1+x))}
}
for(index in 1:4){
  t=2
  q = c()
  for(p in d)
  {
    n = 5 + as.integer(p^gamma)
    #n = 30
    m = n
    var1 = varf(p,0.5)
    var2 = varf(p,0.1)   
    list = foreach(itr=1:ITR,.packages = c("MTSKNN","cluster"))%dopar%{
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
            l = l+phi_x((V[i,k]-V[j,k])^2,index)
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
    q = c(q,A)
  }
  power[,index] = q
}
for(t in 3:8){
  q = c()
  for(p in d)
  {
    n = 5 + as.integer(p^gamma)
    #n = 50
    m = n
    var1 = varf(p,0.5)
    var2 = varf(p,0.1)   
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
            l = l+phi_x((V[i,k]-V[j,k])^2,1)
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
    q = c(q,A)
  }
  power[,t+2] = q
}

stopCluster(cl)



#Example 12: Few coordinate difference normal mean change
d = 2^(1:8)
gen = function(n,d,a)
{
  l = as.integer(d^a)
  x = c()
  for(i in 1:l){
    x = c(x,rnorm(1,mean = 2, sd = 1))
  }
  x = c(x, rnorm(d-l))
  X = x
  for(i in 2:n){
    x = c()
    for(i in 1:l){
      x = c(x,rnorm(1,mean = 2, sd = 1))
    }
    x = c(x, rnorm(d-l))
    X= rbind(X,x)
  } 
  return(X)
}


library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)


ITR = 500
alpha = 0.05
R = 500
t = 5
delta = 0.7
gamma = 0.5
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d){
  n = 5 + as.integer(p^gamma)
  m = n
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = gen(n,p,delta)
    Y = matrix(rnorm(m*p),ncol = p, nrow = m)
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
  for(k in 2:ITR){P = c(P,list[[k]])}
  A = mean(P)
  print(c(n,delta,gamma,A))
}

stopCluster(cl)

#Example 13: Few coordinate difference normal scale change
d = 2^(1:8)
gen = function(n,d,a)
{
  l = as.integer(d^a)
  x = c()
  for(i in 1:l){
    x = c(x,rnorm(1,mean = 0, sd = sqrt(5)))
  }
  x = c(x, rnorm(d-l))
  X = x
  for(i in 2:n){
    x = c()
    for(i in 1:l){
      x = c(x,rnorm(1,mean = 0, sd = sqrt(5)))
    }
    x = c(x, rnorm(d-l))
    X= rbind(X,x)
  } 
  return(X)
}


library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)


ITR = 500
alpha = 0.05
R = 500
t = 8
delta = 0.7
gamma = 0.5
phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
#phi = function(x) return(log(1+x))
for(p in d){
  n = 5+as.integer(p^gamma)
  m = 5+as.integer(p^gamma)
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = gen(n,p,delta)
    Y = matrix(rnorm(m*p),ncol = p, nrow = m)
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
  for(k in 2:ITR){P = c(P,list[[k]])}
  A = mean(P)
  print(c(delta,gamma,A))
}

stopCluster(cl)

#Example 14: Few coordinate difference in distribution
d = 2^(1:8)
gen = function(n,d,a)
{
  l = as.integer(d^a)
  x = c()
  for(i in 1:l){
    x = c(x,rt(1,df = 4))
  }
  x = c(x, rnorm(d-l,mean = 0,sd = sqrt(2)))
  X = x
  for(i in 2:n){
    x = c()
    for(i in 1:l){
      x = c(x,rt(1,df = 4))
    }
    x = c(x, rnorm(d-l,mean = 0,sd = sqrt(2)))
    X= rbind(X,x)
  } 
  return(X)
}


library(doParallel)
cl=makeCluster(min(120,detectCores()))
registerDoParallel(cl)


ITR = 500
alpha = 0.05
R = 500
t = 8
delta = 0.7
gamma = 0.5
#phi = function(x) return(x)
#phi = function(x) return(sqrt(x))
#phi = function(x) return(1-exp(-x/2))
phi = function(x) return(log(1+x))
for(p in d){
  n = 5+as.integer(p^gamma)
  m = 5+as.integer(p^gamma)
  list = foreach(itr=1:ITR,.packages = c("MTSKNN","vegan","cluster"))%dopar%{
    set.seed(itr)
    X = gen(n,p,delta)
    Y = matrix(rnorm(m*p,mean = 0, sd = sqrt(2)),ncol = p, nrow = m)
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
  for(k in 2:ITR){P = c(P,list[[k]])}
  A = mean(P)
  print(c(delta,gamma,A))
}

stopCluster(cl)
