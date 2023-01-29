library(HDLSSkST) # our proposed tests
library(energy)  # for Energy
library(LPKsample) # for GLP
library(MTSKNN) # for KNN
library(MASS)



#-------------------REZA(2020)--------------------------------------
statfun <- function(n1,n2,X){
  I1 <- X[1:n1, ]
  I2 <- X[(n1+1):(n1+n2), ]
  
  m1=n1*(n1-1)/2
  m2=n2*(n2-1)/2
  
  m12=n1*n2
  
  sample1=t(I1)
  sample2=t(I2)
  
  
  count=1
  ipd1=rep(0,m1)
  for(i in 1:(n1-1)){
    for(j in (i+1):n1){
      ipd1[count]=sqrt((t(sample1[,i]-sample1[,j]) %*% (sample1[,i]-sample1[,j])))
      count <- count+1
    }
  }
  
  count=1
  ipd2=rep(0,m2)
  for(i in 1:(n2-1)){
    for(j in (i+1):n2){
      ipd2[count]=sqrt((t(sample2[,i]-sample2[,j]) %*% (sample2[,i]-sample2[,j])))
      count <- count+1
    }
  }
  
  
  count=1
  ipd12=rep(0,m12)
  for(i in 1:n1){
    for(j in 1:n2){
      ipd12[count]=sqrt((t(sample1[,i]-sample2[,j]) %*% (sample1[,i]-sample2[,j])))
      count <- count+1
    }
  }
  
  N = m1+m2
  ipd1bar=sum(ipd1)/m1
  ipd2bar=sum(ipd2)/m2
  
  ipd12bar=sum(ipd12)/m12
  
  
  dbar = (m1*ipd1bar + m2*ipd2bar)/N
  
  TSS = (ipd12bar-dbar)^2
  
  MD = (ipd1bar - ipd12bar)^2 + (ipd2bar -ipd12bar)^2 
  
  MV = (ipd1bar - ipd12bar)^2 + (ipd2bar -ipd12bar)^2 
  stas = c(TSS,MD,MV)
  return(stas)
}

reza.test.TSS = function(X,Y, R = 100)
{
  n1 = nrow(X)
  n2 = nrow(Y)
  V = rbind(X,Y)
  T = statfun(n1,n2,V)
  T1 = sapply(1:R, function(k){
    l = sample(n+m)
    V1 = V[l,]
    return(statfun(n1,n2,V1))
  })
  return((sum(T1[1,]>T[1])+1)/(R+1))
}

BG = function(X,Y, R = 100)
{
  n1 = nrow(X)
  n2 = nrow(Y)
  V = rbind(X,Y)
  T = statfun(n1,n2,V)
  T1 = sapply(1:R, function(k){
    l = sample(n1+n2)
    V1 = V[l,]
    return(statfun(n1,n2,V1))
  })
  p.MD = (sum(T1[2,]>T[2])+1)/(R+1)
  return(p.MD)
}

#----------------------------------------------------------------------------

#------------------------K-NN by Schilling(1986)-----------------------------------------------
mksknn <- function(X,Y, k_n){
  M = rbind(X,Y)
  n1 = nrow(X)
  n2 = nrow(Y)
  labels = c(rep(0,n1),rep(1,n2))
  d <- ncol(M)
  n <- nrow(M)
  Set <- cbind(M,labels)
  tSet <- t(Set)
  
  output <- rep(0,n)   
  C.out <- .C("knn", as.double(tSet), as.integer(n), as.integer(d), as.integer(k_n), as.integer(output))
  counts <- C.out[[5]]
  Tk <- sum(counts)/(n*k_n)
  
  return(Tk) # test statistic
}

knnTest = function(X,Y, k, R = 200){
  V = rbind(X,Y)
  n = nrow(X)
  m = nrow(Y)
  T = mksknn(X,Y,k)
  T1 = lapply(1:R, function(i){
    l = sample(n+m)
    X1 = V[l[1:n],]
    Y1 = V[l[n+1:m],]
    return(mksknn(X1,Y1,k))
  })
  (sum(T1>T)+1)/(R+1)
}

#----------------------------------------------------------------------------

#-----------------Hamiltonian Path based method by Biswas -------------------
pr = function(k,n,m){
  l = k%%2
  if(l==0){
    k1 = k/2
    return(2*choose(m-1,k1-1)*choose(n-1,k1-1)/choose(n+m,n))}else{
      k1 = (k+1)/2  
      return((choose(m-1,k1-1)*choose(n-1,k1-2)+choose(m-1,k1-2)*choose(n-1,k1-1))/
               choose(n+m,n))
    }
}

Pr = function(k,n,m)
{
  l = 2:k
  s = 0
  for(i in l){ s = s+pr(i,n,m)}
  return(s)
}


#Testing
kSHP <- function(M, labels){
  n <- nrow(M)
  D <- dist(M)
  newa1 <- c()
  newa2 <- c()
  for (i in 1:(n-1)){
    newa1 <- c(newa1,rep(i,n-i))
    newa2 <- c(newa2,(i+1):n)
  }
  newa <- cbind(newa1,newa2)
  vector <- as.numeric(D)
  MM <- n*(n-1)/2
  v.sort <- sort.int(vector,index.return=TRUE)
  COST <- v.sort$x
  AA <- v.sort$ix
  deg <- rep(0,n)
  valid <- 0 
  i <- 1
  Link <- matrix(0,0,2)
  while ( i<=MM){ 
    index = AA[i]
    link1 = newa[index,1]
    link2 = newa[index,2]
    flag = 0
    if( (deg[link1]==2)|(deg[link2]==2)){
      flag = 1
    }
    if (flag==0){
      deg[link1]=deg[link1]+1
      deg[link2]=deg[link2]+1
      valid=valid+1
      Link <- rbind(Link,c(link1,link2))
    }
    
    i=i+1
    if (valid>=(n-1)){
      i=MM+1
    }
    
    if (deg[link1]==1 & deg[link2]==1) {
      search=link1
      
      label<-rep(0,max(1,valid))
      
      signal=0
      while(signal==0){
        for(q in 1:max(1,valid)){
          if (label[q]==0){
            if(Link[q,1]==search)
              label[q]=1
            newsearch=Link[q,2]
            end
            if(Link[q,2]==search)
              label[q]=1
            newsearch=Link[q,1]
            end
            if ((Link[q,1]==search) | (Link[q,2]==search))
              q=valid+1
            end
          }
        }
        if (deg[newsearch]==1){ 
          signal=1
        }
        if (newsearch==link2){
          signal=2
        }
        search=newsearch
      }
      if (signal==2){
        flag=1
      }
    }
  }
  runstat=1
  for (i in 1:(n-1)){
    if (labels[Link[i,1]]!=labels[Link[i,2]]){ 
      runstat=runstat+1
    }
  }
  return(runstat)
}
kSHP.test = function(X,Y)
{
  n1 = nrow(X)
  n2 = nrow(Y)
  V = rbind(X,Y)
  labels = c(rep(1,n1),rep(2,n2))
  T = kSHP(V,labels)
  return(Pr(T,n1,n2))
}
#-----------------Baringhaus and Franz-------------------
norm = function(x) sqrt(sum(x^2))
BFTest = function(X,Y,R = 200){
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  M1 = matrix(0,ncol = nrow(V), nrow = nrow(V))
  for(i in 1:nrow(V)){
    for(j in 1:nrow(V)){
      M1[i,j] = norm(V[i,]-V[j,])
    }
  }
  
  T = (n*m/(n+m))*(2*mean(M1[1:n,n+1:m])-mean(M1[1:n,1:n])-mean(M1[n+1:m,n+1:m]))
  
  T1 = unlist(lapply(1:R, function(k){
    l = sample(n+m)
    M2 = M1[l,l]
    return((n*m/(n+m))*(2*mean(M2[1:n,n+1:m])-mean(M2[1:n,1:n])-mean(M2[n+1:m,n+1:m])))
  }))
  return((sum(T1>T)+1)/(R+1))
}

#----------------BG--------------------
norm = function(x) sqrt(sum(x^2))
BGTest = function(X,Y,R = 200){
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  M1 = matrix(0,ncol = nrow(V), nrow = nrow(V))
  for(i in 1:nrow(V)){
    for(j in 1:nrow(V)){
      M1[i,j] = norm(V[i,]-V[j,])
    }
  }
  
  T = (sum(M1[1:n,n+1:m])/n/m-sum(M1[1:n,1:n])/n/(n-1))^2
  +  (sum(M1[1:n,n+1:m])/n/m-sum(M1[n+1:m,n+1:m])/m/(m-1))^2
  
  T1 = unlist(lapply(1:R, function(k){
    l = sample(n+m)
    M2 = M1[l,l]
    return((sum(M2[1:n,n+1:m])/n/m-sum(M2[1:n,1:n])/n/(n-1))^2
           +  (sum(M2[1:n,n+1:m])/n/m-sum(M2[n+1:m,n+1:m])/m/(m-1))^2
    )
  }))
  return((sum(T1>T)+1)/(R+1))
}

#-----------------MMD-------------------
norm = function(x) sqrt(sum(x^2))
mmd = function(X,Y,R = 200)
{
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  M = matrix(0,n+m,n+m)
  for(i in 1:(n+m)){
    for(j in 1:(n+m)){
      M[i,j] = norm(V[i,]-V[j,])
    }
  }
  s = quantile(as.numeric(M[upper.tri(M)]^2),0.5)  
  K = exp(-M^2/(s))
  T =(-2*mean(K[1:n,n+1:m])+mean(K[1:n,1:n])+mean(K[n+1:m,n+1:m]))
  
  T1 = unlist(lapply(1:R, function(k){
    l = sample(n+m)
    M2 = K[l,l]
    return((-2*mean(M2[1:n,n+1:m])+mean(M2[1:n,1:n])+mean(M2[n+1:m,n+1:m])))
  }))
  return((sum(T1>T)+1)/(R+1))
}


#-----------------Friedman and Rafzsky-------------------
library(vegan)
library(cluster)
FRTest = function(X,Y,R=100)
{
  m = nrow(X)
  n = nrow(Y)
  V = rbind(X,Y)
  dist = as.matrix(vegdist(V, binary=FALSE, method="euclidean"))
  child = spantree(dist)$kid
  edges = t(rbind(c(2:(m+n)), child))
  
  T = length(which(edges[1:(m-1),][,2]>m))+
    length(which(edges[m:(m+n-1),][,2]<=m))
  T1 = unlist(lapply(1:R, function(k){
    l = sample(n+m)
    V1 = V[l,]
    dist = as.matrix(vegdist(V1, binary=FALSE, method="euclidean"))
    child = spantree(dist)$kid
    edges = t(rbind(c(2:(m+n)), child))
    return(length(which(edges[1:(m-1),][,2]>m))+
             length(which(edges[m:(m+n-1),][,2]<=m)))
  }))
  return((sum(T1<T)+1)/(R+1))
}
