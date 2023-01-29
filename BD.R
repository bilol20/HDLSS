St.2 = function(R,size){
  n = size[1]
  m = size[2]
  N = sum(size)
  l1 = 1:n
  l2 = 1:m
  L1 = matrix(0, ncol = n, nrow = n)
  L2 = matrix(0, ncol = n, nrow = n)
  L3 = matrix(0, ncol = m, nrow = m)
  L4 = matrix(0, ncol = m, nrow = m)
  for(i in l1){
    lab = as.integer(names(R[[i]]))
    iind = which(lab == i)
    laby = which(lab>n)
    labx = which(lab<=n)
    for(j in l1[-i]){
      l = lab[-iind]
      l = which(l == j)
      R1 = .subset2(R,i)[-c(laby,iind)]
      L1[i,j] = which(as.integer(names(R1))==j)-1
      L2[i,j] = l-L1[i,j]-1
    }
  }
  A = sum((L1/(n-2)-L2/m)^2)/(n*(n-1))
  for(i in l2){
    lab = as.integer(names(R[[i+n]]))
    iind = which(lab == i+n)
    labx = which(lab<=n)
    laby = which(lab>n)
    for(j in l2[-i]){ 
      l = lab[-iind]
      l = which(l == j+n)
      R1 = .subset2(R,i+n)[-c(labx,iind)]
      L3[i,j] = which(as.integer(names(R1))==j+n)-1
      L4[i,j] = l-L3[i,j]-1
    }
  }
  C = sum((L3/(m-2)-L4/(n))^2)/(m*(m-1))
  return(A+C)
}


bd_test_gdist = function(D,P = 200){
  m = nrow(Y)
  n = nrow(X)
  d = ncol(X)
  N = n+m
  R = lapply(1:N, function(x){
    return(sort(D[x,]))
  })
  T = St.2(R,c(n,m))
  T1 = lapply(1:P,function(k){
    l = sample(N)
    l1 = which(l %in% 1:N)
    R1 = lapply(1:N, function(k){
      x = R[[l[k]]]
      q = as.numeric(names(x))
      names(x) = match(q,l)
      return(x)
    })
    return(St.2(R1, c(n,m)))
  })
  (sum(T1>T)+1)/(P+1)
}
