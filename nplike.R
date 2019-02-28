####################Nonparametric Likelihood for Partial Recall#####################
AMLEneglog.likelihood <- function(lamda, S, T1, T2, d, Xvec)
{
  b0 <- lamda[1:4]
  b2 <- lamda[5:8]
  b3 <- lamda[9:12]
  b1 <- (1 - (b0 + b2 + b3))
  if(min(b1)< 0)
  {loglike = -999999}
  
  if(min(b1)>0){
  n2 <- length(d[d == "Day recall"])
  r <- (12 + n2)
  q <- lamda[13:r]
  Tex <- T1[d == "Day recall"]
  
  n <- length(S)
  T_mass <- sort(Tex, decreasing = FALSE)
  a <- mat.or.vec(nr = n, nc = n2)
  for (i in 1:n)
  {
    for (j in 1:n2)
    {
      if (d[i] == "Not happened")
      {
        a[i, j] = ifelse(T_mass[j] > S[i], 1, 0)
      }
      
      if (d[i] == "Day recall")
      {
        k <- length(b1)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b1[l] * ifelse(T_mass[j] == T1[i], 1, 0) * ifelse(T1[i] > (S[i] -
                                                                                 Xvec[l + 1]) && T1[i] <= (S[i] - Xvec[l]), 1, 0)
        }
        conl[k] = b1[k] * ifelse(T_mass[j] == T1[i], 1, 0) * ifelse(T1[i] <=
                                                                      (S[i] - Xvec[k]), 1, 0)
        a[i, j] = sum(conl)
      }
      
      if (d[i] == "No recall")
      {
        k <- length(b0)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b0[l] * ifelse(T_mass[j] > (S[i] - Xvec[l + 1]) &&
                                     T_mass[j] <= (S[i] - Xvec[l]), 1, 0)
        }
        conl[k] = b0[k] * ifelse(T_mass[j] <= (S[i] - Xvec[k]), 1, 0)
        a[i, j] = sum(conl)
      }
      
      
      
      if (d[i] == "Month recall")
      {
        k <- length(b2)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b2[l] * ifelse(T_mass[j] > max((S[i] - Xvec[l + 1]), (T1[i])) &&
                                     T_mass[j] <= min((S[i] - Xvec[l]), (T2[i])), 1, 0)
        }
        conl[k] = b2[k] * ifelse(T_mass[j] <= min((S[i] - Xvec[k]), (T2[i])) &&
                                   T_mass[j] > (T1[i]) , 1, 0)
        a[i, j] = sum(conl)
      }
      
      if (d[i] == "Year recall")
      {
        k <- length(b3)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b3[l] * ifelse(T_mass[j] > max((S[i] - Xvec[l + 1]), (T1[i])) &&
                                     T_mass[j] <= min((S[i] - Xvec[l]), (T2[i])), 1, 0)
        }
        conl[k] = b3[k] * ifelse(T_mass[j] <= min((S[i] - Xvec[k]), (T2[i])) &&
                                   T_mass[j] > (T1[i]) , 1, 0)
        a[i, j] = sum(conl)
      }
    }
  }
  
  
  g <- c()
  for (i in 1:n)
  {
    p <- c()
    for (j in 1:n2)
    {
      p[j] <- a[i, j] * q[j]
    }
    g[i] = sum(p)
  }
  
  l <- c()
  for (i in 1:n)
  {
    if (!is.na(g[i]) && g[i] != 0)
    {
      l[i] <- log(g[i])
    } else{
      l[i] <- 0
    }
    
  }
  
  l <- l[is.finite(l)]
  
  loglike <- sum(l)
  
  if (any(is.finite(loglike)))
  {
    loglike = loglike
  } else{
    loglike <- -999999
  }
  
  }
  return(-loglike)
}

####################Nonparametric Likelihood for Partial Recall#####################
AMLEExneglog.likelihood <- function(lamda1, S, T1, T2, d, Xvec)
{
  b0 <- lamda1[1:4]
  b1 <- (1 - b0)
  n2 <- length(d[d == "Day recall"])
  ro <- (4 + n2)
  q <- lamda[5:ro]
  Tex <- T1[d == "Day recall"]
  
  n <- length(S)
  T_mass <- sort(Tex, decreasing = FALSE)
  a <- mat.or.vec(nr = n, nc = n2)
  for (i in 1:n)
  {
    for (j in 1:n2)
    {
      if (d[i] == "Not happened")
      {
        a[i, j] = ifelse(T_mass[j] > S[i], 1, 0)
      }
      
      if (d[i] == "Day recall")
      {
        k <- length(b1)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b1[l] * ifelse(T_mass[j] == T1[i], 1, 0) * ifelse(T1[i] > (S[i] -
                                                                                 Xvec[l + 1]) && T1[i] <= (S[i] - Xvec[l]), 1, 0)
        }
        conl[k] = b1[k] * ifelse(T_mass[j] == T1[i], 1, 0) * ifelse(T1[i] <=
                                                                      (S[i] - Xvec[k]), 1, 0)
        a[i, j] = sum(conl)
      }
      
      if (d[i] == "No recall" ||
          d[i] == "Month recall" || d[i] == "Year recall")
      {
        k <- length(b0)
        conl <- c()
        for (l in 1:k - 1)
        {
          conl[l] = b0[l] * ifelse(T_mass[j] > (S[i] - Xvec[l + 1]) &&
                                     T_mass[j] <= (S[i] - Xvec[l]), 1, 0)
        }
        conl[k] = b0[k] * ifelse(T_mass[j] <= (S[i] - Xvec[k]), 1, 0)
        a[i, j] = sum(conl)
      }
      
      
    }
  }
  
  
  g <- c()
  for (i in 1:n)
  {
    p <- c()
    for (j in 1:n2)
    {
      p[j] <- a[i, j] * q[j]
    }
    g[i] = sum(p)
  }
  
  l <- c()
  for (i in 1:n)
  {
    if (!is.na(g[i]) && g[i] != 0)
    {
      l[i] <- log(g[i])
    } else{
      l[i] <- 0
    }
    
  }
  l <- l[is.finite(l)]
  loglike <- sum(l)
  
  if (any(is.finite(loglike)))
  {
    loglike = loglike
  } else{
    loglike <- -999999
  }
  
  return(-loglike)
}