library(Rsolnp)
mydata <- read.csv("Data1.csv")           
##############Data Format###############
#Column 1: Exact Age (Years)
#Column 2: mens onset days min (Days)
#Column 3: mens onset days max (Days)
#Column 4: Recall Type code
#Column 5: Recall Type
###################################################################################
names(mydata)
mydata[, 2] <- (mydata[, 2]) / 365.25
mydata[, 3] <- (mydata[, 3]) / 365.25
fact <- factor(mydata[, 5])
S <- mydata[, 1]
T1 <- mydata[, 2]
T2 <- mydata[, 3]
ind <- mydata[, 5]
d <- ind
v <- factor(ind)
levels(v)
k = 4
Xvec <- c(0, 3, 6, 9, 12)
delta <- ifelse(d == "Not happened", 0, 1)
NHp <- (length(d[d == "Not happened"])) / (length(d))
NOp <- (length(d[d == "No recall"])) / (length(d))
EXp <- (length(d[d == "Day recall"])) / (length(d))
MOp <- (length(d[d == "Month recall"])) / (length(d))
YRp <- (length(d[d == "Year recall"])) / (length(d))
n2 <- length(d[d == "Day recall"])
Tex <- T1[d == "Day recall"]
T_mass <- sort(Tex, decreasing = FALSE)
n <- length(S)
r <- 12 + n2
ro <- 4 + n2
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
}
  
  return(-loglike)
}
###########Setting initial Values##########
u <- runif(n2, 0, 1)
qq <- u / sum(u)
b00 <- c(.08, .39, .49, .62)
b22 <- c(.18, .25, .14, .13)
b33 <- c(.09, .18, .16, .17)
b11 <- (1 - (b00 + b22 + b33))
lamda <- c(b00, b22, b33, qq)
lamda1 <- c(b00, qq)
r <- 12 + n2
ro <- 4 + n2

##########Setting constraints###############
equal <- function(lamda) {
  sum(lamda[13:r])
}

inequal <- function(lamda) {
  s1 <- lamda[1] + lamda[5] + lamda[9]
  s2 <- lamda[2] + lamda[6] + lamda[10]
  s3 <- lamda[3] + lamda[7] + lamda[11]
  s4 <- lamda[4] + lamda[8] + lamda[12]
  c(s1, s2, s3, s4)
}


func1 <- function(lamda)
{
  AMLEneglog.likelihood(lamda, S, T1, T2, d, Xvec)
}
equal1 <- function(lamda1) {
  sum(lamda1[5:ro])
}
func2 <- function(lamda1)
{
  AMLEExneglog.likelihood(lamda1, S, T1, T2, d, Xvec)
}


#########Performing Optimization###############
m <- solnp(
  lamda,
  func1,
  eqfun = equal,
  eqB = 1,
  LB = rep(0, r),
  UB = rep(1, r)
)
n <- solnp(
  lamda1,
  func2,
  eqfun = equal1,
  eqB = 1,
  LB = rep(0, ro),
  UB = rep(1, ro)
)


qopt <- m$par[13:r]
qopt1 <- n$par[5:ro]
Tkaplan <- T1[delta == 1]
Tkap <- T1

Tmat <- T_mass



#######Calculating the NP distribution function#################
FN <- function(t)
  
{
  temp <- c()
  for (j in 1:n2)
  {
    if (T_mass[j] <= t)
    {
      temp[j] = qopt[j]
    }
    
    if (T_mass[j] > t)
    {
      temp[j] = 0
    }
  }
  sum(temp)
}
FN1 <- Vectorize(FN)

FN0 <- function(t)
  
{
  temp <- c()
  for (j in 1:n2)
  {
    if (T_mass[j] <= t)
    {
      temp[j] = qopt1[j]
    }
    
    if (T_mass[j] > t)
    {
      temp[j] = 0
    }
  }
  sum(temp)
}
FN2 <- Vectorize(FN0)



###================ For NP distribution function plot ===========
qdata <- read.csv(file.choose())
T_mass <- qdata$T_mass
qopt <- qdata$qopt
FN <- function(t)
  
{
  temp <- c()
  for (j in 1:length(T_mass))
  {
    if (T_mass[j] <= t)
    {
      temp[j] = qopt[j]
    }
    
    if (T_mass[j] > t)
    {
      temp[j] = 0
    }
  }
  sum(temp)
}
FN1 <- Vectorize(FN)

qopt1 <- qdata$qopt1

FN0 <- function(t)
  
{
  temp <- c()
  for (j in 1:length(T_mass))
  {
    if (T_mass[j] <= t)
    {
      temp[j] = qopt1[j]
    }
    
    if (T_mass[j] > t)
    {
      temp[j] = 0
    }
  }
  sum(temp)
}
FN2 <- Vectorize(FN0)
t <- seq(4, 16, l = 100)
plot(t, FN1(t), type = "l", col = "green")
lines(t, FN2(t), type = "l", col = "red")
#################################################################
