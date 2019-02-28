#===================================== Functions =========================================

# ======== Likelihood for Partial Recall MLE =======
fNEW<-function (param) {
  theta<-param[1];eta<-param[2]
  alpha1<-param[3];beta1<-param[4]
  alpha2<-param[5];beta2<-param[6]
  alpha3<-param[7];beta3<-param[8]
  ab<-theta/eta
  lnh<-which(data[,3]==4)
  li4<-log(pweibull(data[lnh,2], theta,eta,lower.tail=FALSE))
  rem<-which(data[,3]==0)
  pi0<-1/(1+exp(alpha1+beta1*(data[rem,2]-data[rem,1]))+exp(alpha2+beta2*(data[rem,2]-data[rem,1]))+exp(alpha3+beta3*(data[rem,2]-data[rem,1])))
  li0<-log(dweibull(data[rem,1], theta, eta)*pi0)
  
  NOR<-which(data[,3]==1)
  ll1<-length(NOR)
  s_vec1<-data[NOR,2]
  lens<-length(s_vec1)
  integ_d1 <- function(x) {
    i1<-ab*((x/eta)^(theta-1))*exp(-(x/eta)^theta)
    d1<-exp(alpha1+beta1*(s-x))
    d2<-exp(alpha2+beta2*(s-x))
    d3<-exp(alpha3+beta3*(s-x))
    inte<-i1*(d1/(1+d1+d2+d3))
    return(inte)
  }
  
  integd1<-array(0,lens)
  for (i in 1:lens) {
    s=s_vec1[i]
    integd1[i]<-integrate(integ_d1, 0, s,subdivisions=1000,stop.on.error = T)$value
  }
  li1<-log(integd1)
  
  MR<-which(data[,3]==2)
  ll2<-length(MR)
  s_vec2<-data[MR,2]
  L1M<-data[MR,4]
  U1M<-data[MR,5]
  lens2<-length(s_vec2)
  integ_d2 <- function(x) {
    i1<-ab*((x/eta)^(theta-1))*exp(-(x/eta)^theta)
    d1<-exp(alpha1+beta1*(s-x))
    d2<-exp(alpha2+beta2*(s-x))
    d3<-exp(alpha3+beta3*(s-x))
    inte<-i1*(d2/(1+d1+d2+d3))
    return(inte)
  }
  integd2<-array(0,lens2)
  for (i in 1:lens2) {
    s=s_vec2[i]
    l=L1M[i]
    u<-U1M[i]
    integd2[i]<-integrate(integ_d2, l, u,subdivisions=1000,stop.on.error = T)$value
  }
  li2<-log(integd2)
  
  YR<-which(data[,3]==3)
  ll3<-length(YR)
  s_vec3<-data[YR,2]
  L1Y<-data[YR,4]
  U1Y<-data[YR,5]
  lens3<-length(s_vec3)
  integ_d3 <- function(x) {
    i1<-ab*((x/eta)^(theta-1))*exp(-(x/eta)^theta)
    d1<-exp(alpha1+beta1*(s-x))
    d2<-exp(alpha2+beta2*(s-x))
    d3<-exp(alpha3+beta3*(s-x))
    inte<-i1*(d3/(1+d1+d2+d3))
    return(inte)
  }
  integd3<-array(0,lens3)
  for (i in 1:lens3) {
    s=s_vec3[i]
    l<-L1Y[i]
    u<-U1Y[i]
    integd3[i]<-integrate(integ_d3, l, u,subdivisions=1000,stop.on.error = T)$value
  }
  li3<-log(integd3)
  
  liken<-c(li0,li1,li2,li3,li4, recursive=TRUE)
  likeNEW<-sum(liken)    
  return(- likeNEW)
}
# ======== Likelihood for Binary Recall MLE =======
fNEW_EXNR<-function (param) {
  theta<-param[1];eta<-param[2]
  alpha<-param[3];beta<-param[4]
  ab<-theta/eta
  lnh<-which(datam[,3]==4)
  li4<-log(pweibull(datam[lnh,2], theta,eta,lower.tail=FALSE))
  rem<-which(datam[,3]==0)
  pi0<-1/(1+exp(alpha+beta*(datam[rem,2]-datam[rem,1])))
  li0<-log(dweibull(datam[rem,1], theta, eta)*pi0)
  
  NOR<-which(datam[,3]==1)
  ll1<-length(NOR)
  s_vec1<-datam[NOR,2]
  lens<-length(s_vec1)
  integ_d1 <- function(x) {
    i1<-ab*((x/eta)^(theta-1))*exp(-(x/eta)^theta)
    d1<-exp(alpha+beta*(s-x))
    inte<-i1*(d1/(1+d1))
    return(inte)
  }
  integd1<-array(0,lens)
  for (i in 1:lens) {
    s=s_vec1[i]
    integd1[i]<-integrate(integ_d1, 0, s,subdivisions=1000,stop.on.error = T)$value
  }
  li1<-log(integd1)
  liken<-c(li0,li1,li4, recursive=TRUE)
  likeNEW_EXN<-sum(liken)    
  return(- likeNEW_EXN)
}


#========== Likelihood function for Current Status MLE ==============
fYNO<-function (ss)   {
  if (ss[1]>0 & ss[2]>0)   {
    alpha<-ss[1];beta<-ss[2]
    lyes<-which(data[,3]!=4)
    lno<-which(data[,3]==4)
    lyn1<-log(pweibull(data[lno,2], alpha,beta,lower.tail=FALSE))
    lyn2<-log(pweibull(data[lyes,2], alpha, beta))
    likeYN<-c(lyn1,lyn2, recursive=TRUE)
    likeyn<-sum(likeYN) 
  }
  else   likeyn=-999999
  #  print(c(-likeyn,ss))
  return(- likeyn)
}
