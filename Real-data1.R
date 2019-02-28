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

#=================== Main Program =================   
mydata<- read.csv("Data1.csv",header=T)
names(mydata)
dimex<-dim(mydata)
group<-mydata[,4]

aam<-numeric(dimex[1])
L1<-numeric(dimex[1])
U1<-numeric(dimex[1])
ind0<-which(group==0)
L1<-mydata[,2]/365.25
U1<-mydata[,3]/365.25

aam[ind0]<-mydata[ind0,2]/365.25
aai<-mydata[,1]
data<-cbind(aam,aai,group,L1,U1)


leEX<-length(which(data[,3]==0))
leNO<-length(which(data[,3]==1))
leM<-length(which(data[,3]==2))
leY<-length(which(data[,3]==3))
leNH<-length(which(data[,3]==4))
leEX;leNO;leM;leY;leNH
sall<-leEX+leNO+leM+leY+leNH

#================= For Binary recall data, we consider monthly and yearly recall as no recall ============
#==== Data set for Binary Recal MLE============
datam<-data
da2<-which(datam[,3]==2)
da3<-which(datam[,3]==3)
datam[da2,3]<-1
datam[da3,3]<-1

#=== For those individual who had the event happend, Calculate S-T ===
STda<-mydata[which(mydata[,4]!=4),]
TAVG<-(STda[,3]+STda[,2])/2
STdat<-cbind(TAVG/365.25,STda)
STdata<-cbind(STdat[,1],STdat[,2],STdat[,5],STdat[,2]-STdat[,1])
colnames(STdata) <- c("aam","aai","group","S-T")

#============== Choosing the initial values of parameter for optimization algorithm =============
## partition for pi0,pi1,pi2,pi3
g1<-which(STdata[,4]<=2)
ST1<-STdata[g1,]
p01<-length(which(ST1[,3]==0))/dim(ST1)[1]
p11<-length(which(ST1[,3]==1))/dim(ST1)[1]
p21<-length(which(ST1[,3]==2))/dim(ST1)[1]
p31<-length(which(ST1[,3]==3))/dim(ST1)[1]

g2<-which(STdata[,4]>2 & STdata[,4]<=3.5)
ST2<-STdata[g2,]
p02<-length(which(ST2[,3]==0))/dim(ST2)[1]
p12<-length(which(ST2[,3]==1))/dim(ST2)[1]
p22<-length(which(ST2[,3]==2))/dim(ST2)[1]
p32<-length(which(ST2[,3]==3))/dim(ST2)[1]

g3<-which(STdata[,4]>3.5 & STdata[,4]<=5)
ST3<-STdata[g3,]
p03<-length(which(ST3[,3]==0))/dim(ST3)[1]
p13<-length(which(ST3[,3]==1))/dim(ST3)[1]
p23<-length(which(ST3[,3]==2))/dim(ST3)[1]
p33<-length(which(ST3[,3]==3))/dim(ST3)[1]

g4<-which(STdata[,4]>5 & STdata[,4]<=6.5)
ST4<-STdata[g4,]
p04<-length(which(ST4[,3]==0))/dim(ST4)[1]
p14<-length(which(ST4[,3]==1))/dim(ST4)[1]
p24<-length(which(ST4[,3]==2))/dim(ST4)[1]
p34<-length(which(ST4[,3]==3))/dim(ST4)[1]

g5<-which(STdata[,4]>6.5 & STdata[,4]<=8)
ST5<-STdata[g5,]
p05<-length(which(ST5[,3]==0))/dim(ST5)[1]
p15<-length(which(ST5[,3]==1))/dim(ST5)[1]
p25<-length(which(ST5[,3]==2))/dim(ST5)[1]
p35<-length(which(ST5[,3]==3))/dim(ST5)[1]

g6<-which(STdata[,4]>8 & STdata[,4]<=12.3)
ST6<-STdata[g6,]
p06<-length(which(ST6[,3]==0))/dim(ST6)[1]
p16<-length(which(ST6[,3]==1))/dim(ST6)[1]
p26<-length(which(ST6[,3]==2))/dim(ST6)[1]
p36<-length(which(ST6[,3]==3))/dim(ST6)[1]

# define x and y to do regression of log pi/pi0=alpha+beta*x
x<-c(1,2.75,4.25,5.75,7.25,10.15)
p11<-0.01
y1<-c(log(p11/p01),log(p12/p02),log(p13/p03),log(p14/p04),log(p15/p05),log(p16/p06))
out1<-lm(y1~x)
alpha1_ini<-summary(out1)$coeff[[1]]
beta1_ini<-summary(out1)$coeff[[2]]

y2<-c(log(p21/p01),log(p22/p02),log(p23/p03),log(p24/p04),log(p25/p05),log(p26/p06))
out2<-lm(y2~x)
alpha2_ini<-summary(out2)$coeff[[1]]
beta2_ini<-summary(out2)$coeff[[2]]

y3<-c(log(p31/p01),log(p32/p02),log(p33/p03),log(p34/p04),log(p35/p05),log(p36/p06))
out3<-lm(y3~x)
alpha3_ini<-summary(out3)$coeff[[1]]
beta3_ini<-summary(out3)$coeff[[2]]

ini_par<-c(9,15,alpha1_ini,beta1_ini,alpha2_ini,beta2_ini,alpha3_ini,beta3_ini)
#========== Optimization functions for three different methods =============

OP_NEW<-optim(ini_par, fNEW, method="L-BFGS-B", lower=c(1,1,-20,-20,-20,-20,-20,-20), upper=c(25,25,20,20,20,20,20,20), hessian=T,control = list(maxit=1000000,trace = T,REPORT = 1))
OP_NEW_EXNR<-optim(ini_par[1:4],fNEW_EXNR, method="L-BFGS-B", lower=c(1,1,-20,-20), upper=c(25,25,20,20), hessian=T,control = list(maxit=1000000,trace = T,REPORT = 1))
OP_YN<-optim(c(9,15), fYNO, method="BFGS",hessian=T,control = list(maxit=1000000,trace = T,REPORT = 1))

parnew<-OP_NEW$par
parnew_exnr<-OP_NEW_EXNR$par
paryn<-OP_YN$par

parnew
parnew_exnr
paryn

hsnMNEW<-OP_NEW$hessian
hsnMNEWEXNR<-OP_NEW_EXNR$hessian
hsnMYN<-OP_YN$hessian

dmNEW<-solve(hsnMNEW)
dmNEWEXNR<-solve(hsnMNEWEXNR)
dmYN<-solve(hsnMYN)

#===== Survival plot  ======
theta1<-paryn[1]
eta1<-paryn[2]
theta2<-parnew_exnr[1]
eta2<-parnew_exnr[2]
theta3<-parnew[1]
eta3<-parnew[2]

median1<-eta1*(log(2)^(1/theta1))
median2<-eta2*(log(2)^(1/theta2))
median3<-eta3*(log(2)^(1/theta3))

age<-seq(6,18,.1)
su1<-pweibull(age, theta1,eta1,lower.tail=FALSE)
su2<-pweibull(age, theta2,eta2,lower.tail=FALSE)
su3<-pweibull(age, theta3,eta3,lower.tail=FALSE)

probability_nomenarche<-su1


postscript (file="SURV1.eps", width = 8.0, height = 6.0)
plot(age,probability_nomenarche,cex.lab=1, cex.axis=1, cex.main=1,col="white",ylab="probability of no menarche")
 lines(age,su1,lty=4)
lines(age,su2,col=4)
lines(age,su3,lty=5,col=2)      
legend(6,.5, c("Current Status MLE","Binary Recall MLE", "Partial Recall MLE"), col = c(1,4,2),
       lty = c(4,1,5), pch = c(NA,NA,NA),merge = TRUE ,bty="n",cex=.9)
dev.off()

