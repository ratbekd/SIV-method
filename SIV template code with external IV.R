rm(list=ls()) #Removes all items in Environment!
##Libraries
library(stargazer)
library(dplyr)
library(lmtest)
library(haven)
library(zoo)
library(sandwich)
library(rMR)
library(ivreg)
##Functions used
# Two-sample Anderson-Darling statistic
ad2_stat <- function(x, y) {
  # Sample sizes
  n <- length(x)
  m <- length(y)
  
  # Pooled sample and pooled ecdf
  z <- c(x, y)
  z <- z[-which.max(z)] # Exclude the largest point
  H <- rank(z) / (n + m)
  
  # Statistic computation via ecdf()
  (n * m / (n + m)^2) * sum((ecdf(x)(z) - ecdf(y)(z))^2 / ((1 - H) * H))
  
}
### The function to test absolute non-decreasing trend and sign change for locus E(ee'siv)
# check_initial_abs_increase_and_sign_change <- function(x) {
#   # Step 1: Check if the absolute values start by increasing (non-decreasing trend)
#   abs_x <- abs(x)
#   initial_increasing <- all(diff(abs_x[1:min(20, length(abs_x))]) >= 0)  # Check first 3 points or as many as available
#   
#   # Step 2: Check if there is a sign change
#   signs <- sign(x)  # Get sign (-1, 0, or 1)
#   sign_changes <- any(diff(signs) != 0, na.rm = TRUE)  # Detect if sign changes
#   
#   # Step 3: Return 1 if both conditions are met, otherwise 0
#   return(as.integer(sign_changes))
# }
check_sign_change <- function(x) {
  # Step 2: Check if there is a sign change
  signs <- sign(x)  # Get sign (-1, 0, or 1)
  sign_changes <- any(diff(signs) != 0, na.rm = TRUE)  # Detect if sign changes
  
  # Step 3: Return 1 if both conditions are met, otherwise 0
  return(as.integer(sign_changes))
}


check_initial_abs_increase <- function(x) {
  # Step 1: Check if the absolute values start by increasing (non-decreasing trend)
  abs_x <- abs(x)
  initial_increasing <- all(diff(abs_x[1:min(20, length(abs_x))]) >= 0)  # Check first 3 points or as many as available
  
  # Step 3: Return 1 if both conditions are met, otherwise 0
  return(as.integer(initial_increasing))
}
find_first_sign_change <- function(x) {
  sign_changes <- which(diff(sign(x)) != 0)  # Find indices where sign changes
  if (length(sign_changes) > 0) {
    return(sign_changes[1] + 1)  # Return the first occurrence (adjust for diff)
  } else {
    return(NA)  # Return NA if no sign change
  }
}

####### DATA 
library(wooldridge)
#
mydata<-wooldridge::mroz
mydata<-mydata[complete.cases(mydata), ]
attach(mydata)

############Setting the regression variables. You need to input your variables here
H0 <- data.frame(hours, lwage,educ, age,kidslt6, kidsge6, nwifeinc)  
# Vector of variables used in the regression. Firts var is outcome, the second var is endogenous regressor
IV <- data.frame(exper,expersq) ####IV variables need to be iserted here. If you do not have external IVs use the other version of the code.
Y <- as.character(colnames(H0))[1] ###OUTCOME variable
X <-as.character(colnames(H0))[2] ###ENDOEGENOUS variable
H<- as.character(colnames(H0))[-(1:2)] ##### EXOGENOUS variables
#print(H)
formula_str <- paste(paste0(Y," ~ ", X,"+"), paste(H, collapse = " + "))## Construct formula as a string
formula <- as.formula(formula_str)# Convert to formula object
#print(formula)
iv_str <- as.character(colnames(IV))# Construct the external IVs as a string
ivs <- paste(paste(iv_str, collapse = " + "))# convert into part of the instruments object
instruments_str <- paste0(" ~ ",ivs, " + ", paste(H, collapse = " + "))# Construct instruments as a string
instruments <- as.formula(instruments_str)# Convert to instrument object
# determine the number of rows
N<-nrow(mydata)
#####Traditional methods
### OLS estimation
ols1 <- lm(formula, data = mydata)
summ.ols1 <- summary(ols1, vcov. = function(x) vcovHC(x, type="HC1"),
                     diagnostics=T)
##IV regression
iv1<-ivreg(formula, instruments, data=mydata)
summ.iv1 <- summary(iv1,      diagnostics=T)

#####SIV Method

###Basic vectors for  SIV calculation###
y1<-hours ### the outcome variabe
## Factoring out the effects of other exogenous variables
formula_str <- paste(paste0(Y," ~ "), paste(H, collapse = " + "))## Construct formula as a string
formula <- as.formula(formula_str)# Convert to formula object
fity<-lm(formula, data=mydata)
y<-resid(fity)
x1<-lwage### the endogenous variable
## Factoring out the effects of other exogenous variables
formula_str <- paste(paste0(X," ~ "), paste(H, collapse = " + "))## Construct formula as a string
formula <- as.formula(formula_str)# Convert to formula object
fitx<-lm(formula, data=mydata)
x<-resid(fitx)
#saving the transformed x and y
mydata$x<-(x-mean(x))
mydata$y<-(y-mean(y))
y0<-y
x0<-x
V=0
### Generating a vector orthogonal to x
fity<-lm(y0~(x0), data=mydata)
V<-resid(fity)
V<-(V-mean(V))/sd(V)
V<-V*sd(x0)
mydata$R<-V
############ Determining the sign of cor(x,u)
k=0
j=1
signc=matrix(ncol = 3, nrow = 2)
signc[1,1] <-1
signc[2,1] <--1
for (j in 1:2) {
  if(j<2){k=1}else{k=-1}#the assumed sign for cor(x,u)
  
  # IV regression
  data <- mydata
  dd<-4#end value for delta
  d<-0.01# starting value for delta
  delt<-0.01# step to change delta
  i<-1 ### starting value for a counter
  ## placeholders for variables
  m1=0
  while (d<dd){# we compute m1 until siv is close to become perpendicular to x
    data$siv<-(data$x-k*d*data$R) ### SIV
    rls<-(lm(x~siv, data=data))
    s.rls<-summary(rls,vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
    data$ev21<-resid(s.rls)
    m1[i]<-cov(data$ev21^2,data$siv)
    d<-d+delt
    i<-i+1
  }
  plot(m1)
  m=m1
  
  signc[j,3] <- check_initial_abs_increase(m)#result
  if(signc[j,3]!=1){
    index <- find_first_sign_change(m)
    if(is.na(index)){index=0}
    m=m1[index:i]
  }else{m=m1}
  signc[j,2] <-check_sign_change(m)
  signc[j,3] <- check_initial_abs_increase(m)#result
}

ch=0
for(j in 1:2){cat("the assumed sign for cor(x,u):", signc[j,1], if(signc[j,2]*signc[j,3]==0)
{"is FALSE"
}else{"is TRUE"},  "\n") 
  ch[j] <- signc[j,2]*signc[j,3]
} # TRUE}}
k <- signc[which.max(ch),1]
#s<-1 # the assumed sign for cor(x,u)
if(k!=0){
############# Initial settings for SIV
vvar=0
d0i=0
d0ri=0
d0rni=0
b2=0
b2r=0
b2rp=0
b2t=0
N<-nrow(mydata)
reps=2
S=round(N*.999)
l=1
fitc <- matrix(ncol = 1, nrow = reps)
sumb2=matrix(ncol = 1, nrow = reps)
fitcr <- matrix(ncol = 1, nrow = reps)
sumb2r=matrix(ncol = 1, nrow = reps)
fitcrn <- matrix(ncol = 1, nrow = reps)
sumb2rn=matrix(ncol = 1, nrow = reps)
fitct <- matrix(ncol = 1, nrow = reps)
sumb2t=matrix(ncol = 1, nrow = reps)
lowbp=0
upbp=0

##### Bootstrap sampling loop. You may use data <- mydata instead of data <- mydata[sample(1:N, S),  TRUE]
#if you just want see how it works for the original sample data.
while (l<reps){
  # IV regression
  set.seed(3*(l))   # a different seed for each sub sample  
  data <- mydata[sample(1:N, S),  TRUE]
  ####Computation of m1 
  
  dd<-3#end value for delta
  d<-0.01# starting value for delta
  delt<-0.01# step to change delta
  i<-1 ### starting value for a counter
  ## placeholders for variables
  m1=0
  st=0
  ev22=0
  dv=0
  dv2=0
  x4=0
  l1=0
  l2=0
  while (d<dd){# we compute m1 until siv is close to become perpendicular to x
    data$siv<-(data$x-k*d*data$R) ### SIV
    ####OLS estimates
    rls<-(lm(x~siv, data=data))
    s.rls<-summary(rls,vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
    data$ev21<-resid(s.rls)
    ### FGLS estimates
    ehatsq <- resid(rls)^2
    sighatsq.ols  <- lm(log(ehatsq)~siv, data=data)
    data$vari <- sqrt(exp(fitted(sighatsq.ols)))
    vvar[i] <- var(data$vari)
    fgls <-lm(x~siv, weights=1/vari, data=data)
    data$ev22<-resid(fgls)
    ##   ## Homoscedastic estimate for a simple SIV case
    m1[i]<-(cor(data$ev21^2,data$siv))
    ### Parametric computations for heteroscedastic case
    n=length(data$ev21)
    l1 <- summary(lm((ev21^2)~siv,data=data))
    l2 <-summary( lm((ev22^2)~siv,data=data))
    ssr1 <- sumsq(predict(lm((ev21^2)~siv,data=data))-mean(data$ev21^2))
    sse1=sumsq(data$ev21)
    x1= (ssr1/2)/(sse1/n^2)^2
    ssr2 <- sumsq(predict(lm((ev22^2)~siv,data=data))-mean(data$ev22^2))
    sse2=sumsq(data$ev22)
    x2= (ssr2/2)/(sse2/n^2)^2
    dv[i] <-pchisq(x2, df =1,lower.tail=FALSE)-pchisq(x1, df =1,lower.tail=FALSE)#
    x3 <- x1/x2#sumsq(predict(lm((ev22^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data)))
    dv2[i] <- pf(x3, df = 1, df2 = 1, lower.tail = TRUE)
    # Non-parametric CDF computations for heterscedastic case x3 <- x1/x2#sumsq(predict(lm((ev22^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data))/predict(lm((ev21^2)~siv,data=data)))
      dv2[i] <- pf(x3, df = 1, df2 = 1, lower.tail = TRUE)
    samp1 <- (predict(lm((ev21^2)~siv,data=data)))^2
    samp2 <- (predict(lm((ev22^2)~siv,data=data)))^2
    xx0 <- samp1
    yy0 <- samp2
    ad0 <- ad2_stat(x = xx0, y = yy0)
    x4[i] <- 1-ad0  
    st[i]<-d
    d<-d+delt
    i<-i+1
  }
  ### updating the formula for regressions
  formula_str <- paste(paste0(Y," ~ ", X,"+"), paste(H, collapse = " + "))## Construct formula as a string
  formula <- as.formula(formula_str)# Convert to formula object
  ### update with your own instruments: instruments<-~ siv+all exogenous variables
  iv_str <- paste(paste0(" ~ ", "siv","+"), paste(H, collapse = " + "))## Construct formula as a string
  instruments <- as.formula(iv_str)# Convert to formula object
  #instruments<-~siv+educ+ age+kidslt6+ kidsge6+ nwifeinc
  # formula <-hours~lwage+educ+ age+kidslt6+ kidsge6+ nwifeinc
  #### DT condition of homoscedatic case
  d0 <- (which.min(abs(m1)))*delt
  d0i[l] <- d0
  data$siv<-(data$x-k*d0*data$R)
  iv2<-ivreg(formula, instruments, data=data)
  summ.iv2 <- summary(iv2, diagnostics=T)#, vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
  ## saving the estimation parameters for each sample
  fitc[l] <- iv2$coefficients[2]
  sumb2[l] <-  summ.iv2$coefficients[2,2]
  
  ### DT point for heteroscedastic case- parametric approach
  d0r <-  which.min(dv2)*delt
  d0ri[l] <- d0r
  data$siv<-(data$x-k*d0r*data$R)
  iv3<-ivreg(formula, instruments, data=data)
  summ.iv3 <- summary(iv3, diagnostics=T)#,  vcov. = function(x) vcovHC(x, type="HC1"), diagnostics=T)
  ## saving the estimation paramters for each sample
  fitcr[l] <- iv3$coefficients[2]
  sumb2r[l] <-  summ.iv3$coefficients[2,2]
  
  #### DT point for heteroscedastic case- non-parametric approach
  d0rn <- which.min(x4)*delt
  d0rni[l] <- d0rn
  data$siv<-(data$x-k*d0rn*data$R)
  iv4<-ivreg(formula, instruments, data=data)
  summ.iv4 <- summary(iv4, diagnostics=T)#, 
  ## saving the estimation paramters for each sample
  fitcrn[l] <- iv4$coefficients[2]
  sumb2rn[l] <-  summ.iv4$coefficients[2,2]
  l <- l+1
}

## Distribution of sample paramters
# the simple homogenous case
fitc <- fitc[complete.cases(fitc)]
sumb2<- sumb2[complete.cases(sumb2)]
alpha <- 0.05 # chosen significance level
b2[1] <- mean(fitc)# the endogenous parameter beta
df <- df.residual(iv2)# degrees of freedom for t-stat
seb2 <-mean(sumb2)# the average standard error of the parameter beta
tc <- qt(1-alpha/2, df) ## t-statistic
lowbp[1] <- b2[1]-tc*seb2  # lower bound for beta
upbp[1] <- b2[1]+tc*seb2   # upper bound for beta

#### The parametric heterogenous case
fitcr <- fitcr[complete.cases(fitcr)]
sumb2r<- sumb2r[complete.cases(sumb2r)]
alpha <- 0.05 # chosen significance level
b2[2] <- mean(fitcr)
df <- df.residual(iv3)
seb2r <-mean(sumb2r)
tc <- qt(1-alpha/2, df)
lowbp[2] <- b2[2]-tc*seb2r  # lower bound
upbp[2] <- b2[2]+tc*seb2r   # upper bound

#################The non-parametric heterogenous case
fitcrn <- fitcrn[complete.cases(fitcrn)]
sumb2rn<- sumb2rn[complete.cases(sumb2rn)]
alpha <- 0.05 # chosen significance level
b2[3] <- mean(fitcrn)
df <- df.residual(iv4)
seb2rn <-mean(sumb2rn)
tc <- qt(1-alpha/2, df)
lowbp[3] <- b2[3]-tc*seb2rn  # lower bound
upbp[3] <- b2[3]+tc*seb2rn   # upper bound

### Table for CI of beta
mv<-data.frame(lowbp,b2,upbp)
colnames(mv)<-c("low beta","mean b2", "high beta") 
rownames(mv)<- c("SIV","SIVRr","SIVRn")#, "nearc4")
### 
(mv)

### final satge estimations
###### Simple homogenous assumption case
d0i <-  d0i[complete.cases(d0i)]
d0m <- mean(d0i)
mydata$siv<-(mydata$x-k*d0m*mydata$R)
iv2<-ivreg(formula, instruments, data=mydata)
summ.iv2 <- summary(iv2, diagnostics=T)#
############Paramteric heterogenous case
d0ri <-  d0ri[complete.cases(d0ri)]
d0rm <- mean(d0ri)
mydata$siv<-(mydata$x-k*d0rm*mydata$R)
iv3<-ivreg(formula, instruments, data=mydata)
summ.iv3 <- summary(iv2, diagnostics=T)#

#### Non-parameteric heterogenous case
d0rni <-  d0rni[complete.cases(d0rni)]
d0rnm <- mean(d0rni)
mydata$siv<-(mydata$x-k*d0rnm*mydata$R)
iv4<-ivreg(formula, instruments, data=mydata)
summ.iv4 <- summary(iv2, diagnostics=T)#
}else{
  print("NO endogeneity problem. All SIV estimates are the same as the OLS")
  d0m=0.001
mydata$siv<-(mydata$x-k*d0m*mydata$R)
iv2<-ivreg(formula, instruments, data=mydata)
summ.iv2 <- summary(iv2, diagnostics=T)#
  d0rm=0.001
  mydata$siv<-(mydata$x-k*d0rm*mydata$R)
  iv3<-ivreg(formula, instruments, data=mydata)
  summ.iv3 <- summary(iv2, diagnostics=T)#
  
  d0rnm=0.001
  mydata$siv<-(mydata$x-k*d0rm*mydata$R)
  iv4<-ivreg(formula, instruments, data=mydata)
  summ.iv4 <- summary(iv4, diagnostics=T)#
  }
# The estimation output
stargazer(ols1, iv1, iv2, iv3, iv4,  # Include iv4 in the list of models
          type = "text",
          omit = "reg",
          dep.var.caption = "Work hours",
          dep.var.labels.include = FALSE,
          model.numbers = FALSE,
          model.names = FALSE,
          column.labels = c("OLS", "IV", "SIV", "RSIV-p", "RSIV-n"),  # Add a label for iv4
          no.space = TRUE,
          add.lines = list(
            c("Weak instruments", "",
              round(summ.iv1$diagnostics[1, "p-value"], 2),
              round(summ.iv2$diagnostics[1, "p-value"], 2),
              round(summ.iv3$diagnostics[1, "p-value"], 2),
              round(summ.iv4$diagnostics[1, "p-value"], 2)),  # Add p-value for iv4
            c("Wu-Hausman", "",
              round(summ.iv1$diagnostics[2, "p-value"], 2),
              round(summ.iv2$diagnostics[2, "p-value"], 2),
              round(summ.iv3$diagnostics[2, "p-value"], 2),
              round(summ.iv4$diagnostics[2, "p-value"], 2)),  # Add p-value for iv4
            c("Sargan", "",
              round(summ.iv1$diagnostics[3, "p-value"], 2),
              round(summ.iv2$diagnostics[3, "p-value"], 2),
              round(summ.iv3$diagnostics[3, "p-value"], 3),
              round(summ.iv4$diagnostics[3, "p-value"], 3))  # Add p-value for iv4
          ),
          multicolumn = FALSE)
###End of the main code
