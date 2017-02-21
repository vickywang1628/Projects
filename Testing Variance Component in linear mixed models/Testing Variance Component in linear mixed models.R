# ST 758 Homework 7 chong Wang
# clean the memory
rm(list=ls())
library(RLRsim)
library(Matrix)
library(MASS)  
library(nlme)  
library(RLRsim)

fun <- function(x){
  z <- rep(1,x)
  return(z)
}

lrt.test <- function(p,sigma){
  total <- 50
  n <- sum(p)
  z1 <- bdiag(sapply(p,fun))
  v1 <- crossprod(t(z1),t(z1))
  v2 <- diag(1,n)
  mu <- as.vector(rep(1,n))
  group <- rowSums(t(c(1:length(p))*t(z1)))
  v <- sigma*v1+v2 
  #set.seed(16) 
  ######################################
  rej1 <- 0
  rej2 <- 0
  rej3 <- 0
  for(i in 1:total){ 
    set.seed(758*i)
    y <- mvrnorm(1, mu, v)
    data <- data.frame(y,group)
    names(data)<-c("y", "group")
    data$group<-as.factor(data$group)
    fit0 <- lm(y  ~ 1, data=data, method='qr')
    fita <- lme(y ~1,data=data,random=~1|group,method="ML")
    fitaR <- lme(y ~1,data=data,random=~1|group,method="REML")
    #Test 1####################################################
    pvalue <- 1/2 - 1/2*pchisq(-2*logLik(fit0)[1] + 
                                 2*logLik(fita)[1],df=1,lower.tail=TRUE)
    if (pvalue < 0.05) {
      rej1 <- rej1 + 1
    }
    power1 <- rej1/total
    sd1 <- sqrt(power1*(1-power1)/total)
    #Test 2#################################################### 
    LRT <- suppressMessages(exactLRT(m=fita,m0=fit0))
    if (LRT$p < 0.05) {
      rej2 <- rej2 + 1
    }
    power2 <- rej2/total
    sd2 <-sqrt( power2*(1-power2)/total)
    #Test 3#################################################### 
    RLRT <- exactRLRT(m=fitaR,mA=fitaR,m0=fit0)
    if (RLRT$p < 0.05) {
      rej3 <- rej3 + 1
    }
    power3 <- rej3/total
    sd3 <- sqrt(power3*(1-power3)/total)
 }#########################################
  power <- c(power1,power2,power3)
  sd <- c(sd1,sd2,sd3)
  return(list(power=power,sd=sd))
}

sigma <- c(0,0.1,0.2,0.5,1,2,5)
#pattern 1 power
pp1 <- NULL
sd1 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(3,5,7),sigma[i])$power
  pp1 <- cbind(pp1,temp1)
  temp2 <- lrt.test(c(3,5,7),sigma[i])$sd
  sd1 <- cbind(sd1,temp2)
}

plot(sigma,pp1[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 1: Power Curve")
lines(sigma,pp1[1,],col="red",lwd=1,lty=2)
lines(sigma,pp1[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(0.5,1,2),lty=c(1,2,2),cex=0.65)
#use the multiple paired t test to test the significance of the power difference
power1 <- c(pp1[1,],pp1[2,],pp1[3,])
group1<-c(rep("test1",7),rep("test2",7),rep("test3",7))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)


#pattern 2 power
pp2 <- NULL
sd2 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(1,5,9),sigma[i])$power
  pp2 <- cbind(pp2,temp1)
  temp2 <- lrt.test(c(1,5,9),sigma[i])$sd
  sd2 <- cbind(sd2,temp2)
}

#power curve
plot(sigma,pp2[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 2: Power Curve")
lines(sigma,pp2[1,],col="red",lwd=1,lty=2)
lines(sigma,pp2[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power2 <- c(pp2[1,2:7],pp2[2,2:7],pp2[3,2:7])
group2<-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 4 power
pp4 <- NULL
sd4 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(3,3,5,5,7,7),sigma[i])$power
  pp4 <- cbind(pp4,temp1)
  temp2 <- lrt.test(c(3,3,5,5,7,7),sigma[i])$sd
  sd4 <- cbind(sd4,temp2)
}

#power curve
plot(sigma,pp4[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 4: Power Curve")
lines(sigma,pp4[1,],col="red",lwd=1,lty=2)
lines(sigma,pp4[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power4 <- c(pp4[1,2:7],pp4[2,2:7],pp4[3,2:7])
group4<-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 5 power
pp5 <- NULL
sd5 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(1,1,5,5,9,9),sigma[i])$power
  pp5 <- cbind(pp5,temp1)
  temp2 <- lrt.test(c(1,1,5,5,9,9),sigma[i])$sd
  sd5 <- cbind(sd5,temp2)
}

#power curve
plot(sigma,pp5[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 5: Power Curve")
lines(sigma,pp5[1,],col="red",lwd=1,lty=2)
lines(sigma,pp5[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power5 <- c(pp5[1,2:7],pp5[2,2:7],pp5[3,2:7])
group5<-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 7 power
pp7 <- NULL
sd7 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(1,1,1,1,13,13),sigma[i])$power
  pp7 <- cbind(pp7,temp1)
  temp2 <- lrt.test(c(1,1,1,1,13,13),sigma[i])$sd
  sd7 <- cbind(sd7,temp2)
}

#power curve
plot(sigma,pp7[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 7: Power Curve")
lines(sigma,pp7[1,],col="red",lwd=1,lty=2)
lines(sigma,pp7[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power7 <- c(pp7[1,2:7],pp7[2,2:7],pp7[3,2:7])
group7 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 8 power
pp8 <- NULL
sd8 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(3,3,3,5,5,5,7,7,7),sigma[i])$power
  pp8 <- cbind(pp8,temp1)
  temp2 <- lrt.test(c(3,3,3,5,5,5,7,7,7),sigma[i])$sd
  sd8 <- cbind(sd8,temp2)
}

#power curve
plot(sigma,pp8[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 8: Power Curve")
lines(sigma,pp8[1,],col="red",lwd=1,lty=2)
lines(sigma,pp8[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power8 <- c(pp8[1,2:7],pp8[2,2:7],pp8[3,2:7])
group8 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 9 power
pp9 <- NULL
sd9 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(1,1,1,5,5,5,9,9,9),sigma[i])$power
  pp9 <- cbind(pp9,temp1)
  temp2 <- lrt.test(c(1,1,1,5,5,5,9,9,9),sigma[i])$sd
  sd9 <- cbind(sd9,temp2)
}

#power curve
plot(sigma,pp9[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 9: Power Curve")
lines(sigma,pp9[1,],col="red",lwd=1,lty=2)
lines(sigma,pp9[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power9 <- c(pp9[1,2:7],pp9[2,2:7],pp9[3,2:7])
group9 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 11 power
pp11 <- NULL
sd11 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(1,1,1,1,1,1,1,19,19),sigma[i])$power
  pp11 <- cbind(pp11,temp1)
  temp2 <- lrt.test(c(1,1,1,1,1,1,1,19,19),sigma[i])$sd
  sd11 <- cbind(sd11,temp2)
}

#power curve
plot(sigma,pp11[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 11: Power Curve")
lines(sigma,pp11[1,],col="red",lwd=1,lty=2)
lines(sigma,pp11[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power11 <- c(pp11[1,2:7],pp11[2,2:7],pp11[3,2:7])
group11 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 12 power
pp12 <- NULL
sd12 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(2,10,18),sigma[i])$power
  pp12 <- cbind(pp12,temp1)
  temp2 <- lrt.test(c(2,10,18),sigma[i])$sd
  sd12 <- cbind(sd12,temp2)
}

#power curve
plot(sigma,pp12[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 12: Power Curve")
lines(sigma,pp12[1,],col="red",lwd=1,lty=2)
lines(sigma,pp12[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power12 <- c(pp12[1,2:7],pp12[2,2:7],pp12[3,2:7])
group12 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)

#pattern 13 power
pp13 <- NULL
sd13 <- NULL
for(i in 1:7){
  temp1 <- lrt.test(c(3,15,27),sigma[i])$power
  pp13 <- cbind(pp13,temp1)
  temp2 <- lrt.test(c(3,15,27),sigma[i])$sd
  sd13 <- cbind(sd13,temp2)
}

#power curve
plot(sigma,pp13[3,],type="l",xlab=expression(lambda),ylab="power",lwd=0.5,lty=1,main="Pattern 13: Power Curve")
lines(sigma,pp13[1,],col="red",lwd=1,lty=2)
lines(sigma,pp13[2,],lwd=2,lty=2)
legend("bottomright",c("LRT0","LRT1","RLRT"),lwd=c(1,2,0.5),lty=c(2,2,1),cex=0.65)

#use the multiple paired t test to test the significance of the power difference
power13 <- c(pp13[1,2:7],pp13[2,2:7],pp13[3,2:7])
group13 <-c(rep("test1",6),rep("test2",6),rep("test3",6))
pairwise.t.test(power1,group1,p.adjust="bonf",paired=TRUE)





