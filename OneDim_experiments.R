rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
} #This function generates a data set of size n from the one dimensional TSIV model if the matrix 'coeff' is a 4*4 matrix with correct nonzero entries. 
#The input 'error' is an n*4 matrix containing the error terms

TimeOperator <- function(vector,timepoint) {
  if (timepoint==0) {
    vector <- vector[8:length(vector)]
  }
  else if (timepoint==1) {
    vector <- vector[7:(length(vector)-1)]
  }
  else if (timepoint==2) {
    vector <- vector[6:(length(vector)-2)]
  }
  else if (timepoint==3) {
    vector <- vector[5:(length(vector)-3)]
  }
  else if (timepoint==4) {
    vector <- vector[4:(length(vector)-4)]
  }
  else if (timepoint==5) {
    vector <- vector[3:(length(vector)-5)]
  }
  else if (timepoint==6) {
    vector <- vector[2:(length(vector)-6)]
  }
  else if (timepoint==7) {
    vector <- vector[1:(length(vector)-7)]
  }
  return(vector)
} #This function shortens a data vector according to the time lag it represents

#The following four values/matrix can be changed to create the desired combination
epsilon <- 0.05 #Epsilon is the bound we consider when investigating consistency 
f<- 1           #f is a parameter that influences the number of time steps, n, included in the simulated data
s <- 100        # s is the number of simulations made pr. number of time step
A <- matrix(c(0.5,-0.9,0,0,0,0.5,0,-0.5,0,-0.8,-0.8,-0.8,0,0,0,-0.5),ncol=4) #The data generating coefficient matrix from model M1


P_cov <- numeric(14) #These are vectors that will contain the estimated probability for the different estimators
P_2sls <- numeric(14)
P_direkte <- numeric(14)
P_covI1I5 <- numeric(14)
P_covI1I3 <- numeric(14)
P_covI1I2 <- numeric(14)
P_covI1I7 <- numeric(14)
P_covI1I7deltaY <- numeric(14) 
S <-c(1600,3200,4800,6400,8000,9600,11200,12800,14400,16000,17600,19200,20800,22400)*f #S defines the different number of time steps that will be used
for (j in S ) {
  B_2sls <- numeric(s) #These are vectors that for each time step in S will contain 's' etatimates from the given estimator
  B_cov <- numeric(s)
  B_direkte <- numeric(s)
  B_covI1I5X <- numeric(s)
  B_covI1I3X <- numeric(s)
  B_covI1I2X <- numeric(s)
  B_covI1I7X <- numeric(s)
  B_deltaY <- numeric(s)
  for (i in 1:s) {
    E <- matrix(rnorm(j,0,0.1),ncol=4) #All variables are assumed to have standard deviation 0.1
    sim1data <- rSim(A,E) #Here we simulate a data set.
    Y0 <- TimeOperator(sim1data[,4],0) #Here the data vectors are shortened according to their time lag. 
    Y1 <- TimeOperator(sim1data[,4],1)
    X1 <- TimeOperator(sim1data[,2],1) 
    X2 <- TimeOperator(sim1data[,2],2) 
    X3 <- TimeOperator(sim1data[,2],3)
    I1 <- TimeOperator(sim1data[,1],1)
    I2 <- TimeOperator(sim1data[,1],2) 
    I3 <- TimeOperator(sim1data[,1],3)
    I4 <- TimeOperator(sim1data[,1],4)
    I5 <- TimeOperator(sim1data[,1],5)
    I6 <- TimeOperator(sim1data[,1],6)
    I7 <- TimeOperator(sim1data[,1],7)
    #TS-2SLS estimate - stage 1 
    model <- lm(X1 ~ 0+I2+I3)
    r <- coefficients(model)[1]*I2
    #TS-2SLS estimate stage 2
    model2 <- lm(Y0 ~ 0+r+I3)
    B_2sls[i] = (coefficients(model2)[1])
    #TS-IR estimate including I2-I3
    B_cov[i] <- (cov(Y0,I3)-(cov(Y0,I2)*cov(Y1,I3))/cov(Y1,I2))/(cov(X1,I3)-(cov(X1,I2)*cov(Y1,I3))/cov(Y1,I2))
    #TS-IR estimate including I1 to I5
    Cov2_matrixI1I5 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(X1,I4),cov(X1,I5),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3),cov(Y1,I4),cov(Y1,I5)),ncol=2)
    Cov2_vectorI1I5 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3),cov(Y0,I4),cov(Y0,I5))
    cov2_est_vectorI1I5 <- qr.solve(Cov2_matrixI1I5,Cov2_vectorI1I5)
    B_covI1I5X[i] <-cov2_est_vectorI1I5[1]  
    #TS-IR estimate including I1 to I3
    Cov2_matrixI1I3 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3)),ncol=2)
    Cov2_vectorI1I3 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3))
    cov2_est_vectorI1I3 <- qr.solve(Cov2_matrixI1I3,Cov2_vectorI1I3)
    B_covI1I3X[i] <-cov2_est_vectorI1I3[1] 
    #TS-IR estimate including I1 to I2
    Cov2_matrixI1I2 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(Y1,I1),cov(Y1,I2)),ncol=2)
    Cov2_vectorI1I2 <- c(cov(Y0,I1),cov(Y0,I2))
    cov2_est_vectorI1I2 <- qr.solve(Cov2_matrixI1I2,Cov2_vectorI1I2)
    B_covI1I2X[i] <-cov2_est_vectorI1I2[1] 
    #Naive estimate
    model3 <- lm(Y0~0+X1+Y1)
    B_direkte[i] <- coefficients(model3)[1]
    #TS-IR estimate including I1 to I7 of beta and deltaY
    Cov2_matrixI1I7 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(X1,I4),cov(X1,I5),cov(X1,I6),cov(X1,I7),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3),cov(Y1,I4),cov(Y1,I5),cov(Y1,I6),cov(Y1,I7)),ncol=2)
    Cov2_vectorI1I7 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3),cov(Y0,I4),cov(Y0,I5),cov(Y0,I6),cov(Y0,I7))
    cov2_est_vectorI1I7 <- qr.solve(Cov2_matrixI1I7,Cov2_vectorI1I7)
    B_covI1I7X[i] <-cov2_est_vectorI1I7[1] 
    B_deltaY[i] <- cov2_est_vectorI1I7[2]
  }
  k <- j/(1600*f) #Creating an index from 1 to 14
  P_2sls[k] <-mean(abs(B_2sls-A[4,2])>epsilon) #Here we calculate the estimated probability that the estimate is further from beta than epsilon
  P_cov[k] <- mean(abs(B_cov-A[4,2])>epsilon)
  P_direkte[k] <- mean(abs(B_direkte-A[4,2])>epsilon)
  P_covI1I5[k] <- mean(abs(B_covI1I5X-A[4,2])>epsilon)
  P_covI1I3[k] <- mean(abs(B_covI1I3X-A[4,2])>epsilon)
  P_covI1I2[k] <- mean(abs(B_covI1I2X-A[4,2])>epsilon)
  P_covI1I7[k] <- mean(abs(B_covI1I7X-A[4,2])>epsilon)
  P_covI1I7deltaY[k] <- mean(abs(B_deltaY-A[4,4])>epsilon)
}
N <- c(1600,3200,4800,6400,8000,9600,11200,12800,14400,16000,17600,19200,20800,22400)/4*f #N is a list of the number of time steps used. We devide by 4 since the error matrix has dimension n*4
H <- cbind(N,P_direkte,P_2sls,P_covI1I2,P_cov, P_covI1I3,P_covI1I5,P_covI1I7,P_covI1I7deltaY) #Here we collect all the estimated probabilities

#Below we create the different figures
library(tidyverse)
library(reshape2)
library(ggplot2)
library(tikzDevice)
data2 <- as.data.frame(H)
names(data2) <- c("N","Naive", "2SLS", "I1-I2","I2-I3","I1-I3","I1-I5","I1-I7","I1-I7Y")
data2sls <- data.frame(data2$N,data2$`2SLS`)
dataIR <- data.frame(data2$N,data2$`I1-I2`,data2$`I2-I3`,data2$`I1-I3`,data2$`I1-I5`,data2$`I1-I7`)
data2slsnaive <- data.frame(data2$N,data2$`2SLS`,data2$Naive)
dataIRnaive <- data.frame(data2$N,data2$`I1-I3`,data2$`I1-I7`,data2$Naive)
data2slsIR <- data.frame(data2$N,data2$`2SLS`, data2$`I1-I3`,data2$`I1-I7`)
datanaiveIR2sls <- data.frame(data2$N,data2$`I1-I3`,data2$`I1-I7`,data2$`2SLS`,data2$Naive)
datanaive <- data.frame(data2$N,data2$Naive)
dataY <- data.frame(data2$N,data2$`I1-I7`,data2$`I1-I7Y`)


names(data2sls) <-c("N", "2SLS") 
names(dataIR) <- c("N", "I1-I2","I2-I3","I1-I3","I1-I5","I1-I7")
names(data2slsnaive) <- c("N", "2SLS","Naive")
names(dataIRnaive) <- c("N","I1-I3","I1-I7","Naive")
names(data2slsIR) <- c("N", "2SLS","I1-I3","I1-I7")
names(datanaiveIR2sls) <- c("N","I1-I3","I1-I7","2SLS","Naive")
names(datanaive) <- c("N","Naive")
names(dataY) <- c("N","I1-I7","I1-I7Y")

redata2 <- reshape2::melt(data2, id.var = "N")
redata2sls <- reshape2::melt(data2sls, id.var = "N")
redataIR <- reshape2::melt(dataIR, id.var = "N")
redata2slsnaive <- reshape2::melt(data2slsnaive, id.var = "N")
redataIRnaive <- reshape2::melt(dataIRnaive, id.var = "N")
redata2slsIR <- reshape2::melt(data2slsIR, id.var = "N")
redatanaiveIR2sls <- reshape2::melt(datanaiveIR2sls, id.var = "N")
redatanaive <- reshape2::melt(datanaive, id.var = "N")
redataY <- reshape2::melt(dataY, id.var = "N")

P_direkte_low <- data2$Naive-1.96*sqrt(data2$Naive*(1-data2$Naive)/s)
P_direkte_high <- data2$Naive+1.96*sqrt(data2$Naive*(1-data2$Naive)/s)

P_2sls_low <- data2$`2SLS`-1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)
P_2sls_high <- data2$`2SLS`+1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)

P_I1I2_low <- data2$`I1-I2`-1.96*sqrt(data2$`I1-I2`*(1-data2$`I1-I2`)/s)
P_I1I2_high <- data2$`I1-I2`+1.96*sqrt(data2$`I1-I2`*(1-data2$`I1-I2`)/s)

P_I2I3_low <- data2$`I2-I3`-1.96*sqrt(data2$`I2-I3`*(1-data2$`I2-I3`)/s)
P_I2I3_high <- data2$`I2-I3`+1.96*sqrt(data2$`I2-I3`*(1-data2$`I2-I3`)/s)

P_I1I3_low <- data2$`I1-I3`-1.96*sqrt(data2$`I1-I3`*(1-data2$`I1-I3`)/s)
P_I1I3_high <- data2$`I1-I3`+1.96*sqrt(data2$`I1-I3`*(1-data2$`I1-I3`)/s)

P_I1I5_low <- data2$`I1-I5`-1.96*sqrt(data2$`I1-I5`*(1-data2$`I1-I5`)/s)
P_I1I5_high <- data2$`I1-I5`+1.96*sqrt(data2$`I1-I5`*(1-data2$`I1-I5`)/s)

P_I1I7_low <- data2$`I1-I7`-1.96*sqrt(data2$`I1-I7`*(1-data2$`I1-I7`)/s)
P_I1I7_high <- data2$`I1-I7`+1.96*sqrt(data2$`I1-I7`*(1-data2$`I1-I7`)/s)

P_I1I7Y_low <- data2$`I1-I7Y`-1.96*sqrt(data2$`I1-I7Y`*(1-data2$`I1-I7Y`)/s)
P_I1I7Y_high <- data2$`I1-I7Y`+1.96*sqrt(data2$`I1-I7Y`*(1-data2$`I1-I7Y`)/s)

LowY <- c(P_I1I7_low,P_I1I7Y_low)
HighY <- c(P_I1I7_high,P_I1I7Y_high)

Low <- c(P_direkte_low,P_2sls_low,P_I1I2_low,P_I2I3_low,P_I1I3_low,P_I1I5_low,P_I1I7_low,P_I1I7Y_low)
High <- c(P_direkte_high,P_2sls_high,P_I1I2_high,P_I2I3_high,P_I1I3_high,P_I1I5_high,P_I1I7_high,P_I1I7Y_high)

Low2sls <- P_2sls_low
High2sls <- P_2sls_high

LowIR <- c(P_I1I2_low,P_I2I3_low,P_I1I3_low,P_I1I5_low,P_I1I7_low)
HighIR <- c(P_I1I2_high,P_I2I3_high,P_I1I3_high,P_I1I5_high,P_I1I7_high)

Low2slsnaive <- c(P_2sls_low,P_direkte_low)
High2slsnaive <- c(P_2sls_high,P_direkte_high)

LowIRnaive <- c(P_I1I3_low,P_I1I7_low,P_direkte_low)
HighIRnaive <- c(P_I1I3_high,P_I1I7_high,P_direkte_high)

Low2slsIR <- c(P_2sls_low,P_I1I3_low,P_I1I7_low)
High2slsIR <- c(P_2sls_high,P_I1I3_high,P_I1I7_high)

LownaiveIR2sls <- c(P_I1I3_low,P_I1I7_low,P_2sls_low,P_direkte_low)
HighnaiveIR2sls <-c(P_I1I3_high,P_I1I7_high,P_2sls_high,P_direkte_high)

#All 
ggplot(data=redata2,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =Low, ymax = High))

#2sls+ naive+ IR
ggplot(data=redatanaiveIR2sls,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =LownaiveIR2sls , ymax = HighnaiveIR2sls))

#2sls
ggplot(data=redata2sls,aes(x=N,y=value))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =Low2sls, ymax = High2sls))

#IR
ggplot(data=redataIR,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =LowIR, ymax = HighIR))

#2sls + naive
ggplot(data=redata2slsnaive,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =Low2slsnaive, ymax = High2slsnaive))

#IR + naive
ggplot(data=redataIRnaive,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =LowIRnaive, ymax = HighIRnaive))

#deltaY
ggplot(data=redataY,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1.1)+geom_errorbar(aes(ymin =LowY, ymax = HighY))

#naive
ggplot(data=redatanaive,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =P_direkte_low, ymax = P_direkte_high))

#2sls + IR
ggplot(data=redata2slsIR,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1.1)+geom_errorbar(aes(ymin =Low2slsIR, ymax = High2slsIR))

