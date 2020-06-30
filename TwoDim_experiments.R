library(quadprog)
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
} #This function generates a data set of size n from the TSIV model with two X variables and two I variables if the matrix 'coeff' is a 6*6 matrix with correct nonzero entries. 
#The input 'error' is an n*6 matrix containing the error terms

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
epsilon <- 0.01 #Epsilon is the bound we consider when investigating consistency 
f<- 12           #f must be divisible by 6. It  is a parameter that influences the number of time steps, n, included in the simulated data
s <- 100        # s is the number of simulations made pr. number of time step
A <- matrix(c(0.6,0.2,-0.4,0.6,0,0,0.5,-0.7,-0.8,-0.5,0,0,0,0,0.7,0.2,0,0.6,0,0,-0.4,0.7,0,0.3,0,0,0.4,0.6,0.5,-0.4,0,0,0,0,0,0.5),ncol=6) #The data generating coefficient matrix. This one is the one from model M2


P_covI4 <- numeric(10) #These are vectors that will contain the estimated probability for the different estimators
P_covI2I3 <- numeric(10)
P_2sls <- numeric(10)
P_direkte <- numeric(10)
P_covI1I3 <- numeric(10)
P_covI1I2 <- numeric(10)
P_covI1I4 <- numeric(10)
P_covI1I7 <- numeric(10)
S <- c(2100,4200,6300,8400,10500,12600,14700,16800,18900,21000)*f #S defines the different number of time steps that will be used
for (j in S) {
  B_direkteX1 <- numeric(n) #These are vectors that for each time step in S will contain 's' etatimates from the given estimator
  B_direkteX2 <- numeric(n)
  B_covX1 <- numeric(n)
  B_covX2 <- numeric(n)
  B_covnyX1 <- numeric(n)
  B_covnyX2 <- numeric(n)
  B_2slsX1 <- numeric(n)
  B_2slsX2 <- numeric(n)
  B_covI4X1 <- numeric(n)
  B_covI4X2 <- numeric(n)
  B_covI1I3X1 <- numeric(n)
  B_covI1I3X2 <- numeric(n)
  B_covI1I4X1 <- numeric(n)
  B_covI1I4X2 <- numeric(n)
  B_covI1I7X1 <- numeric(n)
  B_covI1I7X2 <- numeric(n)
  B_covI1I2X1 <- numeric(n)
  B_covI1I2X2 <- numeric(n)
  for (i in 1:n) {
    E <- matrix(rnorm(j,0,0.1),ncol=6) #All variables are assumed to have standard deviation 0.1
    sim1data <- rSim(A,E) #Here we simulate a data set.
    Y0 <- TimeOperator(sim1data[,6],0) #Here the data vectors are shortened according to their time lag.
    Y1 <- TimeOperator(sim1data[,6],1)
    X2_1 <- TimeOperator(sim1data[,4],1)
    X2_2 <- TimeOperator(sim1data[,4],2) #X1_2 means X1 at time t-2
    X2_3 <- TimeOperator(sim1data[,4],3)
    I1_1 <- TimeOperator(sim1data[,1],1)
    I1_2 <- TimeOperator(sim1data[,1],2) 
    I2_1 <- TimeOperator(sim1data[,2],1)
    I2_2 <- TimeOperator(sim1data[,2],2)
    I1_3 <- TimeOperator(sim1data[,1],3)
    I2_3 <- TimeOperator(sim1data[,2],3)
    I1_4 <- TimeOperator(sim1data[,1],4)
    I2_4 <- TimeOperator(sim1data[,2],4)
    I1_5 <- TimeOperator(sim1data[,1],5)
    I2_5 <- TimeOperator(sim1data[,2],5)
    I1_6 <- TimeOperator(sim1data[,1],6)
    I2_6 <- TimeOperator(sim1data[,2],6)
    I1_7 <- TimeOperator(sim1data[,1],7)
    I2_7 <- TimeOperator(sim1data[,2],7)
    X1_1 <- TimeOperator(sim1data[,3],1)
    X1_2 <- TimeOperator(sim1data[,3],2)
    X1_3 <- TimeOperator(sim1data[,3],3)
    X1_4 <- TimeOperator(sim1data[,3],4)
    #TS-IR estimate including I1 - I2
    Cov2_matrixI1I2 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2)),ncol=3)
    Cov2_vectorI1I2 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2))
    cov2_est_vectorI1I2 <- qr.solve(Cov2_matrixI1I2,Cov2_vectorI1I2)
    B_covI1I2X1[i] <-cov2_est_vectorI1I2[1] 
    B_covI1I2X2[i] <-cov2_est_vectorI1I2[2]
    #TS-IR estimate including I2 - I3
    Cov2_matrixny <- matrix(c(cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3)),ncol=3)
    Cov2_vectorny <- c(cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3))
    cov2_est_vectorny <- qr.solve(Cov2_matrixny,Cov2_vectorny)
    B_covnyX1[i] <-cov2_est_vectorny[1] 
    B_covnyX2[i] <-cov2_est_vectorny[2] 
    #TS-IR estimate including I1 - I3
    Cov2_matrixI1I3 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3)),ncol=3)
    Cov2_vectorI1I3 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3))
    cov2_est_vectorI1I3 <- qr.solve(Cov2_matrixI1I3,Cov2_vectorI1I3)
    B_covI1I3X1[i] <-cov2_est_vectorI1I3[1] 
    B_covI1I3X2[i] <-cov2_est_vectorI1I3[2] 
    #TS-IR estimate including  I2- I4
    Cov2_matrixI4 <- matrix(c(cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X1_1,I1_4),cov(X1_1,I2_4),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(X2_1,I1_4),cov(X2_1,I2_4),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3),cov(Y1,I1_4),cov(Y1,I2_4)),ncol=3)
    Cov2_vectorI4 <- c(cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3),cov(Y0,I1_4),cov(Y0,I2_4))
    cov2_est_vectorI4 <- qr.solve(Cov2_matrixI4,Cov2_vectorI4)
    B_covI4X1[i] <-cov2_est_vectorI4[1] 
    B_covI4X2[i] <-cov2_est_vectorI4[2] 
    #TS-IR estimate including  I1- I4
    Cov2_matrixI1I4 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X1_1,I1_4),cov(X1_1,I2_4),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(X2_1,I1_4),cov(X2_1,I2_4),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3),cov(Y1,I1_4),cov(Y1,I2_4)),ncol=3)
    Cov2_vectorI1I4 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3),cov(Y0,I1_4),cov(Y0,I2_4))
    cov2_est_vectorI1I4 <- qr.solve(Cov2_matrixI1I4,Cov2_vectorI1I4)
    B_covI1I4X1[i] <-cov2_est_vectorI1I4[1] 
    B_covI1I4X2[i] <-cov2_est_vectorI1I4[2] 
    #TS-IR estimate including  I1- I7
    Cov2_matrixI1I7 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X1_1,I1_4),cov(X1_1,I2_4),cov(X1_1,I1_5),cov(X1_1,I2_5),cov(X1_1,I1_6),cov(X1_1,I2_6),cov(X1_1,I1_7),cov(X1_1,I2_7),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(X2_1,I1_4),cov(X2_1,I2_4),cov(X2_1,I1_5),cov(X2_1,I2_5),cov(X2_1,I1_6),cov(X2_1,I2_6),cov(X2_1,I1_7),cov(X2_1,I2_7),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3),cov(Y1,I1_4),cov(Y1,I2_4),cov(Y1,I1_5),cov(Y1,I2_5),cov(Y1,I1_6),cov(Y1,I2_6),cov(Y1,I1_7),cov(Y1,I2_7)),ncol=3)
    Cov2_vectorI1I7 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3),cov(Y0,I1_4),cov(Y0,I2_4),cov(Y0,I1_5),cov(Y0,I2_5),cov(Y0,I1_6),cov(Y0,I2_6),cov(Y0,I1_7),cov(Y0,I2_7))
    cov2_est_vectorI1I7 <- qr.solve(Cov2_matrixI1I7,Cov2_vectorI1I7)
    B_covI1I7X1[i] <-cov2_est_vectorI1I7[1] 
    B_covI1I7X2[i] <-cov2_est_vectorI1I7[2] 
    #TS-2SLS - stage 1 
    modelX1 <- lm(X1_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    modelX2 <- lm(X2_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    r_X1 <- coefficients(modelX1)[1]*I1_2+coefficients(modelX1)[2]*I2_2
    r_X2 <- coefficients(modelX2)[1]*I1_2+coefficients(modelX2)[2]*I2_2
    #TS-2SLS stage 2
    model2 <- lm(Y0 ~ 0+ r_X1+r_X2+I1_3+I2_3)
    B_2slsX1[i] = (coefficients(model2)[1])
    B_2slsX2[i] = (coefficients(model2)[2])
    #Naive estimator
    model3 <- lm(Y0~0+X1_1+X2_1+Y1)
    B_direkteX1[i] <- coefficients(model3)[1]
    B_direkteX2[i] <- coefficients(model3)[2]
  }
  k <- j/(2100*f) #creating an index from 1 to 14
  P_covI2I3[k] <- mean(sqrt((B_covnyX1-A[6,3])^2+(B_covnyX2-A[6,4])^2)>epsilon) #Here we calculate the estimated probability that the estimate is further from beta than epsilon
  P_covI1I3[k] <- mean(sqrt((B_covI1I3X1-A[6,3])^2+(B_covI1I3X2-A[6,4])^2)>epsilon)
  P_covI1I2[k] <- mean(sqrt((B_covI1I2X1-A[6,3])^2+(B_covI1I2X2-A[6,4])^2)>epsilon)
  P_covI4[k] <- mean(sqrt((B_covI4X1-A[6,3])^2+(B_covI4X2-A[6,4])^2)>epsilon)
  P_covI1I4[k] <- mean(sqrt((B_covI1I4X1-A[6,3])^2+(B_covI1I4X2-A[6,4])^2)>epsilon)
  P_covI1I7[k] <- mean(sqrt((B_covI1I7X1-A[6,3])^2+(B_covI1I7X2-A[6,4])^2)>epsilon)
  P_2sls[k] <- mean(sqrt((B_2slsX1-A[6,3])^2+(B_2slsX2-A[6,4])^2)>epsilon)
  P_direkte[k] <- mean(sqrt((B_direkteX1-A[6,3])^2+(B_direkteX2-A[6,4])^2)>epsilon)
}

N=c(2100,4200,6300,8400,10500,12600,14700,16800,18900,21000)*f/6 #N is a list of the number of time steps used. We devide by 6 since the error matrix has dimension n*6
H <- cbind(N,P_direkte,P_2sls,P_covI1I2,P_covI2I3,P_covI1I3, P_covI1I4,P_covI1I7) #Here we collect all the estimated probabilities

#Below we create the different figures
library(tidyverse)
library(reshape2)
library(ggplot2)
library(tikzDevice)
data2 <- as.data.frame(H)
names(data2) <- c("N","Naive", "2SLS", "I1-I2","I2-I3","I1-I3","I1-I5","I1-I7")
data2sls <- data.frame(data2$N,data2$`2SLS`)
dataIR <- data.frame(data2$N,data2$`I1-I2`,data2$`I2-I3`,data2$`I1-I3`,data2$`I1-I5`,data2$`I1-I7`)
data2slsnaive <- data.frame(data2$N,data2$`2SLS`,data2$Naive)
dataIRnaive <- data.frame(data2$N,data2$`I1-I3`,data2$`I1-I7`,data2$Naive)
data2slsIR <- data.frame(data2$N,data2$`2SLS`, data2$`I1-I3`,data2$`I1-I7`)
datanaiveIR2sls <- data.frame(data2$N,data2$`I1-I3`,data2$`I1-I7`,data2$`2SLS`,data2$Naive)
datanaive <- data.frame(data2$N,data2$Naive)



names(data2sls) <-c("N", "2SLS") 
names(dataIR) <- c("N", "I1-I2","I2-I3","I1-I3","I1-I5","I1-I7")
names(data2slsnaive) <- c("N", "2SLS","Naive")
names(dataIRnaive) <- c("N","I1-I3","I1-I7","Naive")
names(data2slsIR) <- c("N", "2SLS","I1-I3","I1-I7")
names(datanaiveIR2sls) <- c("N","I1-I3","I1-I7","2SLS","Naive")
names(datanaive) <- c("N","Naive")


redata2 <- reshape2::melt(data2, id.var = "N")
redata2sls <- reshape2::melt(data2sls, id.var = "N")
redataIR <- reshape2::melt(dataIR, id.var = "N")
redata2slsnaive <- reshape2::melt(data2slsnaive, id.var = "N")
redataIRnaive <- reshape2::melt(dataIRnaive, id.var = "N")
redata2slsIR <- reshape2::melt(data2slsIR, id.var = "N")
redatanaiveIR2sls <- reshape2::melt(datanaiveIR2sls, id.var = "N")
redatanaive <- reshape2::melt(datanaive, id.var = "N")


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




Low <- c(P_direkte_low,P_2sls_low,P_I1I2_low,P_I2I3_low,P_I1I3_low,P_I1I5_low,P_I1I7_low)
High <- c(P_direkte_high,P_2sls_high,P_I1I2_high,P_I2I3_high,P_I1I3_high,P_I1I5_high,P_I1I7_high)

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

