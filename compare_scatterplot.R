rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
}#This function generates a data set of size n from the one dimensional TSIV model if the matrix 'coeff' is a 4*4 matrix with correct nonzero entries. 
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

#The following two values can be changed to create the desired combination
epsilon <- 0.05  #Epsilon is the bound we consider when investigating consistency 
s <- 100     # s is the number of simulations made pr.random marrix
n <- 4000 #The number of time steps included in the generated date is n/4

P_2sls <- numeric(100) #These are vectors that will contain the estimated probability for the different estimators
P_covI1I7 <- numeric(100)

for (k in 1:100 ) {
  B_2sls <- numeric(s) #These are vectors that for each random matrix will contain 's' etatimates from the given estimator
  B_covI1I7X <- numeric(s)
  A <- matrix(c(runif(1,-0.99,0.99),runif(1,-2,2),0,0,0,runif(1,-0.99,0.99),0,runif(1,-2,2),0,runif(1,-2,2),runif(1,-0.99,0.99),runif(1,-2,2),0,0,0,runif(1,-0.99,0.99)),ncol=4) #The data generating coefficient matrix chosen randomly
  for (i in 1:s) {
    E <- matrix(rnorm(n,0,0.1),ncol=4)#All variables are assumed to have standard deviation 0.1
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
    #TS-2SLS - stage 1 
    model <- lm(X1 ~ 0+I2+I3)
    r <- coefficients(model)[1]*I2
    #TS-2SLS stage 2
    model2 <- lm(Y0 ~ 0+r+I3)
    B_2sls[i] = (coefficients(model2)[1])
    #TS-IR estimate I1 to I7
    Cov2_matrixI1I7 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(X1,I4),cov(X1,I5),cov(X1,I6),cov(X1,I7),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3),cov(Y1,I4),cov(Y1,I5),cov(Y1,I6),cov(Y1,I7)),ncol=2)
    Cov2_vectorI1I7 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3),cov(Y0,I4),cov(Y0,I5),cov(Y0,I6),cov(Y0,I7))
    cov2_est_vectorI1I7 <- qr.solve(Cov2_matrixI1I7,Cov2_vectorI1I7)
    B_covI1I7X[i] <-cov2_est_vectorI1I7[1] 
  }
  P_2sls[k] <-mean(abs(B_2sls-A[4,2])>epsilon) #Here we calculate the estimated probability that the estimate is further from beta than epsilon
  P_covI1I7[k] <- mean(abs(B_covI1I7X-A[4,2])>epsilon)
}

H <- cbind(P_2sls,P_covI1I7) #Here we collect the estimated probabilities

#Below we create the figure
library(tidyverse)
library(reshape2)
library(ggplot2)
library(tikzDevice)
data2 <- as.data.frame(H)
names(data2) <- c("tsls","I1-I7")
redata2 <- reshape2::melt(data2, id.var = "tsls")

ggplot(data=redata2,aes(x=tsls,y=value))+
  geom_point() + labs(x="est. prob. TS-2SLS", y="est. prob. TS-IR")+
  ggtitle("Scatter plot comparing TS-2SLS and TS-IR")+geom_abline(intercept=0,slope=1)

