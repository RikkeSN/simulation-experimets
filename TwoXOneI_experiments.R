library(quadprog)
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
}#This function generates a data set of size n from the TSIV model with two X variables and one I variable if the matrix 'coeff' is a 5*5 matrix with correct nonzero entries. 
#The input 'error' is an n*5 matrix containing the error terms

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
epsilon <- 0.05  #Epsilon is the bound we consider when investigating consistency 
f <- 10 #f must be divisible by 5. It  is a parameter that influences the number of time steps, n, included in the simulated data
s <- 100 # s is the number of simulations made pr. number of time step
A <- matrix(c(0.6,-0.4,0.6,0,0,0,0.7,0.2,0,0.6,0,-0.4,0.7,0,0.3,0,0.4,0.6,0.5,-0.4,0,0,0,0,0.5),ncol=5) #The data generating coefficient matrix from model M3


P_covI7 <- numeric(10) #This vector will contain the estimated probability for the TS-IR I1-I7 estimator
S <- c(2100,4200,6300,8400,10500,12600,14700,16800,18900,21000)*f #S defines the different number of time steps that will be used
for (j in S) {
  B_covI7X1 <- numeric(s) #These are vectors that for each time step in S will contain 's' etatimates from the given estimator
  B_covI7X2 <- numeric(s)
  for (i in 1:s) {
    E <- matrix(rnorm(j,0,0.1),ncol=5) #All variables are assumed to have standard deviation 0.1
    sim1data <- rSim(A,E) #Here we simulate a data set.
    Y0 <- TimeOperator(sim1data[,5],0) 
    Y1 <- TimeOperator(sim1data[,5],1)
    X2_1 <- TimeOperator(sim1data[,3],1)
    X2_2 <- TimeOperator(sim1data[,3],2) #X1_2 means X1 at time t-2
    X2_3 <- TimeOperator(sim1data[,3],3)
    I1_0 <- TimeOperator(sim1data[,1],0)
    I1_1 <- TimeOperator(sim1data[,1],1)
    I1_2 <- TimeOperator(sim1data[,1],2) 
    I1_3 <- TimeOperator(sim1data[,1],3)
    I1_4 <- TimeOperator(sim1data[,1],4)
    I1_5 <- TimeOperator(sim1data[,1],5)
    I1_6 <- TimeOperator(sim1data[,1],6)
    I1_7 <- TimeOperator(sim1data[,1],7)
    X1_1 <- TimeOperator(sim1data[,2],1)
    X1_2 <- TimeOperator(sim1data[,2],2)
    X1_3 <- TimeOperator(sim1data[,2],3)
    X1_4 <- TimeOperator(sim1data[,2],4)
    #TS-IR estimate including I1-I7
    Cov2_matrixI4 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I1_2),cov(X1_1,I1_3),cov(X1_1,I1_4),cov(X1_1,I1_5),cov(X1_1,I1_6),cov(X1_1,I1_7),cov(X2_1,I1_1),cov(X2_1,I1_2),cov(X2_1,I1_3),cov(X2_1,I1_4),cov(X2_1,I1_5),cov(X2_1,I1_6),cov(X2_1,I1_7),cov(Y1,I1_1),cov(Y1,I1_2),cov(Y1,I1_3),cov(Y1,I1_4),cov(Y1,I1_5),cov(Y1,I1_6),cov(Y1,I1_7)),ncol=3)
    Cov2_vectorI4 <- c(cov(Y0,I1_1),cov(Y0,I1_2),cov(Y0,I1_3),cov(Y0,I1_4),cov(Y0,I1_5),cov(Y0,I1_6),cov(Y0,I1_7))
    cov2_est_vectorI4 <- qr.solve(Cov2_matrixI4,Cov2_vectorI4)
    B_covI7X1[i] <-cov2_est_vectorI4[1] 
    B_covI7X2[i] <-cov2_est_vectorI4[2] 
  }
  k <- j/(2100*f) #Creating an index from 1 to 14
  P_covI7[k] <- mean(sqrt((B_covI7X1-A[5,2])^2+(B_covI7X2-A[5,3])^2)>epsilon) #Here we calculate the estimated probability that the estimate is further from beta than epsilon
}

N=c(2100,4200,6300,8400,10500,12600,14700,16800,18900,21000)*f/5 #N is a list of the number of time steps used. We devide by 5 since the error matrix has dimension n*5
H <- cbind(N,P_covI7) #Here we collect the estimated probability

#Below we create the fiugres
library(tidyverse)
library(reshape2)
library(ggplot2)
library(tikzDevice)
data2 <- as.data.frame(H)
names(data2) <- c("N","I1-I7")
redata2<- reshape2::melt(data2, id.var = "N")


P_cov_low <- data2$`I1-I7`-1.96*sqrt(data2$`I1-I7`*(1-data2$`I1-I7`)/s)
P_cov_high <- data2$`I1-I7`+1.96*sqrt(data2$`I1-I7`*(1-data2$`I1-I7`)/s)

ggplot(data=redata2,aes(x=N,y=value,color=variable))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability", color="Estimates")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =P_cov_low, ymax = P_cov_high))

