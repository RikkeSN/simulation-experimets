library(quadprog)
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
} #This function generates a data set of size n from the TSIV model with two X variable and two I variable if the matrix 'coeff' is a 6*6 matrix with correct nonzero entries. 
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
}#This function shortens a data vector according to the time lag it represents

#The following two values can be changed to create the desired combination
epsilon <- 0.03 #Epsilon is the bound we consider when investigating consistency 
s <- 200 # s is the number of simulations made pr. matrix A

P_2sls <- numeric(40) #This vector will contain the estimated probability for the different value of the eigenvalue lambda
Xakse <- numeric(40) #This vector will contain the eigenvalue lambda
for (k in 1:20) {
  #The matrix A is the coefficient matrix. We change the smallest eigenvale for alpha; for each k we add k*0.06 times the diagonal in alpha
  A <- matrix(c(0.6,0.2,-0.4,0.3,0,0,0.5,-0.7,-0.8,0.6,0,0,0,0,0.7,0.2,0,0.6,0,0,-0.4,0.7,0,0.3,0,0,0.4,0.6,0.5,-0.4,0,0,0,0,0,0.5),ncol=6)+matrix(c(0,0,k*0.06,0,0,0,0,0,0,k*0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),ncol=6)
  B_2slsX1 <- numeric(s) #These are vectors that for each eigenvalue lambda will contain 's' etatimates from the TS-2SLS estimator
  B_2slsX2 <- numeric(s)
  for (i in 1:s) {
    E <- matrix(rnorm(36000,0,0.1),ncol=6) #All variables are assumed to have standard deviation 0.1
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
    #2SLS - stage 1 
    modelX1 <- lm(X1_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    modelX2 <- lm(X2_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    r_X1 <- coefficients(modelX1)[1]*I1_2+coefficients(modelX1)[2]*I2_2
    r_X2 <- coefficients(modelX2)[1]*I1_2+coefficients(modelX2)[2]*I2_2
    #2SLS stage 2
    model2 <- lm(Y0 ~ 0+ r_X1+r_X2+I1_3+I2_3)
    B_2slsX1[i] = (coefficients(model2)[1])
    B_2slsX2[i] = (coefficients(model2)[2])
  }
  P_2sls[k] <- mean(sqrt((B_2slsX1-A[6,3])^2+(B_2slsX2-A[6,4])^2)>epsilon)  #Here we calculate the estimated probability that the estimate is further from beta than epsilon
  eig <- eigen(A[3:4,1:2]) #Her we calculate the eigenvalues and vectors for alpha
  Xakse[k] <- eig$values[2] #And her we choose the numerically smallest eigenvalue.
}
#We repeat it all again the only difference is that we subtract k*0.06 to get negative eigenvalues.
for (k in 21:40) {
  A <- matrix(c(0.6,0.2,-0.4,0.3,0,0,0.5,-0.7,-0.8,0.6,0,0,0,0,0.7,0.2,0,0.6,0,0,-0.4,0.7,0,0.3,0,0,0.4,0.6,0.5,-0.4,0,0,0,0,0,0.5),ncol=6)+matrix(c(0,0,-(k-20)*0.06,0,0,0,0,0,0,-(k-20)*0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),ncol=6)
  B_2slsX1 <- numeric(s)
  B_2slsX2 <- numeric(s)
  for (i in 1:s) {
    E <- matrix(rnorm(36000,0,0.1),ncol=6)
    sim1data <- rSim(A,E)
    Y0 <- TimeOperator(sim1data[,6],0) 
    Y1 <- TimeOperator(sim1data[,6],1)
    X2_1 <- TimeOperator(sim1data[,4],1)
    X2_2 <- TimeOperator(sim1data[,4],2) 
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
    #2SLS - stage 1 
    modelX1 <- lm(X1_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    modelX2 <- lm(X2_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
    r_X1 <- coefficients(modelX1)[1]*I1_2+coefficients(modelX1)[2]*I2_2
    r_X2 <- coefficients(modelX2)[1]*I1_2+coefficients(modelX2)[2]*I2_2
    #2SLS stage 2
    model2 <- lm(Y0 ~ 0+ r_X1+r_X2+I1_3+I2_3)
    B_2slsX1[i] = (coefficients(model2)[1])
    B_2slsX2[i] = (coefficients(model2)[2])
  }
  P_2sls[k] <- mean(sqrt((B_2slsX1-A[6,3])^2+(B_2slsX2-A[6,4])^2)>epsilon)
  eig <- eigen(A[3:4,1:2])
  Xakse[k] <- eig$values[2]
}


H <- cbind(Xakse,P_2sls) #We collect the list of eigenvalues and the list of estimated probabilities


#Below we create the plot
library(tidyverse)
library(reshape2)
library(ggplot2)
library(tikzDevice)
data2 <- as.data.frame(H)
names(data2) <- c("eigen","2SLS")
redata2 <- reshape2::melt(data2, id.var = "eigen")


P_2sls_low <- data2$`2SLS`-1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)
P_2sls_high <- data2$`2SLS`+1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)

ggplot(data=redata2,aes(x=eigen,y=value))+
  geom_point() + 
  geom_line() + labs(x="Eigenvalue", y="Est. probability")+
  ggtitle("Estimated probability")+ylim(0,1.1)+geom_errorbar(aes(ymin =P_2sls_low, ymax = P_2sls_high))
