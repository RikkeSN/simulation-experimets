#In this file are two functions. The function 'Estimate_1D' can calculate the TS_2SLS and the TS-IR estimate for beta if X is one-dimensional and I is one-dimensional
# The function 'Estimate_2D' can calculate the TS-2SLS and the TS-IR estimate for beta if X is two-dimensional and I is two-dimensional


#The first function:
#The input 'data' must be a data matrix from a time series with three columns; I,X, and Y in that order
Estimate_1D <- function(data) {
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
  }
  Y0 <- TimeOperator(data[,3],0) #Here the data vectors are shortened according to their time lag. 
  Y1 <- TimeOperator(data[,3],1)
  X1 <- TimeOperator(data[,2],1) 
  X2 <- TimeOperator(data[,2],2) 
  X3 <- TimeOperator(data[,2],3)
  I1 <- TimeOperator(data[,1],1)
  I2 <- TimeOperator(data[,1],2) 
  I3 <- TimeOperator(data[,1],3)
  I4 <- TimeOperator(data[,1],4)
  I5 <- TimeOperator(data[,1],5)
  I6 <- TimeOperator(data[,1],6)
  I7 <- TimeOperator(data[,1],7)
  
  #TS-2SLS - stage 1 
  model <- lm(X1 ~ 0+I2+I3)
  r <- coefficients(model)[1]*I2
  #TS-2SLS stage 2
  model2 <- lm(Y0 ~ 0+r+I3)
  B_2sls <- (coefficients(model2)[1])
  
  #TS-IR estimate including I1 to I3
  Cov2_matrixI1I3 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3)),ncol=2)
  Cov2_vectorI1I3 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3))
  cov2_est_vectorI1I3 <- qr.solve(Cov2_matrixI1I3,Cov2_vectorI1I3)
  B_covI1I3X <-cov2_est_vectorI1I3[1] 
  
  #TS-IR estimate ID I1 to I7
  Cov2_matrixI1I7 <- matrix(c(cov(X1,I1),cov(X1,I2),cov(X1,I3),cov(X1,I4),cov(X1,I5),cov(X1,I6),cov(X1,I7),cov(Y1,I1),cov(Y1,I2),cov(Y1,I3),cov(Y1,I4),cov(Y1,I5),cov(Y1,I6),cov(Y1,I7)),ncol=2)
  Cov2_vectorI1I7 <- c(cov(Y0,I1),cov(Y0,I2),cov(Y0,I3),cov(Y0,I4),cov(Y0,I5),cov(Y0,I6),cov(Y0,I7))
  cov2_est_vectorI1I7 <- qr.solve(Cov2_matrixI1I7,Cov2_vectorI1I7)
  B_covI1I7X <-cov2_est_vectorI1I7[1] 
  
  Beta<- c(B_2sls,B_covI1I3X,B_covI1I7X)
  names <- c("beta_2sls","beta_TS-IR_I1-I3", "beta_TS-IR_I1-I7")
  beta <- cbind(names,Beta)
  return(beta)
}
#The function returns three numbers:
#The first is the estimate for beta based on the TS-2SLS estimator
#The second is the estimate for beta based on the TS-IR estimator I1-I3
#The third is the estimate for beta based on the TS-IR estimator I1-I7

#Testing the function with simulated data:
rSim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
} 
A <- matrix(c(0.5,-0.9,0,0,0,0.5,0,-0.5,0,-0.8,-0.8,-0.8,0,0,0,-0.5),ncol=4) 
E <- matrix(rnorm(400000,0,0.1),ncol=4) #All variables are assumed to have standard deviation 0.1
sim1data <- rSim(A,E) 
sim1data_ny <- sim1data[,-3]
Estimate_1D(sim1data_ny)

###################################

#The second function 'Estimate_2D'
#The input 'data' must be a data matrix from a time series with 5 columns; I1, I2, X1, X2, and Y in that order

Estimate_2D <- function(data) {
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
  }
  Y0 <- TimeOperator(data[,5],0) 
  Y1 <- TimeOperator(data[,5],1)
  X2_1 <- TimeOperator(data[,4],1)
  X2_2 <- TimeOperator(data[,4],2) 
  X2_3 <- TimeOperator(data[,4],3)
  I1_1 <- TimeOperator(data[,1],1)
  I1_2 <- TimeOperator(data[,1],2) 
  I2_1 <- TimeOperator(data[,2],1)
  I2_2 <- TimeOperator(data[,2],2)
  I1_3 <- TimeOperator(data[,1],3)
  I2_3 <- TimeOperator(data[,2],3)
  I1_4 <- TimeOperator(data[,1],4)
  I2_4 <- TimeOperator(data[,2],4)
  I1_5 <- TimeOperator(data[,1],5)
  I2_5 <- TimeOperator(data[,2],5)
  I1_6 <- TimeOperator(data[,1],6)
  I2_6 <- TimeOperator(data[,2],6)
  I1_7 <- TimeOperator(data[,1],7)
  I2_7 <- TimeOperator(data[,2],7)
  X1_1 <- TimeOperator(data[,3],1)
  X1_2 <- TimeOperator(data[,3],2)
  X1_3 <- TimeOperator(data[,3],3)
  X1_4 <- TimeOperator(data[,3],4)
  #Overidentified estimate including I1 - I3
  Cov2_matrixI1I3 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3)),ncol=3)
  Cov2_vectorI1I3 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3))
  cov2_est_vectorI1I3 <- qr.solve(Cov2_matrixI1I3,Cov2_vectorI1I3)
  B_covI1I3X1 <-cov2_est_vectorI1I3[1] 
  B_covI1I3X2 <-cov2_est_vectorI1I3[2] 
  #Overidentified estimate including  I1- I7
  Cov2_matrixI1I7 <- matrix(c(cov(X1_1,I1_1),cov(X1_1,I2_1),cov(X1_1,I1_2),cov(X1_1,I2_2),cov(X1_1,I1_3),cov(X1_1,I2_3),cov(X1_1,I1_4),cov(X1_1,I2_4),cov(X1_1,I1_5),cov(X1_1,I2_5),cov(X1_1,I1_6),cov(X1_1,I2_6),cov(X1_1,I1_7),cov(X1_1,I2_7),cov(X2_1,I1_1),cov(X2_1,I2_1),cov(X2_1,I1_2),cov(X2_1,I2_2),cov(X2_1,I1_3),cov(X2_1,I2_3),cov(X2_1,I1_4),cov(X2_1,I2_4),cov(X2_1,I1_5),cov(X2_1,I2_5),cov(X2_1,I1_6),cov(X2_1,I2_6),cov(X2_1,I1_7),cov(X2_1,I2_7),cov(Y1,I1_1),cov(Y1,I2_1),cov(Y1,I1_2),cov(Y1,I2_2),cov(Y1,I1_3),cov(Y1,I2_3),cov(Y1,I1_4),cov(Y1,I2_4),cov(Y1,I1_5),cov(Y1,I2_5),cov(Y1,I1_6),cov(Y1,I2_6),cov(Y1,I1_7),cov(Y1,I2_7)),ncol=3)
  Cov2_vectorI1I7 <- c(cov(Y0,I1_1),cov(Y0,I2_1),cov(Y0,I1_2),cov(Y0,I2_2),cov(Y0,I1_3),cov(Y0,I2_3),cov(Y0,I1_4),cov(Y0,I2_4),cov(Y0,I1_5),cov(Y0,I2_5),cov(Y0,I1_6),cov(Y0,I2_6),cov(Y0,I1_7),cov(Y0,I2_7))
  cov2_est_vectorI1I7 <- qr.solve(Cov2_matrixI1I7,Cov2_vectorI1I7)
  B_covI1I7X1 <-cov2_est_vectorI1I7[1] 
  B_covI1I7X2 <-cov2_est_vectorI1I7[2] 
  #2SLS - stage 1 
  modelX1 <- lm(X1_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
  modelX2 <- lm(X2_1 ~ 0+I1_2+I2_2+I1_3+I2_3)
  r_X1 <- coefficients(modelX1)[1]*I1_2+coefficients(modelX1)[2]*I2_2
  r_X2 <- coefficients(modelX2)[1]*I1_2+coefficients(modelX2)[2]*I2_2
  #2SLS stage 2
  model2 <- lm(Y0 ~ 0+ r_X1+r_X2+I1_3+I2_3)
  B_2slsX1 = (coefficients(model2)[1])
  B_2slsX2 = (coefficients(model2)[2])
  
  Beta_1 <- c(B_2slsX1,B_covI1I3X1,B_covI1I7X1)
  Beta_2 <- c(B_2slsX2,B_covI1I3X2,B_covI1I7X2)
  names <- c("beta_2sls","beta_TS-IR_I1-I3", "beta_TS-IR_I1-I7")
  beta <- cbind(names,Beta_1,Beta_2)
  return(beta)
}

#The function returns three numbers for beta_1 and three numbers for beta_2:
#The first is the estimate for beta_i based on the TS-2SLS estimator
#The second is the estimate for beta_i based on the TS-IR estimator I1-I3
#The third is the estimate for beta_i based on the TS-IR estimator I1-I7


#Testing the function:
A <- matrix(c(0.6,0.2,-0.4,0.6,0,0,0.5,-0.7,-0.8,-0.5,0,0,0,0,-0.7,0.2,0,0.6,0,0,-0.4,0.7,0,0.3,0,0,0.4,0.6,0.5,-0.4,0,0,0,0,0,0.5),ncol=6)
E <- matrix(rnorm(600000,0,0.1),ncol=6)
sim1data <- rSim(A,E)
sim1data_2 <- sim1data[,-5]
Estimate_2D(sim1data_2)
