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
}#This function shortens a data vector according to the time lag it represents

#In this first section we create the histograms 
s <- 100 # s is the number of simulations

B_2sls <- numeric(s) #These are vectors that will contain 's' etatimates from the given estimator
B_direkte <- numeric(s)
for (i in 1:s) {
  A <-  matrix(c(0.5,0.5,0,0,0,0.5,0,0.5,0,0.5,0.5,0.5,0,0,0,0.5),ncol=4) #The data generating model with all nonzero variables being 0.5
  E <- matrix(rnorm(22400,0,0.1),ncol=4) #All error terms have standard deviation 0.1
  sim1data <- rSim(A,E) #Here we estimate a data set
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
  #Naive estimator
  model3 <- lm(Y0~0+X1+Y1)
  B_direkte[i] <- coefficients(model3)[1]
  #2SLS estimator - stage 1 
  model <- lm(X1 ~ 0+I2)
  r <- coefficients(model)[1]*I2
  #2SLS estimator stage 2
  model2 <- lm(Y0 ~ 0+r)
  B_2sls[i] = (coefficients(model2)[1])
}

data_2sls <- data.frame(
  weight=B_2sls)
data_naive<- data.frame(
  weight=B_direkte)

ggplot(data_2sls, aes(x=weight)) + geom_histogram(color="blue", fill="lightblue")+geom_vline(aes(xintercept=0.5),color="blue", linetype="dashed", size=1)+geom_vline(aes(xintercept=0.667),
                                                                                                                                                                   color="red", linetype="dashed", size=1)+
  ggtitle("histogram of the intuitive 2sls estimator")

ggplot(data_naive, aes(x=weight)) + geom_histogram(color="blue", fill="lightblue")+geom_vline(aes(xintercept=0.5),color="blue", linetype="dashed", size=1)+geom_vline(aes(xintercept=0.549),
                                                                                                                                                                     color="red", linetype="dashed", size=1)+
  ggtitle("histogram of the naive estimator")

##########################


#In this section we create the estimated probability plots
epsilon <- 0.01  #Epsilon is the bound we consider when investigating consistency 
f <- 7 #f is a parameter that influences the number of time steps, n, included in the simulated data
s <- 100 # s is the number of simulations made pr. number of time step


P_2sls <- numeric(14) #These vectors will contain the estimated probability for the different estimators
P_direkte <- numeric(14)
S <-c(1600,3200,4800,6400,8000,9600,11200,12800,14400,16000,17600,19200,20800,22400)*f #S defines the different number of time steps that will be used
for (j in S ) {
  B_2sls <- numeric(s) #These are vectors that for each time step in S will contain 's' etatimates from the given estimator
  B_direkte <- numeric(s)
  for (i in 1:s) {
    A <-  matrix(c(0.5,0.5,0,0,0,0.5,0,0.5,0,0.5,0.5,0.5,0,0,0,0.5),ncol=4) #The data generating coefficient matrix with all nonzero values 0.5
    E <- matrix(rnorm(j,0,0.1),ncol=4) #All error terms have standard deviation 0.1
    sim1data <- rSim(A,E) #Here we simulate data
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
    #2SLS estimator - stage 1 
    model <- lm(X1 ~ 0+I2)
    r <- coefficients(model)[1]*I2
    #2SLS estimator stage 2
    model2 <- lm(Y0 ~ 0+r)
    B_2sls[i] = (coefficients(model2)[1])
    #Naive estimator
    model3 <- lm(Y0~0+X1+Y1)
    B_direkte[i] <- coefficients(model3)[1]
  }
  k <- j/(1600*f) #Creating an index from 1 to 14
  P_2sls[k] <-mean(abs(B_2sls-0.667)>epsilon)  #Here we calculate the estimated probability that the estimate is further from beta than epsilon
  P_direkte[k] <- mean(abs(B_direkte-0.549)>epsilon)
}
N <- c(1600,3200,4800,6400,8000,9600,11200,12800,14400,16000,17600,19200,20800,22400)/4*f #N is a list of the number of time steps used. We devide by 4 since the error matrix has dimension n*4
H <- cbind(N,P_direkte,P_2sls) #Collecting the number of time steps and the estimated probability vectors


#below we create the plots
data2 <- as.data.frame(H)
names(data2) <- c("N","Naive", "2SLS")
data2sls <- data.frame(data2$N,data2$`2SLS`)
datanaive <- data.frame(data2$N,data2$Naive)


names(data2sls) <-c("N", "2SLS") 
names(datanaive) <- c("N","Naive")

redata2 <- reshape2::melt(data2, id.var = "N")
redata2sls <- reshape2::melt(data2sls, id.var = "N")
redatanaive <- reshape2::melt(datanaive, id.var = "N")

P_direkte_low <- data2$Naive-1.96*sqrt(data2$Naive*(1-data2$Naive)/s)
P_direkte_high <- data2$Naive+1.96*sqrt(data2$Naive*(1-data2$Naive)/s)

P_2sls_low <- data2$`2SLS`-1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)
P_2sls_high <- data2$`2SLS`+1.96*sqrt(data2$`2SLS`*(1-data2$`2SLS`)/s)

#2sls
 ggplot(data=redata2sls,aes(x=N,y=value))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability")+
  ggtitle("Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =P_2sls_low, ymax = P_2sls_high))

#naive
ggplot(data=redatanaive,aes(x=N,y=value))+
  geom_point() + 
  geom_line() + labs(x="Number of time steps", y="Est. probability")+
  ggtitle("1Estimated probability")+ylim(0,1)+geom_errorbar(aes(ymin =P_direkte_low, ymax = P_direkte_high))

