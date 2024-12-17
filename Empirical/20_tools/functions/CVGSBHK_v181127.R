cvgs.bhk2015 <- function(Y,ZZ,k,lambdavec,lambda2=0.001,dimX,plot=0,sweigths=1,eweigths=1,nf=dimX){
  #what takes time is to compute the K matrix each time. Using the FWL trick 
  #and having the k loop being the outter loop rather than lambda (as GA) function, 
  #we can speed things up dramatically. However, CV interpretation will change. 
  #Instead of finding the best model, it will find the best improvment wrt to constant 
  #paramater VAR.
  
  set.seed(1071) # Seed for replicability
  
  
  # k - number of folds
  Y<- as.matrix(Y)
  T<-dim(ZZ)[1]
  data <- cbind(Y,ZZ)
  dfTrain<- data
  
  index <- sample(1:k, nrow(data), replace = TRUE)
  
  #randomseq <- sample(nrow(data))
  #dfTrain <- data[randomseq,]
  #folds <- cut( seq(1,nrow(data)), breaks = k, labels=FALSE)
  
  seqY = 1:dim(Y)[2] #for MV model
  PMSE <- array(0,dim=c(k,length(lambdavec)))
  
  # Do the cross validation for K fold
  for(i in 1:k){
    #Segment your data by fold using the which() function 
    testIndexes <- which(index == i, arr.ind=TRUE)
    testData <- dfTrain[testIndexes, ]
    trainData <- dfTrain[-testIndexes, ]
    
    Lambda_T<-diag(nrow(trainData))
    if(length(eweigths)>1){Lambda_T<- diag(as.numeric((eweigths[-testIndexes])))}
    
    #Correction for dropouts (ineficient for now)
    bindex = index
    bindex[index==i]=0
    bindex[bindex!=0]=1
    DOsigfactor <- (1+cumul_zeros(bindex))
    
    #Kmat (Train)
    Z = trainData[,-seqY]
    ncolZ = ncol(Z)
    X = Z[,(ncolZ-dimX+1):ncolZ]
    MX = diag(nrow(X)) - X%*%solve(crossprod(X)+lambda2*diag(ncol(X)))%*%t(X)
    MXZ=MX%*%Z
    
    for(m in 1:nf){
        begin=(m-1)*(T-1)+1
        end=(m)*(T-1)
        if(length(sweigths)>1){MXZ[,begin:end]= (sweigths[m])*MXZ[,begin:end]}
        for(tt in 1:(T-1)){
        MXZ[,(m-1)*(T-1)+tt]= MXZ[,(m-1)*(T-1)+tt]*DOsigfactor[tt+1]
        }
    }
    
    Kmat = MXZ%*%t(MXZ)
    
    #kmat (for TEST)
    z = testData[,-seqY]
    ncolz = ncol(z)
    x = z[,(ncolz-dimX+1):ncolz]
    Mx = diag(nrow(x)) - x%*%solve(crossprod(x)+lambda2*diag(ncol(x)))%*%t(x)
    mxz=Mx%*%z
    if(length(sweigths)==1){kmat = mxz%*%t(MXZ)}
    else{
      for(m in 1:length(sweigths)){
        begin=(m-1)*(T-1)+1
        end=(m)*(T-1)
        mxz[,begin:end]= sweigths[m]*mxz[,begin:end]
      }
      kmat = mxz%*%t(MXZ)
    }
    
    #OOS pred
    for(j in 1:length(lambdavec)){
      #print(j)
      pred <- kmat%*%solve(Kmat+lambdavec[j]*Lambda_T,MX%*%trainData[,seqY])
      PMSE[i,j] <- mean((pred-Mx%*%testData[,1:dim(Y)[2]])^2)
    }
  }
  
  #Find min
  score <- colMeans(PMSE)
  lambdastarpos <- which(score == min(score))
  finalmse <- score[lambdastarpos]
  lambdastar <- lambdavec[lambdastarpos]
  
  if(length(lambdavec)>1){
  #one SE rule
  SE = array(NA,dim=c(1,length(lambdavec)))
  for(j in 1:length(lambdavec)){
    SE[1,j]=sd(PMSE[,j])/sqrt(k)
  }
  se <- SE[lambdastarpos] #mean(SE)
  scoreub = score + as.numeric(SE) 
  scorelb = score - as.numeric(SE) 
  
  lambda1sepos <- lambdastarpos
  repeat { 
    lambda1sepos <- lambda1sepos + 1
    if(lambda1sepos>=length(score)) {break}
    if(score[lambda1sepos]>(finalmse+se)) {break}
  }
  lambda1se <- lambdavec[lambda1sepos]
  
  lambda2sepos <- lambdastarpos
  repeat { 
    lambda2sepos <- lambda2sepos + 1
    if(lambda2sepos>=length(score)) {break}
    if(score[lambda2sepos]>(finalmse+2*se)) {break}
  }
  lambda2se <- lambdavec[lambda2sepos]
  
  #Plot
  if(plot==1){
    limit2 <- lambda2sepos
    repeat { 
      limit2 <- limit2 + 1
      if(limit2>=length(score)) {break}
      if(score[limit2]>(finalmse+20*se)) {break}
    }
    limit1 <- lambdastarpos
    repeat { 
      limit1 <- limit1 - 1
      if(limit1<=1) {break}
      if(score[limit1]>(finalmse+20*se)) {break}
    }
  ts.plot(score[limit1:limit2]) #,scoreub,scorelb)
  abline(h=finalmse+se, col="blue")
  abline(h=finalmse+2*se, col="purple")
  }
  return(list(minimizer=lambdastar,minima=finalmse,
              minimizer1se=lambda1se,minimizer2se=lambda2se))
  }
  else{return(list(minimizer=lambdastar,minima=finalmse))}
}


cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the 
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}
